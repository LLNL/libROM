/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Compile with: make parametric_tw_csv
//
// Generate CSV or HDF database on heat conduction with either
// heat_conduction_csv.sh or heat_conduction_hdf.sh (HDF is more efficient).
//
// =============================================================================
//
// Parametric serial DMD command (for HDF version, append -hdf):
//   parametric_tw_csv -o hc_parametric_serial -rdim 16 -dtc 0.01 -offline
//   parametric_tw_csv -o hc_parametric_serial -rdim 16 -dtc 0.01 -online
//
// Final-time prediction error (Last line in run/hc_parametric_serial/dmd_par5_prediction_error.csv):
//   0.0012598331433506
//
// Parametric time windowing DMD command (for HDF version, append -hdf):
//   parametric_tw_csv -o hc_parametric_tw -nwinsamp 25 -dtc 0.01 -offline
//   parametric_tw_csv -o hc_parametric_tw -nwinsamp 25 -dtc 0.01 -online
//
// Final-time prediction error (Last line in run/hc_parametric_tw/dmd_par5_prediction_error.csv):
//   0.0006507358659606
//
// =============================================================================
//
// Description: Parametric time windowing DMD on general CSV datasets.
//
// User specify file locations and names by -list LIST_DIR -train-set TRAIN_LIST -test-set TEST_LIST -data DATA_DIR -var VAR_NAME -o OUT_DIR
//
// File structure:
// 1. LIST_DIR/TRAIN_LIST.csv             -- each row specifies one training DATASET
// 2. LIST_DIR/TEST_LIST.csv              -- each row specifies one testing DATASET
// 3. LIST_DIR/DATASET.csv                -- each row specifies the suffix of one STATE in DATASET
// 4. DATA_DIR/dim.csv                    -- specifies the dimension of VAR_NAME
// 5. DATA_DIR/DATASET/tval.csv           -- specifies the time instances
// 6. DATA_DIR/DATASET/STATE/VAR_NAME.csv -- each row specifies one value of VAR_NAME of STATE
// 7. DATA_DIR/DATASET/TEMPORAL_IDX.csv   -- (optional) specifies the first and last temporal index in DATASET
// 8. DATA_DIR/SPATIAL_IDX.csv            -- (optional) each row specifies one spatial index of VAR_NAME
// 9. run/OUT_DIR/indicator_val.csv       -- (optional) each row specifies one indicator endpoint value

#include "mfem.hpp"
#include "algo/DMD.h"
#include "algo/AdaptiveDMD.h"
#include "algo/NonuniformDMD.h"
#include "linalg/Vector.h"
#include "linalg/Matrix.h"
#include "utils/HDFDatabase.h"
#include "utils/CSVDatabase.h"
#include <cmath>
#include <iostream>
#include <limits>

#ifndef _WIN32
#include <sys/stat.h>  // mkdir
#else
#include <direct.h>    // _mkdir
#define mkdir(dir, mode) _mkdir(dir)
#endif

using namespace std;
using namespace mfem;


string GetFilenameHDF(const char *data_dir, const char *sim_name,
                      string const& par_dir, const char *hdf_name,
                      const int myid, const int mode)
{
    switch(mode) {
    case 1:
        cout << string(data_dir) + "/" + sim_name + par_dir + "_" + to_string(
                 myid) + ".hdf" << endl;
        return string(data_dir) + "/" + sim_name + par_dir + "_"
               + to_string(myid) + ".hdf";
    default:
        cout << string(data_dir) + "/" + sim_name + "/" + string(hdf_name) + "_"
             + par_dir + "_" + to_string(myid) + ".hdf" << endl;
        return string(data_dir) + "/" + sim_name + par_dir + "/" + string(hdf_name)
               + "_" + to_string(myid) + ".hdf";
    }
}

int main(int argc, char *argv[])
{
    const int infty = numeric_limits<int>::max();
    int precision = 16;
    cout.precision(precision);

    // 1. Initialize MPI.
    MPI_Session mpi;
    int num_procs = mpi.WorldSize();
    int myid = mpi.WorldRank();

    // 2. Parse command-line options.
    bool offline = false;
    bool online = false;
    bool predict = false;
    double t_final = -1.0;
    double dtc = 0.0;
    double ddt = 0.0;
    int numWindows = 0;
    int windowNumSamples = infty;
    int windowOverlapSamples = 0;
    bool offset_indicator = false;
    const char *rbf = "G";
    const char *interp_method = "LS";
    double admd_closest_rbf_val = 0.9;
    double pdmd_closest_rbf_val = 0.9;
    double ef = 0.9999;
    int rdim = -1;
    const char *list_dir = "dmd_list";
    const char *data_dir = "dmd_data";
    const char *sim_name = "sim";
    const char *var_name = "sol";
    const char *train_list = "dmd_train_parametric";
    const char *test_list = "dmd_test";
    const char *temporal_idx_list = "temporal_idx";
    const char *spatial_idx_list = "spatial_idx";
    const char *hdf_name = "dmd";
    const char *snap_pfx = "step";
    const char *basename = "";
    bool save_csv = false;
    bool save_hdf = false;
    bool csvFormat = true;
    bool useWindowDims = false;
    int subsample = 0;
    int eval_subsample = 0;
    int fileNameMode = 0;  // Mode for HDF filenames

    OptionsParser args(argc, argv);
    args.AddOption(&offline, "-offline", "--offline", "-no-offline", "--no-offline",
                   "Enable or disable the offline phase.");
    args.AddOption(&online, "-online", "--online", "-no-online", "--no-online",
                   "Enable or disable the online phase.");
    args.AddOption(&predict, "-predict", "--predict", "-no-predict", "--no-predict",
                   "Enable or disable DMD prediction.");
    args.AddOption(&t_final, "-tf", "--t-final",
                   "Final time.");
    args.AddOption(&dtc, "-dtc", "--dtc",
                   "Fixed (constant) dt.");
    args.AddOption(&ddt, "-ddt", "--dtime-step",
                   "Desired Time step.");
    args.AddOption(&numWindows, "-nwin", "--numwindows",
                   "Number of DMD windows.");
    args.AddOption(&windowNumSamples, "-nwinsamp", "--numwindowsamples",
                   "Number of samples in DMD windows.");
    args.AddOption(&windowOverlapSamples, "-nwinover", "--numwindowoverlap",
                   "Number of samples for DMD window overlap.");
    args.AddOption(&offset_indicator, "-os", "--offset-indicator", "-no-os",
                   "--no-offset-indicator",
                   "Enable or distable the option of offset indicator.");
    args.AddOption(&rbf, "-rbf", "--radial-basis-function",
                   "Radial basis function used in interpolation. Options: \"G\", \"IQ\", \"IMQ\".");
    args.AddOption(&interp_method, "-interp", "--interpolation-method",
                   "Method of interpolation. Options: \"LS\", \"IDW\", \"LP\".");
    args.AddOption(&admd_closest_rbf_val, "-acrv", "--admd-crv",
                   "Adaptive DMD closest RBF value.");
    args.AddOption(&pdmd_closest_rbf_val, "-pcrv", "--pdmd-crv",
                   "Parametric DMD closest RBF value.");
    args.AddOption(&ef, "-ef", "--energy-fraction",
                   "Energy fraction for DMD.");
    args.AddOption(&rdim, "-rdim", "--rdim",
                   "Reduced dimension for DMD.");
    args.AddOption(&list_dir, "-list", "--list-directory",
                   "Location of training and testing data list.");
    args.AddOption(&data_dir, "-data", "--data-directory",
                   "Location of training and testing data.");
    args.AddOption(&hdf_name, "-hdffile", "--hdf-file",
                   "Base of name of HDF file for training and testing data.");
    args.AddOption(&sim_name, "-sim", "--sim-name",
                   "Name of simulation.");
    args.AddOption(&var_name, "-var", "--variable-name",
                   "Name of variable.");
    args.AddOption(&train_list, "-train-set", "--training-set-name",
                   "Name of the training datasets within the list directory.");
    args.AddOption(&test_list, "-test-set", "--testing-set-name",
                   "Name of the testing datasets within the list directory.");
    args.AddOption(&temporal_idx_list, "-t-idx", "--temporal-index",
                   "Name of the file indicating bound of temporal indices.");
    args.AddOption(&spatial_idx_list, "-x-idx", "--spatial-index",
                   "Name of the file indicating spatial indices.");
    args.AddOption(&snap_pfx, "-snap-pfx", "--snapshot-prefix",
                   "Prefix of snapshots.");
    args.AddOption(&basename, "-o", "--outputfile-name",
                   "Name of the sub-folder to dump files within the run directory.");
    args.AddOption(&save_csv, "-save", "--save", "-no-save", "--no-save",
                   "Enable or disable prediction result output (files in CSV format).");
    args.AddOption(&save_hdf, "-save-hdf", "--save-hdf", "-no-save-hdf",
                   "--no-save-hdf",
                   "Enable or disable prediction result output (files in HDF format).");
    args.AddOption(&csvFormat, "-csv", "--csv", "-hdf", "--hdf",
                   "Use CSV or HDF format for input files.");
    args.AddOption(&useWindowDims, "-wdim", "--wdim", "-no-wdim", "--no-wdim",
                   "Use DMD dimensions for each window, input from a CSV file.");
    args.AddOption(&subsample, "-subs", "--subsample",
                   "Subsampling factor for training snapshots.");
    args.AddOption(&eval_subsample, "-esubs", "--eval_subsample",
                   "Subsampling factor for evaluation.");
    args.AddOption(&fileNameMode, "-hdfmode", "--hdfmodefilename",
                   "HDF filename mode.");
    args.Parse();
    if (!args.Good())
    {
        if (mpi.Root())
        {
            args.PrintUsage(cout);
        }
        return 1;
    }
    if (mpi.Root())
    {
        args.PrintOptions(cout);
    }

    CAROM_VERIFY(!(offline && online) && (offline || online));
    CAROM_VERIFY(!(dtc > 0.0 && ddt > 0.0));
    double dt_est = max(ddt, dtc);

    if (t_final > 0.0)
    {
        save_csv = true;
    }

    string outputPath = "run";
    if (string(basename) != "") {
        outputPath += "/" + string(basename);
    }

    if (mpi.Root()) {
        const char path_delim = '/';
        string::size_type pos = 0;
        do {
            pos = outputPath.find(path_delim, pos+1);
            string subdir = outputPath.substr(0, pos);
            mkdir(subdir.c_str(), 0777);
        }
        while (pos != string::npos);
    }

    CAROM::CSVDatabase csv_db;
    CAROM::Database *db = NULL;
    CAROM::HDFDatabase *hdf_db = NULL;
    string prefix = "";
    string suffix = "";
    if (csvFormat)
    {
        db = new CAROM::CSVDatabase();
        suffix = ".csv";
    }
    else
    {
        db = new CAROM::HDFDatabase();
    }

    if (save_hdf)
    {
        hdf_db = new CAROM::HDFDatabase();
        hdf_db->create("prediction" + to_string(myid) + ".hdf");
    }

    string variable = string(var_name);
    int nelements = 0;

    if (csvFormat)
        db->getIntegerArray(string(data_dir) + "/dim.csv", &nelements, 1);
    else
    {
        db->open(GetFilenameHDF(data_dir, sim_name, "0", hdf_name,
                                myid, fileNameMode), "r");
        nelements = db->getDoubleArraySize("step0sol");
        db->close();
    }

    CAROM_VERIFY(nelements > 0);
    if (myid == 0)
    {
        cout << "Variable " << var_name << " has dimension " << nelements
             << "." << endl;
    }

    int dim = nelements;
    vector<int> idx_state;
    csv_db.getIntegerVector(string(data_dir) + "/" + string(
                                spatial_idx_list) + ".csv", idx_state, false);
    if (idx_state.size() > 0)
    {
        dim = idx_state.size();
        if (myid == 0)
        {
            cout << "Restricting on " << dim << " entries out of " << nelements
                 << "." << endl;
        }
    }

    vector<double> indicator_val;
    if (online || numWindows > 0)
    {
        csv_db.getDoubleVector(string(outputPath) + "/indicator_val.csv",
                               indicator_val, false);
        if (indicator_val.size() > 0)
        {
            if (numWindows > 0)
            {
                CAROM_VERIFY(numWindows == indicator_val.size());
            }
            else
            {
                numWindows = indicator_val.size();
            }
            if (myid == 0)
            {
                cout << "Read indicator range partition with " << numWindows
                     << " windows." << endl;
            }
        }
    }

    vector<string> training_par_list, testing_par_list; // DATASET info
    vector<string> par_dir_list; // DATASET name
    vector<CAROM::Vector*> par_vectors; // DATASET param

    vector<double> indicator_init, indicator_last; // DATASET indicator range

    csv_db.getStringVector(string(list_dir) + "/" + train_list + ".csv",
                           training_par_list, false);
    int npar = training_par_list.size();
    CAROM_VERIFY(npar > 0);
    if (myid == 0)
    {
        cout << "Loading " << npar << " training datasets." << endl;
    }

    int dpar = -1;
    for (int idx_dataset = 0; idx_dataset < npar; ++idx_dataset)
    {
        stringstream par_ss(training_par_list[idx_dataset]); // training DATASET
        vector<string> par_info;
        string par_entry;
        while (getline(par_ss, par_entry, ','))
        {
            par_info.push_back(par_entry);
        }

        dpar = par_info.size() - 1;
        CAROM::Vector* curr_par = new CAROM::Vector(dpar, false);

        if (idx_dataset == 0)
        {
            CAROM_VERIFY(dpar > 0);
            if (myid == 0)
            {
                cout << "Dimension of parametric space = " << dpar << "." << endl;
            }
        }
        else
        {
            CAROM_VERIFY(dpar == par_info.size() - 1);
        }

        string par_dir = par_info[0];
        par_dir_list.push_back(par_dir);

        if (!csvFormat)
        {
            db->open(GetFilenameHDF(data_dir, sim_name, par_dir, hdf_name,
                                    myid, fileNameMode), "r");
        }

        int snap_bound_size = 0;
        if (csvFormat)
        {
            prefix = string(data_dir) + "/" + par_dir + "/";
            snap_bound_size = csv_db.getLineCount(prefix + temporal_idx_list + suffix);
        }
        else
        {
            db->open(GetFilenameHDF(data_dir, sim_name, par_dir, hdf_name,
                                    myid, fileNameMode), "r");
            db->getInteger("snap_bound_size", snap_bound_size);
        }

        int num_snap_orig = 0;
        if (csvFormat)
        {
            num_snap_orig = csv_db.getLineCount(string(list_dir) + "/" + par_dir + ".csv");
        }
        else
        {
            db->open(GetFilenameHDF(data_dir, sim_name, par_dir, hdf_name,
                                    myid, fileNameMode), "r");
            db->getInteger("numsnap", num_snap_orig);
        }

        CAROM_VERIFY(num_snap_orig > 0);
        vector<string> snap_list(num_snap_orig);
        vector<int> snap_index_list(num_snap_orig);
        if (csvFormat)
        {
            csv_db.getStringVector(string(list_dir) + "/" + par_dir + ".csv",
                                   snap_list, false);
        }
        else
        {
            db->getIntegerArray("snap_list", snap_index_list.data(),
                                num_snap_orig);
        }

        vector<int> snap_bound(snap_bound_size);
        if (csvFormat) prefix = string(data_dir) + "/" + par_dir + "/";
        db->getIntegerArray(prefix + temporal_idx_list + suffix,
                            snap_bound.data(), snap_bound_size);
        if (snap_bound.size() > 0)
        {
            snap_bound[0] -= 1;
            snap_bound[1] -= 1;
            if (myid == 0)
            {
                cout << "Restricting on snapshot #" << snap_bound[0] << " to #"
                     << snap_bound[1] << "." << endl;
            }
        }
        else
        {
            snap_bound.push_back(0);
            snap_bound.push_back(num_snap_orig-1);
        }

        for (int par_order = 0; par_order < dpar; ++par_order)
        {
            curr_par->item(par_order) = stod(par_info[par_order+1]);
        }
        par_vectors.push_back(curr_par);

        if (offline)
        {
            vector<double> tvec(num_snap_orig);
            if (csvFormat) prefix = string(data_dir) + "/" + par_dir + "/";
            db->getDoubleArray(prefix + "tval" + suffix, tvec.data(), num_snap_orig);

            if (offset_indicator)
            {
                indicator_init.push_back(0.0);
                indicator_last.push_back(tvec[snap_bound[1]]-tvec[snap_bound[0]]);
            }
            else
            {
                indicator_init.push_back(tvec[snap_bound[0]]);
                indicator_last.push_back(tvec[snap_bound[1]]);
            }
        }

        db->close();
    }

    CAROM_VERIFY(windowOverlapSamples < windowNumSamples);
    if (offline && indicator_val.size() == 0)
    {
        double indicator_min = *min_element(indicator_init.begin(),
                                            indicator_init.end());
        double indicator_max = *max_element(indicator_last.begin(),
                                            indicator_last.end());
        numWindows = (windowNumSamples < infty) ? round((indicator_max -
                     indicator_min) / (dt_est * windowNumSamples)) : 1;
        for (int window = 0; window < numWindows; ++window)
        {
            indicator_val.push_back(indicator_min + dt_est * windowNumSamples * window);
        }
        if (myid == 0)
        {
            cout << "Created new indicator range partition with " << numWindows
                 << " windows." << endl;
            csv_db.putDoubleVector(string(outputPath) + "/indicator_val.csv",
                                   indicator_val, numWindows);
        }
    }

    CAROM_VERIFY(numWindows > 0);
    if (myid == 0)
    {
        if (numWindows > 1)
        {
            cout << "Using time windowing DMD with " << numWindows << " windows." << endl;
        }
        else
        {
            cout << "Using serial DMD." << endl;
        }
    }

    vector<int> windowDim;
    if (useWindowDims)
    {
        csv_db.getIntegerVector("run/maxDim.csv", windowDim);
        CAROM_VERIFY(windowDim.size() == numWindows);
    }

    StopWatch dmd_training_timer, dmd_preprocess_timer, dmd_prediction_timer;
    vector<vector<CAROM::DMD*>> dmd;
    vector<CAROM::DMD*> dmd_curr_par;
    double* sample = new double[dim];

    if (offline)
    {
        dmd_training_timer.Start();

        vector<int> maxDim;
        maxDim.assign(numWindows, 0);
        vector<int> minSamp;
        minSamp.assign(numWindows, -1);

        for (int idx_dataset = 0; idx_dataset < npar; ++idx_dataset)
        {
            string par_dir = par_dir_list[idx_dataset];
            if (myid == 0)
            {
                cout << "Loading samples for " << par_dir << " to train DMD." << endl;
            }

            dmd_curr_par.assign(numWindows, nullptr);
            for (int window = 0; window < numWindows; ++window)
            {
                if (ddt > 0.0)
                {
                    dmd_curr_par[window] = new CAROM::AdaptiveDMD(dim, ddt, string(rbf),
                            string(interp_method), admd_closest_rbf_val);
                }
                else if (dtc > 0.0)
                {
                    dmd_curr_par[window] = new CAROM::DMD(dim, dtc, true);
                }
                else
                {
                    dmd_curr_par[window] = new CAROM::NonuniformDMD(dim);
                }
            }
            dmd.push_back(dmd_curr_par);
            dmd_curr_par.clear();

            int num_snap_orig = 0;
            if (csvFormat)
            {
                num_snap_orig = csv_db.getLineCount(string(list_dir) + "/" + par_dir + ".csv");
            }
            else
            {
                db->open(GetFilenameHDF(data_dir, sim_name, par_dir, hdf_name,
                                        myid, fileNameMode), "r");
                db->getInteger("numsnap", num_snap_orig);
            }

            CAROM_VERIFY(num_snap_orig > 0);
            vector<string> snap_list(num_snap_orig);
            vector<int> snap_index_list(num_snap_orig);
            if (csvFormat)
            {
                csv_db.getStringVector(string(list_dir) + "/" + par_dir + ".csv", snap_list,
                                       false);
            }
            else
            {
                db->getIntegerArray("snap_list", snap_index_list.data(),
                                    num_snap_orig);
            }

            vector<double> tvec(num_snap_orig);
            db->getDoubleArray(prefix + "tval" + suffix, tvec.data(), num_snap_orig);

            int snap_bound_size = 0;
            if (csvFormat)
            {
                prefix = string(data_dir) + "/" + par_dir + "/";
                snap_bound_size = csv_db.getLineCount(prefix + temporal_idx_list + suffix);
            }
            else
            {
                db->open(GetFilenameHDF(data_dir, sim_name, par_dir, hdf_name,
                                        myid, fileNameMode), "r");
                db->getInteger("snap_bound_size", snap_bound_size);
            }

            vector<int> snap_bound(snap_bound_size);
            db->getIntegerArray(prefix + temporal_idx_list + suffix,
                                snap_bound.data(), snap_bound_size);

            if (snap_bound.size() > 0)
            {
                snap_bound[0] -= 1;
                snap_bound[1] -= 1;
                CAROM_VERIFY(snap_bound.size() == 2);
                if (myid == 0)
                {
                    cout << "Restricting on snapshot #" << snap_bound[0]
                         << " to #" << snap_bound[1] << "." << endl;
                }
            }
            else
            {
                snap_bound.push_back(0);
                snap_bound.push_back(num_snap_orig-1);
            }

            int curr_window = 0;
            int overlap_count = 0;
            for (int idx_snap = snap_bound[0]; idx_snap <= snap_bound[1]; ++idx_snap)
            {
                if (subsample > 1 && (idx_snap % subsample != 0)
                        && idx_snap > snap_bound[0] && idx_snap < snap_bound[1])
                    continue;

                string snap = snap_pfx;
                double tval = tvec[idx_snap];
                if (csvFormat)
                {
                    snap += snap_list[idx_snap]; // STATE
                    string data_filename = string(data_dir) + "/" + par_dir + "/" + snap +
                                           "/" + variable + ".csv"; // path to VAR_NAME.csv
                    db->getDoubleArray(data_filename, sample, nelements, idx_state);
                }
                else
                {
                    snap += to_string(snap_index_list[idx_snap]); // STATE
                    db->getDoubleArray(snap + "sol", sample, nelements, idx_state);
                }

                dmd[idx_dataset][curr_window]->takeSample(sample,
                        tval - offset_indicator * tvec[snap_bound[0]]);
                if (overlap_count > 0)
                {
                    dmd[idx_dataset][curr_window-1]->takeSample(sample,
                            tval - offset_indicator * tvec[snap_bound[0]]);
                    overlap_count -= 1;
                }
                // a rough estimate to correct the precision of the indicator range partition
                double indicator_snap = tval - offset_indicator * tvec[snap_bound[0]]
                                        + dt_est / 100.0;

                if (curr_window+1 < numWindows && idx_snap+1 <= snap_bound[1]
                        && indicator_snap > indicator_val[curr_window+1])
                {
                    overlap_count = windowOverlapSamples;
                    curr_window += 1;
                    dmd[idx_dataset][curr_window]->takeSample(sample,
                            tval - offset_indicator * tvec[snap_bound[0]]);
                }
            }

            if (myid == 0)
            {
                cout << "Loaded " << snap_bound[1] - snap_bound[0] + 1
                     << " samples for " << par_dir << "." << endl;
            }

            for (int window = 0; window < numWindows; ++window)
            {
                if (useWindowDims)
                    rdim = windowDim[window];

                if (rdim != -1)
                {
                    if (myid == 0)
                    {
                        cout << "Creating DMD model #" << window << " with rdim: " << rdim << endl;
                    }
                    dmd[idx_dataset][window]->train(rdim);
                }
                else if (ef > 0.0)
                {
                    if (myid == 0)
                    {
                        cout << "Creating DMD model #" << window
                             << " with energy fraction: " << ef << endl;
                    }
                    dmd[idx_dataset][window]->train(ef);
                }
                if (window > 0 && predict)
                {
                    if (myid == 0)
                    {
                        cout << "Projecting initial condition at t = "
                             << indicator_val[window] + offset_indicator * tvec[snap_bound[0]]
                             << " for DMD model #" << window << endl;
                    }
                    CAROM::Vector* init_cond = dmd[idx_dataset][window-1]->predict(
                                                   indicator_val[window]);
                    dmd[idx_dataset][window]->projectInitialCondition(init_cond);
                    delete init_cond;
                }

                // Make a directory for this window, only on the first parameter.
                if (idx_dataset == 0)
                {
                    string outWindow = outputPath + "/window" + to_string(window);
                    mkdir(outWindow.c_str(), 0777);
                }

                dmd[idx_dataset][window]->save(outputPath + "/window"
                                               + to_string(window) + "/par"
                                               + to_string(idx_dataset));

                if (myid == 0)
                {
                    dmd[idx_dataset][window]->summary(outputPath + "/window"
                                                      + to_string(window) + "/par"
                                                      + to_string(idx_dataset));

                    cout << "Window " << window << ", DMD " << idx_dataset << " dim "
                         << dmd[idx_dataset][window]->getDimension() << endl;
                }

                const int dim_w = std::min(dmd[idx_dataset][window]->getDimension(),
                                           dmd[idx_dataset][window]->getNumSamples()-1);
                maxDim[window] = std::max(maxDim[window], dim_w);

                if (minSamp[window] < 0 ||
                        dmd[idx_dataset][window]->getNumSamples() < minSamp[window])
                    minSamp[window] = dmd[idx_dataset][window]->getNumSamples();
            } // escape for-loop over window

            db->close();

            if (!online && !predict)
            {
                for (int window = 0; window < numWindows; ++window)
                {
                    delete dmd[idx_dataset][window];
                    dmd[idx_dataset][window] = nullptr;
                }
            }
        } // escape for-loop over idx_dataset
        dmd_training_timer.Stop();

        // Limit maxDim by minSamp-1
        for (int window = 0; window < numWindows; ++window)
            maxDim[window] = std::min(maxDim[window], minSamp[window]-1);

        // Write maxDim to a CSV file
        if (myid == 0) csv_db.putIntegerArray("run/maxDim.csv", maxDim.data(),
                                                  numWindows);
    } // escape if-statement of offline

    CAROM::Vector* curr_par = new CAROM::Vector(dpar, false);

    if (online)
    {
        par_dir_list.clear();

        dmd_preprocess_timer.Start();
        csv_db.getStringVector(string(list_dir) + "/" + test_list + ".csv",
                               testing_par_list, false);
        npar = testing_par_list.size();
        CAROM_VERIFY(npar > 0);

        if (myid == 0)
        {
            cout << "Loading " << npar << " testing datasets." << endl;
        }

        dmd_curr_par.assign(numWindows, nullptr);
        dmd.assign(numWindows, dmd_curr_par);

        int num_tests = 0;
        vector<double> prediction_time, prediction_error;

        for (int idx_dataset = 0; idx_dataset < npar; ++idx_dataset)
        {
            stringstream par_ss(testing_par_list[idx_dataset]); // testing DATASET
            vector<string> par_info;
            string par_entry;
            while (getline(par_ss, par_entry,','))
            {
                par_info.push_back(par_entry);
            }

            CAROM_VERIFY(dpar == par_info.size() - 1);

            string par_dir = par_info[0];
            par_dir_list.push_back(par_dir);
            if (myid == 0)
            {
                cout << "Interpolating DMD models for dataset " << par_dir << endl;
            }

            int num_snap_orig = 0;
            if (csvFormat)
            {
                num_snap_orig = csv_db.getLineCount(string(list_dir) + "/" + par_dir + ".csv");
            }
            else
            {
                db->open(GetFilenameHDF(data_dir, sim_name, par_dir, hdf_name,
                                        myid, fileNameMode), "r");
                db->getInteger("numsnap", num_snap_orig);
            }

            CAROM_VERIFY(num_snap_orig > 0);
            vector<string> snap_list(num_snap_orig);
            vector<int> snap_index_list(num_snap_orig);
            if (csvFormat)
            {
                csv_db.getStringVector(string(list_dir) + "/" + par_dir + ".csv", snap_list,
                                       false);
            }
            else
            {
                db->getIntegerArray("snap_list", snap_index_list.data(),
                                    num_snap_orig);
            }

            for (int par_order = 0; par_order < dpar; ++par_order)
            {
                curr_par->item(par_order) = stod(par_info[par_order+1]);
            }

            vector<double> tvec(num_snap_orig);
            db->getDoubleArray(prefix + "tval" + suffix, tvec.data(), num_snap_orig);

            int snap_bound_size = 0;
            if (csvFormat)
            {
                prefix = string(data_dir) + "/" + par_dir + "/";
                snap_bound_size = csv_db.getLineCount(prefix + temporal_idx_list + suffix);
            }
            else
            {
                db->open(GetFilenameHDF(data_dir, sim_name, par_dir, hdf_name,
                                        myid, fileNameMode), "r");
                db->getInteger("snap_bound_size", snap_bound_size);
            }

            vector<int> snap_bound(snap_bound_size);
            db->getIntegerArray(prefix + temporal_idx_list + suffix,
                                snap_bound.data(), snap_bound_size);

            if (snap_bound.size() > 0)
            {
                snap_bound[0] -= 1;
                snap_bound[1] -= 1;
                CAROM_VERIFY(snap_bound.size() == 2);
                if (myid == 0)
                {
                    cout << "Restricting on snapshot #" << snap_bound[0]
                         << " to " << snap_bound[1] << "." << endl;
                }
            }
            else
            {
                snap_bound.push_back(0);
                snap_bound.push_back(num_snap_orig-1);
            }

            string snap = snap_pfx;
            if (csvFormat)
            {
                snap += snap_list[snap_bound[0]]; // STATE
                string data_filename = string(data_dir) + "/" + par_dir + "/" + snap +
                                       "/" + variable + ".csv"; // path to VAR_NAME.csv
                db->getDoubleArray(data_filename, sample, nelements, idx_state);
            }
            else
            {
                snap += to_string(snap_index_list[snap_bound[0]]); // STATE
                db->getDoubleArray(snap + "sol", sample, nelements, idx_state);
            }

            for (int window = 0; window < numWindows; ++window)
            {
                std::vector<std::string> dmd_paths;
                for (int idx_trainset = 0; idx_trainset < par_vectors.size(); ++idx_trainset)
                {
                    dmd_paths.push_back(outputPath + "/window" + to_string(window)
                                        + "/par" + to_string(idx_trainset));
                }

                if (par_vectors.size() > 1)
                {
                    if (myid == 0)
                    {
                        cout << "Interpolating DMD model #" << window << endl;
                    }

                    CAROM::getParametricDMD(dmd[idx_dataset][window], par_vectors,
                                            dmd_paths, curr_par, string(rbf),
                                            string(interp_method), pdmd_closest_rbf_val);
                }
                else if (par_vectors.size() == 1 && dmd_paths.size() == 1)
                {
                    if (myid == 0)
                    {
                        cout << "Loading local DMD model #" << window << endl;
                    }
                    dmd[idx_dataset][window] = new CAROM::DMD(dmd_paths[0]);
                }

                if (myid == 0)
                {
                    cout << "Projecting initial condition at t = "
                         << indicator_val[window] + offset_indicator * tvec[snap_bound[0]]
                         << " for DMD model #" << window << endl;
                }
                CAROM::Vector* init_cond = nullptr;
                if (window == 0)
                {
                    init_cond = new CAROM::Vector(dim, true);
                    for (int i = 0; i < dim; ++i)
                    {
                        init_cond->item(i) = sample[i];
                    }
                }
                else
                {
                    init_cond = dmd[idx_dataset][window-1]->predict(indicator_val[window]);
                }
                dmd[idx_dataset][window]->projectInitialCondition(init_cond);

                const double norm_init_cond = init_cond->norm();
                if (myid == 0)
                {
                    cout << "Initial condition norm " << norm_init_cond
                         << " for parameter " << idx_dataset << ", window "
                         << window - 1 << endl;
                }

                delete init_cond;

                if (window > 0 && indicator_val[window] < t_final)
                {
                    // To save memory, delete dmd[idx_dataset][window] for windows
                    // not containing t_final.
                    if (myid == 0)
                    {
                        cout << "Deleting DMD for parameter " << idx_dataset << ", window "
                             << window - 1 << endl;
                    }

                    delete dmd[idx_dataset][window-1];
                    dmd[idx_dataset][window-1] = nullptr;
                }

            } // escape for-loop over window

            db->close();
        } // escape for-loop over idx_dataset
        dmd_preprocess_timer.Stop();
    } // escape if-statement of online

    if (online || predict)
    {
        int num_tests = 0;
        vector<double> prediction_time, prediction_error;

        for (int idx_dataset = 0; idx_dataset < npar; ++idx_dataset)
        {
            string par_dir = par_dir_list[idx_dataset];
            if (myid == 0)
            {
                cout << "Predicting solution for " << par_dir << " using DMD." << endl;
            }
            int num_snap_orig = 0;
            if (csvFormat)
            {
                num_snap_orig = csv_db.getLineCount(string(list_dir) + "/" + par_dir + ".csv");
            }
            else
            {
                db->open(GetFilenameHDF(data_dir, sim_name, par_dir, hdf_name,
                                        myid, fileNameMode), "r");
                db->getInteger("numsnap", num_snap_orig);
            }
            CAROM_VERIFY(num_snap_orig > 0);
            vector<string> snap_list(num_snap_orig);
            vector<int> snap_index_list(num_snap_orig);
            if (csvFormat)
            {
                csv_db.getStringVector(string(list_dir) + "/" + par_dir + ".csv", snap_list,
                                       false);
            }
            else
            {
                db->getIntegerArray("snap_list", snap_index_list.data(),
                                    num_snap_orig);
            }

            vector<double> tvec(num_snap_orig);
            db->getDoubleArray(prefix + "tval" + suffix, tvec.data(), num_snap_orig);

            int snap_bound_size = 0;
            if (csvFormat)
            {
                prefix = string(data_dir) + "/" + par_dir + "/";
                snap_bound_size = csv_db.getLineCount(prefix + temporal_idx_list + suffix);
            }
            else
            {
                db->open(GetFilenameHDF(data_dir, sim_name, par_dir, hdf_name,
                                        myid, fileNameMode), "r");
                db->getInteger("snap_bound_size", snap_bound_size);
            }

            vector<int> snap_bound(snap_bound_size);
            db->getIntegerArray(prefix + temporal_idx_list + suffix,
                                snap_bound.data(), snap_bound_size);

            if (snap_bound.size() > 0)
            {
                snap_bound[0] -= 1;
                snap_bound[1] -= 1;
                CAROM_VERIFY(snap_bound.size() == 2);
                if (myid == 0)
                {
                    cout << "Restricting on snapshot #" << snap_bound[0]
                         << " to " << snap_bound[1] << "." << endl;
                }
            }
            else
            {
                snap_bound.push_back(0);
                snap_bound.push_back(num_snap_orig-1);
            }

            int num_snap = snap_bound[1] - snap_bound[0] + 1;
            int curr_window = 0;
            for (int idx_snap = snap_bound[0]; idx_snap <= snap_bound[1]; ++idx_snap)
            {
                if (eval_subsample > 1 && (idx_snap % eval_subsample != 0)
                        && idx_snap > snap_bound[0] && idx_snap < snap_bound[1])
                    continue;

                string snap = snap_pfx;
                double tval = tvec[idx_snap];
                if (csvFormat)
                {
                    snap += snap_list[idx_snap]; // STATE
                    string data_filename = string(data_dir) + "/" + par_dir + "/" + snap +
                                           "/" + variable + ".csv"; // path to VAR_NAME.csv
                    db->getDoubleArray(data_filename, sample, nelements, idx_state);
                }
                else
                {
                    snap += to_string(snap_index_list[idx_snap]); // STATE
                    db->getDoubleArray(snap + "sol", sample, nelements, idx_state);
                }

                if (myid == 0)
                {
                    cout << "State " << snap << " read." << endl;
                }

                if (t_final > 0.0) // Actual prediction without true solution for comparison
                {
                    num_tests += 1;
                    while (curr_window+1 < numWindows &&
                            t_final - offset_indicator * tvec[snap_bound[0]] >
                            indicator_val[curr_window+1])
                    {
                        curr_window += 1;
                    }
                    if (myid == 0)
                    {
                        cout << "Predicting DMD solution at t = " << t_final
                             << " using DMD model #" << curr_window << endl;
                    }

                    dmd_prediction_timer.Start();
                    CAROM::Vector* result = dmd[idx_dataset][curr_window]->predict(
                                                t_final - offset_indicator * tvec[snap_bound[0]]);
                    dmd_prediction_timer.Stop();

                    if (myid == 0)
                    {
                        csv_db.putDoubleArray(outputPath + "/" + par_dir +
                                              "_final_time_prediction.csv",
                                              result->getData(), dim);
                    }

                    if (save_hdf)
                    {
                        hdf_db->putDoubleArray(snap, result->getData(), dim);
                    }

                    idx_snap = snap_bound[1]+1; // escape for-loop over idx_snap
                    delete result;
                }
                else // Verify DMD prediction results against dataset
                {
                    while (curr_window+1 < numWindows &&
                            tval - offset_indicator * tvec[snap_bound[0]] >
                            indicator_val[curr_window+1])
                    {
                        curr_window += 1;
                    }
                    if (myid == 0)
                    {
                        cout << "Predicting DMD solution #" << idx_snap
                             << " at t = " << tval << " using DMD model #"
                             << curr_window << endl;
                    }

                    dmd_prediction_timer.Start();
                    CAROM::Vector* result = dmd[idx_dataset][curr_window]->predict(
                                                tval - offset_indicator * tvec[snap_bound[0]]);
                    dmd_prediction_timer.Stop();

                    // Calculate the relative error between the DMD final solution and the true solution.
                    Vector dmd_solution(result->getData(), result->dim());
                    Vector true_solution(sample, dim);
                    Vector diff(true_solution.Size());
                    subtract(dmd_solution, true_solution, diff);

                    double tot_diff_norm = sqrt(InnerProduct(MPI_COMM_WORLD, diff, diff));
                    double tot_true_solution_norm = sqrt(InnerProduct(MPI_COMM_WORLD, true_solution,
                                                         true_solution));
                    double rel_error = tot_diff_norm / tot_true_solution_norm;

                    prediction_time.push_back(tval);
                    prediction_error.push_back(rel_error);

                    if (myid == 0)
                    {
                        cout << "Norm of true solution at t = " << tval
                             << " is " << tot_true_solution_norm << endl;
                        cout << "Absolute error of DMD solution at t = " << tval
                             << " is " << tot_diff_norm << endl;
                        cout << "Relative error of DMD solution at t = " << tval
                             << " is " << rel_error << endl;

                        if (idx_snap == snap_bound[1])  // Final step
                        {
                            std::ofstream ofs ("log.txt", std::ofstream::app);
                            ofs << windowNumSamples << "," << subsample << "," << rel_error << endl;
                            ofs.close();
                        }

                        if (save_csv)
                        {
                            csv_db.putDoubleArray(outputPath + "/" + par_dir + "_" + snap +
                                                  "_prediction.csv", result->getData(), dim);
                            if (dim < nelements)
                            {
                                csv_db.putDoubleArray(outputPath + "/" + par_dir + "_"
                                                      + snap + "_state.csv", sample, dim);
                            }
                        }
                    }

                    if (save_hdf)
                    {
                        hdf_db->putDoubleArray(snap, result->getData(), dim);
                    }

                    delete result;
                }
            }
            if (myid == 0 && t_final <= 0.0)
            {
                csv_db.putDoubleVector(outputPath + "/" + par_dir + "_prediction_time.csv",
                                       prediction_time, num_snap);
                csv_db.putDoubleVector(outputPath + "/" + par_dir + "_prediction_error.csv",
                                       prediction_error, num_snap);
            }
            prediction_time.clear();
            prediction_error.clear();
            num_tests = (t_final > 0.0) ? num_tests + 1 : num_tests + num_snap;

            db->close();
        }

        CAROM_VERIFY(num_tests > 0);

        if (myid == 0)
        {
            printf("Elapsed time for training DMD: %e second\n",
                   dmd_training_timer.RealTime());
            printf("Elapsed time for preprocessing DMD: %e second\n",
                   dmd_preprocess_timer.RealTime());
            printf("Total elapsed time for predicting DMD: %e second\n",
                   dmd_prediction_timer.RealTime());
            printf("Average elapsed time for predicting DMD: %e second\n",
                   dmd_prediction_timer.RealTime() / num_tests);
        }
    } // escape case (online || predict)

    if (save_hdf)
    {
        hdf_db->close();
        delete hdf_db;
    }

    delete[] sample;
    delete curr_par;
    for (int idx_dataset = 0; idx_dataset < npar; ++idx_dataset)
    {
        for (int window = 0; window < numWindows; ++window)
        {
            delete dmd[idx_dataset][window];
        }
    }

    return 0;
}
