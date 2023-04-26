/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Compile with: make local_dw_csv
//
// Generate CSV database on heat conduction with: heat_conduction_csv.sh
//
// =================================================================================
//
// Description: Local solution manifold decomposition DMD on general CSV datasets.
//
// User specify file locations and names by -list LIST_DIR -train-set TRAIN_LIST -test-set TEST_LIST -data DATA_DIR -var VAR_NAME -o OUT_DIR
//
// File structure:
//  1. LIST_DIR/TRAIN_LIST.csv             -- each row specifies one training DATASET
//  2. LIST_DIR/TEST_LIST.csv              -- each row specifies one testing DATASET
//  3. LIST_DIR/DATASET.csv                -- each row specifies one STATE in DATASET
//  4. DATA_DIR/dim.csv                    -- specifies the dimension of VAR_NAME
//  5. DATA_DIR/DATASET/tval.csv           -- specifies the time instances
//  6. DATA_DIR/DATASET/STATE/VAR_NAME.csv -- each row specifies one value of VAR_NAME of STATE
//  7. DATA_DIR/SPATIAL_IDX.csv            -- (optional) each row specifies one spatial index of VAR_NAME
//  8. DATA_DIR/DATASET/TEMPORAL_IDX.csv   -- (optional) specifies the first and last temporal index in DATASET
//  9. run/OUT_DIR/INDICATOR_IDX.csv       -- (optional) each row specifies one indicator index
// 10. run/OUT_DIR/INDICATOR_VAL.csv       -- (optional) each row specifies one indicator endpoint value

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
    bool train = true;
    double t_final = -1.0;
    double dtc = 0.0;
    double ddt = 0.0;
    int numWindows = 1;
    int windowNumSamples = infty;
    int windowOverlapSamples = 0;
    const char *rbf = "G";
    const char *interp_method = "LS";
    double admd_closest_rbf_val = 0.9;
    double ef = 0.9999;
    int rdim = -1;
    const char *list_dir = "hc_list";
    const char *data_dir = "hc_data";
    const char *var_name = "sol";
    const char *train_list = "hc_train_local";
    const char *test_list = "hc_test";
    const char *temporal_idx_list = "temporal_idx";
    const char *spatial_idx_list = "spatial_idx";
    const char *indicator_idx_list = "";
    const char *indicator_val_list = "indicator_val";
    const char *window_endpoint_option = "right";
    const char *basename = "";
    bool save_csv = false;

    OptionsParser args(argc, argv);
    args.AddOption(&train, "-train", "--train", "-no-train", "--no-train",
                   "Enable or disable DMD training.");
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
    args.AddOption(&rbf, "-rbf", "--radial-basis-function",
                   "Radial basis function used in interpolation. Options: \"G\", \"IQ\", \"IMQ\".");
    args.AddOption(&interp_method, "-interp", "--interpolation-method",
                   "Method of interpolation. Options: \"LS\", \"IDW\", \"LP\".");
    args.AddOption(&admd_closest_rbf_val, "-acrv", "--admd-crv",
                   "Adaptive DMD closest RBF value.");
    args.AddOption(&ef, "-ef", "--energy-fraction",
                   "Energy fraction for DMD.");
    args.AddOption(&rdim, "-rdim", "--rdim",
                   "Reduced dimension for DMD.");
    args.AddOption(&list_dir, "-list", "--list-directory",
                   "Location of training and testing data list.");
    args.AddOption(&data_dir, "-data", "--data-directory",
                   "Location of training and testing data.");
    args.AddOption(&var_name, "-var", "--variable-name",
                   "Name of variable.");
    args.AddOption(&train_list, "-train-set", "--training-set-name",
                   "Name of the training datasets within the list directory.");
    args.AddOption(&test_list, "-test-set", "--testing-set-name",
                   "Name of the testing datasets within the list directory.");
    args.AddOption(&temporal_idx_list, "-t-idx", "--temporal-index",
                   "Name of the file specifying bound of temporal indices.");
    args.AddOption(&spatial_idx_list, "-x-idx", "--spatial-index",
                   "Name of the file specifying spatial indices.");
    args.AddOption(&indicator_idx_list, "-idc-idx", "--indicator-index",
                   "Name of the file specifying indicator indices.");
    args.AddOption(&indicator_val_list, "-idc-val", "--indicator-value",
                   "Name of the file specifying indicator values.");
    args.AddOption(&window_endpoint_option, "-wep", "--wep",
                   "Determine window endpoint using previous state, current state or bisection.");
    args.AddOption(&basename, "-o", "--outputfile-name",
                   "Name of the sub-folder to dump files within the run directory.");
    args.AddOption(&save_csv, "-csv", "--csv", "-no-csv", "--no-csv",
                   "Enable or disable prediction result output (files in CSV format).");
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

    CAROM_VERIFY(!(dtc > 0.0 && ddt > 0.0));
    CAROM_VERIFY(numWindows == 1 || windowNumSamples == infty);
    CAROM_VERIFY(windowOverlapSamples < windowNumSamples);

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

    string variable = string(var_name);
    int nelements = 0;
    csv_db.getIntegerArray(string(data_dir) + "/dim.csv", &nelements, 1);
    CAROM_VERIFY(nelements > 0);
    if (myid == 0)
    {
        cout << "Variable " << var_name << " has dimension " << nelements << "." <<
             endl;
    }

    int dim = nelements;
    vector<int> idx_state;
    csv_db.getIntegerVector(string(data_dir) + "/" +
                            string(spatial_idx_list) + ".csv", idx_state, false);
    if (idx_state.size() > 0)
    {
        dim = idx_state.size();
        if (myid == 0)
        {
            cout << "Restricting on " << dim << " entries out of " << nelements << "." <<
                 endl;
        }
    }

    vector<int> indicator_idx;
    csv_db.getIntegerVector(string(outputPath) + "/" +
                            string(indicator_idx_list) + ".csv", indicator_idx, false);
    if (indicator_idx.size() > 0)
    {
        if (indicator_idx[0] >= 0)
        {
            if (myid == 0)
            {
                cout << "Using " << indicator_idx.size() << " entries out of " << dim
                     << " as indicator." << endl;
            }
        }
        else
        {
            indicator_idx.clear();
        }
    }

    vector<double> indicator_val;
    if (!train || numWindows > 1)
    {
        csv_db.getDoubleVector(string(outputPath) + "/" +
                               string(indicator_val_list) + ".csv", indicator_val, false);
        if (indicator_val.size() > 0)
        {
            CAROM_VERIFY(windowNumSamples == infty);

            if (numWindows > 1)
            {
                CAROM_VERIFY(numWindows == indicator_val.size() - 1);
            }
            else
            {
                numWindows = indicator_val.size() - 1;
            }

            if (indicator_idx.size() == 1)
            {
                indicator_idx.resize(numWindows, indicator_idx[0]);
            }
            CAROM_VERIFY(indicator_idx.size() == indicator_val.size());

            if (myid == 0)
            {
                cout << "Read indicator range partition with " << numWindows << " windows." <<
                     endl;
            }
        }
    }

    int npar = 0;
    int num_train_snap;
    vector<string> training_par_list;
    if (train)
    {
        csv_db.getStringVector(string(list_dir) + "/" + train_list + ".csv",
                               training_par_list, false);
        npar = training_par_list.size();
        CAROM_VERIFY(npar == 1);

        stringstream par_ss(training_par_list[0]); // training DATASET
        string par_dir;
        getline(par_ss, par_dir, ',');

        vector<int> snap_bound;
        csv_db.getIntegerVector(string(data_dir) + "/" + par_dir + "/" +
                                temporal_idx_list + ".csv", snap_bound, false);
        if (snap_bound.size() > 0)
        {
            snap_bound[0] -= 1;
            snap_bound[1] -= 1;
            CAROM_VERIFY(snap_bound.size() == 2);
            num_train_snap = snap_bound[1] - snap_bound[0] + 1;
        }
        else
        {
            num_train_snap = csv_db.getLineCount(string(list_dir) + "/" + par_dir + ".csv");
        }

        if (windowNumSamples < infty)
        {
            numWindows = ceil((((double) num_train_snap) - 1.5) / (double)
                              windowNumSamples);
        }
    }

    if (myid == 0)
    {
        if (numWindows > 1)
        {
            cout << "Using windowed DMD with " << numWindows << " windows." << endl;
        }
        else
        {
            cout << "Using serial DMD." << endl;
        }
    }

    vector<CAROM::DMD*> dmd;
    dmd.assign(numWindows, nullptr);
    if (train)
    {
        if (ddt > 0.0)
        {
            dmd[0] = new CAROM::AdaptiveDMD(dim, ddt, string(rbf),
                                            string(interp_method), admd_closest_rbf_val);
        }
        else if (dtc > 0.0)
        {
            dmd[0] = new CAROM::DMD(dim, dtc);
        }
        else
        {
            dmd[0] = new CAROM::NonuniformDMD(dim);
        }
    }
    else
    {
        for (int window = 0; window < numWindows; ++window)
        {
            if (myid == 0)
            {
                cout << "Loading DMD model #" << window << "." << endl;
            }
            dmd[window] = new CAROM::DMD(outputPath + "/window" + to_string(window));
        }
    }

    StopWatch dmd_training_timer, dmd_preprocess_timer, dmd_prediction_timer;
    double* sample = new double[dim];

    if (train)
    {
        dmd_training_timer.Start();

        stringstream par_ss(training_par_list[0]); // training DATASET
        string par_dir;
        getline(par_ss, par_dir, ',');

        vector<string> snap_list;
        csv_db.getStringVector(string(list_dir) + "/" + par_dir + ".csv", snap_list,
                               false);
        CAROM_VERIFY(snap_list.size() > 0);

        if (myid == 0)
        {
            cout << "Loading samples for " << par_dir << " to train DMD." << endl;
        }

        vector<double> tvec;
        csv_db.getDoubleVector(string(data_dir) + "/" + par_dir + "/tval.csv", tvec,
                               false);
        CAROM_VERIFY(tvec.size() == snap_list.size());

        vector<int> snap_bound;
        csv_db.getIntegerVector(string(data_dir) + "/" + par_dir + "/" +
                                temporal_idx_list + ".csv", snap_bound, false);
        if (snap_bound.size() > 0)
        {
            snap_bound[0] -= 1;
            snap_bound[1] -= 1;
            CAROM_VERIFY(snap_bound.size() == 2);
            if (myid == 0)
            {
                cout << "Restricting on snapshot #" << snap_bound[0] << " to " << snap_bound[1]
                     << "." << endl;
            }
        }
        else
        {
            snap_bound.push_back(0);
            snap_bound.push_back(snap_list.size()-1);
        }

        int curr_window = 0;
        int overlap_count = 0;
        int min_idx_snap = -1;
        int max_idx_snap = -1;
        double curr_indicator_val = -1.0;

        for (int idx_snap = snap_bound[0]; idx_snap <= snap_bound[1]; ++idx_snap)
        {
            string snap = snap_list[idx_snap]; // STATE
            double tval = tvec[idx_snap];

            string data_filename = string(data_dir) + "/" + par_dir + "/" + snap + "/" +
                                   variable + ".csv"; // path to VAR_NAME.csv
            csv_db.getDoubleArray(data_filename, sample, nelements, idx_state);

            if (myid == 0)
            {
                cout << "State #" << idx_snap << " - " << data_filename
                     << " is read." << endl;
            }

            if (min_idx_snap == -1)
            {
                curr_indicator_val = (indicator_idx.size() == 0) ? tval :
                                     sample[indicator_idx[0]];
                if (indicator_val.size() > 0)
                {
                    if (curr_indicator_val >= indicator_val[0])
                    {
                        min_idx_snap = idx_snap;
                        indicator_val[0] = curr_indicator_val;
                        if (myid == 0)
                        {
                            cout << "State #" << idx_snap << " - " << data_filename
                                 << " is the beginning of window 0." << endl;
                        }
                    }
                    else
                    {
                        if (myid == 0)
                        {
                            cout << "State #" << idx_snap << " - " << data_filename
                                 << " is omitted." << endl;
                        }
                        continue;
                    }
                }
                else
                {
                    min_idx_snap = snap_bound[0];
                    max_idx_snap = snap_bound[1];
                    indicator_val.push_back(curr_indicator_val);
                }
            }

            dmd[curr_window]->takeSample(sample, tval);
            if (overlap_count > 0)
            {
                dmd[curr_window-1]->takeSample(sample, tval);
                overlap_count -= 1;
            }

            if (curr_window+1 < numWindows && idx_snap+1 <= snap_bound[1])
            {
                curr_indicator_val = (indicator_idx.size() == 0) ? tval :
                                     sample[indicator_idx[curr_window+1]];
                bool new_window = false;
                if (windowNumSamples < infty)
                {
                    new_window = (idx_snap >= min_idx_snap + (curr_window+1)*windowNumSamples);
                    if (new_window)
                    {
                        indicator_val.push_back(curr_indicator_val);
                    }
                }
                else
                {
                    new_window = (curr_indicator_val >= indicator_val[curr_window+1]);
                    if (new_window)
                    {
                        indicator_val[curr_window+1] = curr_indicator_val;
                    }
                }
                if (new_window)
                {
                    if (myid == 0)
                    {
                        cout << "State #" << idx_snap << " - " << data_filename
                             << " is the beginning of window " << curr_window+1 << "." << endl;
                    }

                    curr_window += 1;
                    overlap_count = windowOverlapSamples;

                    if (ddt > 0.0)
                    {
                        dmd[curr_window] = new CAROM::AdaptiveDMD(dim, ddt, string(rbf),
                                string(interp_method), admd_closest_rbf_val);
                    }
                    else if (dtc > 0.0)
                    {
                        dmd[curr_window] = new CAROM::DMD(dim, dtc);
                    }
                    else
                    {
                        dmd[curr_window] = new CAROM::NonuniformDMD(dim);
                    }
                    dmd[curr_window]->takeSample(sample, tval);
                }
            }

            if (max_idx_snap == -1 && curr_window == numWindows-1)
            {
                curr_indicator_val = (indicator_idx.size() == 0) ? tval :
                                     sample[indicator_idx[numWindows]];
                if (curr_indicator_val >= indicator_val[numWindows])
                {
                    max_idx_snap = idx_snap;
                    indicator_val[numWindows] = curr_indicator_val;
                    if (myid == 0)
                    {
                        cout << "State #" << idx_snap << " - " << data_filename
                             << " is the end of window " << numWindows-1 << "." << endl;
                    }
                    break;
                }
            }
        }

        if (indicator_idx.size() > 0)
        {
            numWindows = curr_window + 1;
            dmd.resize(numWindows);
        }
        CAROM_VERIFY(numWindows == curr_window + 1);

        if (myid == 0)
        {
            cout << "Sampling on " << par_dir << " completed" << endl;
            cout << "Samples loaded " << num_train_snap << endl;
            cout << "Samples used " << max_idx_snap - min_idx_snap + 1 << endl;;
            cout << "Number of windows: " << numWindows << endl;
            if (windowNumSamples < infty)
            {
                csv_db.putDoubleVector(string(outputPath) + "/" + string(
                                           indicator_val_list) + ".csv",
                                       indicator_val, numWindows);
            }
        }

        for (int window = 0; window < numWindows; ++window)
        {
            if (rdim != -1)
            {
                if (myid == 0)
                {
                    cout << "Creating DMD model #" << window << " with rdim: " << rdim << endl;
                }
                dmd[window]->train(rdim);
            }
            else if (ef != -1)
            {
                if (myid == 0)
                {
                    cout << "Creating DMD model #" << window << " with energy fraction: " << ef <<
                         endl;
                }
                dmd[window]->train(ef);
            }
            dmd[window]->save(outputPath + "/window" + to_string(window));
            if (myid == 0)
            {
                dmd[window]->summary(outputPath + "/window" + to_string(window));
            }
        }

        dmd_training_timer.Stop();

        if (ddt > 0.0)
        {
            CAROM::AdaptiveDMD* admd = nullptr;
            CAROM::Vector interp_snap(dim, true);
            vector<double> interp_error;

            for (int window = 0; window < numWindows; ++window)
            {
                admd = dynamic_cast<CAROM::AdaptiveDMD*> (dmd[window]);
                CAROM_VERIFY(admd);
                double t_init = dmd[window]->getTimeOffset();

                dtc = admd->getTrueDt();
                const CAROM::Matrix* f_snapshots = admd->getInterpolatedSnapshots();
                if (myid == 0)
                {
                    cout << "Verifying Adaptive DMD model #" << window <<
                         " against interpolated snapshots." << endl;
                }
                for (int k = 0; k < f_snapshots->numColumns(); ++k)
                {
                    f_snapshots->getColumn(k, interp_snap);
                    CAROM::Vector* result = admd->predict(t_init+k*dtc);

                    Vector dmd_solution(result->getData(), dim);
                    Vector true_solution(interp_snap.getData(), dim);
                    Vector diff(true_solution.Size());
                    subtract(dmd_solution, true_solution, diff);

                    double tot_diff_norm = sqrt(InnerProduct(MPI_COMM_WORLD, diff, diff));
                    double tot_true_solution_norm = sqrt(InnerProduct(MPI_COMM_WORLD, true_solution,
                                                         true_solution));
                    double rel_error = tot_diff_norm / tot_true_solution_norm;
                    interp_error.push_back(rel_error);

                    if (myid == 0)
                    {
                        cout << "Norm of interpolated snapshot #" << k << " is " <<
                             tot_true_solution_norm << endl;
                        cout << "Absolute error of DMD prediction for interpolated snapshot #" << k <<
                             " is " << tot_diff_norm << endl;
                        cout << "Relative error of DMD prediction for interpolated snapshot #" << k <<
                             " is " << rel_error << endl;
                    }
                    delete result;
                }
                if (myid == 0)
                {
                    csv_db.putDoubleVector(outputPath + "/window" + to_string(
                                               window) + "_interp_error.csv",
                                           interp_error, f_snapshots->numColumns());
                }
                interp_error.clear();
            }
        }
    }

    vector<string> testing_par_list;
    csv_db.getStringVector(string(list_dir) + "/" + test_list + ".csv",
                           testing_par_list, false);
    npar = testing_par_list.size();
    CAROM_VERIFY(npar > 0);

    int num_tests = 0;
    vector<double> prediction_time, prediction_error;

    for (int idx_dataset = 0; idx_dataset < npar; ++idx_dataset)
    {
        stringstream par_ss(testing_par_list[idx_dataset]); // testing DATASET
        string par_dir;
        getline(par_ss, par_dir, ',');
        if (myid == 0)
        {
            cout << "Predicting solution for " << par_dir << " using DMD." << endl;
        }
        vector<string> snap_list;
        csv_db.getStringVector(string(list_dir) + "/" + par_dir + ".csv", snap_list,
                               false);
        CAROM_VERIFY(snap_list.size() > 0);

        vector<double> tvec;
        csv_db.getDoubleVector(string(data_dir) + "/" + par_dir + "/tval.csv", tvec,
                               false);
        CAROM_VERIFY(tvec.size() == snap_list.size());

        vector<int> snap_bound;
        csv_db.getIntegerVector(string(data_dir) + "/" + par_dir + "/" +
                                temporal_idx_list + ".csv", snap_bound, false);
        if (snap_bound.size() > 0)
        {
            snap_bound[0] -= 1;
            snap_bound[1] -= 1;
            CAROM_VERIFY(snap_bound.size() == 2);
            if (myid == 0)
            {
                cout << "Restricting on snapshot #" << snap_bound[0] << " to " << snap_bound[1]
                     << "." << endl;
            }
        }
        else
        {
            snap_bound.push_back(0);
            snap_bound.push_back(snap_list.size()-1);
        }

        int min_idx_snap = -1;
        int max_idx_snap = snap_bound[1];
        int curr_window = 0;
        for (int idx_snap = snap_bound[0]; idx_snap <= snap_bound[1]; ++idx_snap)
        {
            string snap = snap_list[idx_snap]; // STATE
            double tval = tvec[idx_snap];
            string data_filename = string(data_dir) + "/" + par_dir + "/" + snap + "/" +
                                   variable + ".csv"; // path to VAR_NAME.csv
            csv_db.getDoubleArray(data_filename, sample, nelements, idx_state);
            if (myid == 0)
            {
                cout << "State #" << idx_snap << " - " << data_filename
                     << " is read." << endl;
            }

            CAROM::Vector* init_cond = nullptr;
            if (min_idx_snap == -1)
            {
                double curr_indicator_val = (indicator_idx.size() == 0) ? tval :
                                            sample[indicator_idx[0]];
                if (curr_indicator_val >= indicator_val[0])
                {
                    min_idx_snap = idx_snap;
                    if (myid == 0)
                    {
                        cout << "State #" << idx_snap << " - " << data_filename
                             << " is the beginning of window 0." << endl;
                    }

                    dmd_preprocess_timer.Start();
                    if (myid == 0)
                    {
                        cout << "Projecting initial condition at t = " << tval <<
                             " for DMD model #0." << endl;
                    }
                    init_cond = new CAROM::Vector(dim, true);
                    for (int i = 0; i < dim; ++i)
                    {
                        init_cond->item(i) = sample[i];
                    }
                    dmd[curr_window]->projectInitialCondition(init_cond, tval);
                    delete init_cond;
                    dmd_preprocess_timer.Stop();
                }
                else
                {
                    if (myid == 0)
                    {
                        cout << "State #" << idx_snap << " - " << data_filename
                             << " is omitted." << endl;
                    }
                    continue;
                }
            }

            if (myid == 0)
            {
                cout << "Predicting DMD solution #" << idx_snap << " at t = " << tval <<
                     " using DMD model #" << curr_window << endl;
            }
            dmd_prediction_timer.Start();
            CAROM::Vector* result = dmd[curr_window]->predict(tval);
            dmd_prediction_timer.Stop();

            double curr_indicator_val = (indicator_idx.size() == 0) ? tval : result->item(
                                            indicator_idx[curr_window+1]);

            while (curr_window+1 < numWindows
                    && curr_indicator_val >= indicator_val[curr_window+1])
            {
                double t_offset;
                if (indicator_idx.size() == 0)
                {
                    t_offset = indicator_val[curr_window+1];
                }
                else
                {
                    if (myid == 0)
                    {
                        cout << "Indicator state index: " << indicator_idx[curr_window+1] << endl;
                    }
                    if (string(window_endpoint_option) == "bisection")
                    {
                        double t_left = tvec[idx_snap-1];
                        double t_right = tval;

                        for (int k = 0; k < 5; ++k)
                        {
                            t_offset = (t_left + t_right) / 2.0;
                            init_cond = dmd[curr_window]->predict(t_offset);
                            curr_indicator_val = init_cond->item(indicator_idx[curr_window+1]);
                            cout << "t_offset: " << t_offset << endl;
                            cout << "Indicator endpoint: " << indicator_val[curr_window+1] << endl;
                            cout << "Current indicator value: " << curr_indicator_val << endl;
                            if (curr_indicator_val >= indicator_val[curr_window+1])
                            {
                                t_right = t_offset;
                            }
                            else
                            {
                                t_left = t_offset;
                            }
                            delete init_cond;
                        }
                        t_offset = (t_left + t_right) / 2.0;
                    }
                    else if (string(window_endpoint_option) == "left")
                    {
                        t_offset = tvec[idx_snap-1];
                    }
                    else if (string(window_endpoint_option) == "right")
                    {
                        t_offset = tvec[idx_snap];
                    }
                    else
                    {
                        cout << "Invalid window endpoint option." << endl;
                        cout << "Using current time as window endpoint." << endl;
                        t_offset = tvec[idx_snap];
                    }
                }

                if (myid == 0)
                {
                    cout << "Indicator endpoint: " << indicator_val[curr_window+1] << endl;
                    cout << "Projecting initial condition at t = " << t_offset
                         << " for DMD model #" << curr_window+1 << endl;
                    cout << "State #" << idx_snap << " - " << data_filename
                         << " is the beginning of window " << curr_window+1 << "." << endl;
                }

                init_cond = dmd[curr_window]->predict(t_offset);
                dmd[curr_window+1]->projectInitialCondition(init_cond, t_offset);
                delete init_cond;

                delete result;
                curr_window += 1;
                result = dmd[curr_window]->predict(tval);
                curr_indicator_val = (indicator_idx.size() == 0) ? tval : result->item(
                                         indicator_idx[curr_window+1]);
            }

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
                cout << "Norm of true solution at t = " << tval << " is " <<
                     tot_true_solution_norm << endl;
                cout << "Absolute error of DMD solution at t = " << tval << " is " <<
                     tot_diff_norm << endl;
                cout << "Relative error of DMD solution at t = " << tval << " is " << rel_error
                     << endl;
                if (save_csv)
                {
                    csv_db.putDoubleArray(outputPath + "/" + par_dir + "_" + snap +
                                          "_prediction.csv", result->getData(), dim);
                    if (dim < nelements)
                    {
                        csv_db.putDoubleArray(outputPath + "/" + par_dir + "_" + snap + "_state.csv",
                                              sample, dim);
                    }
                }
            }
            delete result;

            if (curr_window == numWindows-1
                    && curr_indicator_val >= indicator_val[numWindows])
            {
                max_idx_snap = idx_snap;
                if (myid == 0)
                {
                    cout << "State #" << idx_snap << " - " << data_filename
                         << " is the end of window " << numWindows-1 << "." << endl;
                }
                break;
            }
        }

        int num_snap = max_idx_snap - min_idx_snap + 1;
        if (myid == 0 && t_final <= 0.0)
        {
            csv_db.putDoubleVector(outputPath + "/" + par_dir + "_prediction_time.csv",
                                   prediction_time, num_snap);
            csv_db.putDoubleVector(outputPath + "/" + par_dir + "_prediction_error.csv",
                                   prediction_error, num_snap);
        }
        prediction_time.clear();
        prediction_error.clear();
        num_tests += (t_final > 0.0) ? 1 : num_snap;
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

    delete[] sample;
    for (int window = 0; window < numWindows; ++window)
    {
        delete dmd[window];
    }

    return 0;
}
