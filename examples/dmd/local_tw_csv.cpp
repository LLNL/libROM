/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Compile with: make local_tw_csv
//
// Generate CSV databse on heat conduction with: heat_conduction_csv.sh
//
// =================================================================================
//
// Local serial DMD command:
//   mpirun -np 8 local_tw_csv -o local_csv_serial -rdim 16 -dtc 0.01
//
// Final-time prediction error (Last line in run/local_csv_serial/hc_par5_prediction_error.csv):
//   0.000199922
//
// Local time windowing DMD command:
//   mpirun -np 8 local_tw_csv -o local_csv_tw -rdim 16 -nwinsamp 25 -dtc 0.01
//
// Final-time prediction error (Last line in run/local_csv_tw/hc_par5_prediction_error.csv):
//   0.000334656
//
// =================================================================================
//
// Description: Local time windowing DMD on general CSV datasets.
//
// User specify file locations and names by -list LIST_DIR -train-set TRAIN_LIST -test-set TEST_LIST -data DATA_DIR -var VAR_NAME -o OUT_DIR
//
// File structure:
// 1. LIST_DIR/TRAIN_LIST.csv             -- each row specifies one training DATASET
// 2. LIST_DIR/TEST_LIST.csv              -- each row specifies one testing DATASET
// 3. LIST_DIR/DATASET.csv                -- each row specifies one STATE in DATASET
// 4. DATA_DIR/dim.csv                    -- specifies the dimension of VAR_NAME
// 5. DATA_DIR/DATASET/tval.csv           -- specifies the time instances
// 6. DATA_DIR/DATASET/STATE/VAR_NAME.csv -- each row specifies one value of VAR_NAME of STATE
// 7. DATA_DIR/DATASET/TEMPORAL_IDX.csv   -- (optional) specifies the first and last temporal index in DATASET
// 8. DATA_DIR/SPATIAL_IDX.csv            -- (optional) each row specifies one spatial index of VAR_NAME
// 9. run/OUT_DIR/indicator_val.csv       -- (optional) each row specifies one indicator endpoint value

#include "mfem.hpp"
#include "algo/DMD.h"
#include "algo/AdaptiveDMD.h"
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
    int numWindows = 0;
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
    args.AddOption(&temporal_idx_list, "-t-idx", "--temporal0index",
                   "Name of the file indicating bound of temporal indices.");
    args.AddOption(&spatial_idx_list, "-x-idx", "--spatial-index",
                   "Name of the file indicating spatial indices.");
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

    CAROM_VERIFY((dtc > 0.0 || ddt > 0.0) && !(dtc > 0.0 && ddt > 0.0));

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
        cout << "Variable " << var_name << " has dimension " << nelements << "." << endl;
    }

    int dim = nelements;
    vector<int> idx_state;
    csv_db.getIntegerVector(string(data_dir) + "/" + string(spatial_idx_list) + ".csv", idx_state, false);
    if (idx_state.size() > 0)
    {
        dim = idx_state.size();
        if (myid == 0)
        {
            cout << "Restricting on " << dim << " entries out of " << nelements << "." << endl;
        }
    }

    vector<double> indicator_val;
    if (!train || numWindows > 0)
    {
        csv_db.getDoubleVector(string(outputPath) + "/indicator_val.csv", indicator_val, false);
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
                cout << "Read indicator range partition with " << numWindows << " windows." << endl;
            }
        }
    }

    int npar = 0;
    int num_train_snap;
    vector<string> training_par_list;
    if (train)
    {
        npar = csv_db.getLineCount(string(list_dir) + "/" + train_list + ".csv");
        CAROM_VERIFY(npar == 1);

        csv_db.getStringVector(string(list_dir) + "/" + train_list + ".csv", training_par_list, false);
        stringstream par_ss(training_par_list[0]); // training DATASET
        string par_dir;
        getline(par_ss, par_dir, ',');

        vector<int> snap_bound;
        csv_db.getIntegerVector(string(data_dir) + "/" + par_dir + "/" + temporal_idx_list + ".csv", snap_bound, false);
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

        CAROM_VERIFY(windowOverlapSamples < windowNumSamples);
        numWindows = (windowNumSamples < infty) ? ceil(num_train_snap / windowNumSamples) : 1;
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

    vector<CAROM::DMD*> dmd;
    dmd.assign(numWindows, nullptr);
    for (int window = 0; window < numWindows; ++window)
    {
        if (train)
        {
            if (ddt > 0.0)
            {
                dmd[window] = new CAROM::AdaptiveDMD(dim, ddt, string(rbf), string(interp_method), admd_closest_rbf_val);
            }
            else
            {
                dmd[window] = new CAROM::DMD(dim, dtc);
            }
        }
        else
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
        csv_db.getStringVector(string(list_dir) + "/" + par_dir + ".csv", snap_list, false);

        if (myid == 0)
        {
            cout << "Loading samples for " << par_dir << " to train DMD." << endl;
        }

        vector<double> tvec;
        csv_db.getDoubleVector(string(data_dir) + "/" + par_dir + "/tval.csv", tvec, false);
        CAROM_VERIFY(tvec.size() == snap_list.size());

        vector<int> snap_bound;
        csv_db.getIntegerVector(string(data_dir) + "/" + par_dir + "/" + temporal_idx_list + ".csv", snap_bound, false);
        if (snap_bound.size() > 0)
        {
            snap_bound[0] -= 1;
            snap_bound[1] -= 1;
            CAROM_VERIFY(snap_bound.size() == 2);
            if (myid == 0)
            {
                cout << "Restricting on snapshot #" << snap_bound[0] << " to " << snap_bound[1] << "." << endl;
            }
        }
        else
        {
            snap_bound.push_back(0);
            snap_bound.push_back(snap_list.size()-1);
        }

        int curr_window = 0;
        int overlap_count = 0;
        for (int idx_snap = snap_bound[0]; idx_snap <= snap_bound[1]; ++idx_snap)
        {
            string snap = snap_list[idx_snap]; // STATE
            double tval = tvec[idx_snap];

            if (idx_snap == snap_bound[0])
            {
                indicator_val.push_back(tval);
            }

            string data_filename = string(data_dir) + "/" + par_dir + "/" + snap + "/" + variable + ".csv"; // path to VAR_NAME.csv
            csv_db.getDoubleArray(data_filename, sample, nelements, idx_state);
            dmd[curr_window]->takeSample(sample, tval);
            if (overlap_count > 0)
            {
                dmd[curr_window-1]->takeSample(sample, tval);
                overlap_count -= 1;
            }
            if (curr_window+1 < numWindows && idx_snap+1 <= snap_bound[1])
            {
                bool new_window = false;
                if (windowNumSamples < infty)
                {
                    new_window = (idx_snap >= snap_bound[0] + (curr_window+1)*windowNumSamples);
                }
                else
                {
                    new_window = (tval >= indicator_val[curr_window+1]);
                }
                if (new_window)
                {
                    overlap_count = windowOverlapSamples;
                    curr_window += 1;
                    if (windowNumSamples < infty)
                    {
                        indicator_val.push_back(tval);
                    }
                    dmd[curr_window]->takeSample(sample, tval);
                }
            }
        }

        if (myid == 0)
        {
            cout << "Loaded " << num_train_snap << " samples for " << par_dir << "." << endl;
            if (windowNumSamples < infty)
            {
                cout << "Created new indicator range partition." << endl;
                csv_db.putDoubleVector(string(outputPath) + "/indicator_val.csv", indicator_val, numWindows);
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
                    cout << "Creating DMD model #" << window << " with energy fraction: " << ef << endl;
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
                double t_init = dmd[window]->getTimeOffset();

                dtc = admd->getTrueDt();
                const CAROM::Matrix* f_snapshots = admd->getInterpolatedSnapshots();
                if (myid == 0)
                {
                    cout << "Verifying Adaptive DMD model #" << window << " against interpolated snapshots." << endl;
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
                    double tot_true_solution_norm = sqrt(InnerProduct(MPI_COMM_WORLD, true_solution, true_solution));
                    double rel_error = tot_diff_norm / tot_true_solution_norm;
                    interp_error.push_back(rel_error);

                    if (myid == 0)
                    {
                        cout << "Norm of interpolated snapshot #" << k << " is " << tot_true_solution_norm << endl;
                        cout << "Absolute error of DMD prediction for interpolated snapshot #" << k << " is " << tot_diff_norm << endl;
                        cout << "Relative error of DMD prediction for interpolated snapshot #" << k << " is " << rel_error << endl;
                    }
                    delete result;
                }
                if (myid == 0)
                {
                    csv_db.putDoubleVector(outputPath + "/window" + to_string(window) + "_interp_error.csv",
                                           interp_error, f_snapshots->numColumns());
                }
                interp_error.clear();
            }
        }
    }

    vector<string> testing_par_list;
    npar = csv_db.getLineCount(string(list_dir) + "/" + test_list + ".csv");
    csv_db.getStringVector(string(list_dir) + "/" + test_list + ".csv", testing_par_list, false);

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
        csv_db.getStringVector(string(list_dir) + "/" + par_dir + ".csv", snap_list, false);

        vector<double> tvec;
        csv_db.getDoubleVector(string(data_dir) + "/" + par_dir + "/tval.csv", tvec, false);
        CAROM_VERIFY(tvec.size() == snap_list.size());

        vector<int> snap_bound;
        csv_db.getIntegerVector(string(data_dir) + "/" + par_dir + "/" + temporal_idx_list + ".csv", snap_bound, false);
        if (snap_bound.size() > 0)
        {
            snap_bound[0] -= 1;
            snap_bound[1] -= 1;
            CAROM_VERIFY(snap_bound.size() == 2);
            if (myid == 0)
            {
                cout << "Restricting on snapshot #" << snap_bound[0] << " to " << snap_bound[1] << "." << endl;
            }
        }
        else
        {
            snap_bound.push_back(0);
            snap_bound.push_back(snap_list.size()-1);
        }

        int num_snap = snap_bound[1] - snap_bound[0] + 1;
        int curr_window = 0;
        for (int idx_snap = snap_bound[0]; idx_snap <= snap_bound[1]; ++idx_snap)
        {
            string snap = snap_list[idx_snap]; // STATE
            double tval = tvec[idx_snap];
            string data_filename = string(data_dir) + "/" + par_dir + "/" + snap + "/" + variable + ".csv"; // path to VAR_NAME.csv
            csv_db.getDoubleArray(data_filename, sample, nelements, idx_state);
            if (myid == 0)
            {
                cout << "State " << data_filename << " read." << endl;
            }

            if (idx_snap == 0)
            {
                dmd_preprocess_timer.Start();
                CAROM::Vector* init_cond = nullptr;
                for (int window = 0; window < numWindows; ++window)
                {
                    if (myid == 0)
                    {
                        cout << "Projecting initial condition at t = " << indicator_val[window] << " for DMD model #" << window << endl;
                    }
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
                        init_cond = dmd[window-1]->predict(indicator_val[window]);
                    }
                    dmd[window]->projectInitialCondition(init_cond);
                    delete init_cond;
                }
                dmd_preprocess_timer.Stop();
            }

            if (t_final > 0.0) // Actual prediction without true solution for comparison
            {
                num_tests += 1;
                while (curr_window+1 < numWindows && t_final > indicator_val[curr_window+1])
                {
                    curr_window += 1;
                }
                if (myid == 0)
                {
                    cout << "Predicting DMD solution at t = " << t_final << " using DMD model #" << curr_window << endl;
                }
                dmd_prediction_timer.Start();
                CAROM::Vector* result = dmd[curr_window]->predict(t_final);
                dmd_prediction_timer.Stop();
                if (myid == 0)
                {
                    csv_db.putDoubleArray(outputPath + "/" + par_dir + "_final_time_prediction.csv", result->getData(), dim);
                }
                idx_snap = snap_bound[1]+1; // escape for-loop over idx_snap
                delete result;
            }
            else // Verify DMD prediction results against dataset
            {
                while (curr_window+1 < numWindows && tval > indicator_val[curr_window+1])
                {
                    curr_window += 1;
                }
                if (myid == 0)
                {
                    cout << "Predicting DMD solution #" << idx_snap << " at t = " << tval << " using DMD model #" << curr_window << endl;
                }
                dmd_prediction_timer.Start();
                CAROM::Vector* result = dmd[curr_window]->predict(tval);
                dmd_prediction_timer.Stop();

                // Calculate the relative error between the DMD final solution and the true solution.
                Vector dmd_solution(result->getData(), result->dim());
                Vector true_solution(sample, dim);
                Vector diff(true_solution.Size());
                subtract(dmd_solution, true_solution, diff);

                double tot_diff_norm = sqrt(InnerProduct(MPI_COMM_WORLD, diff, diff));
                double tot_true_solution_norm = sqrt(InnerProduct(MPI_COMM_WORLD, true_solution, true_solution));
                double rel_error = tot_diff_norm / tot_true_solution_norm;

                prediction_time.push_back(tval);
                prediction_error.push_back(rel_error);

                if (myid == 0)
                {
                    cout << "Norm of true solution at t = " << tval << " is " << tot_true_solution_norm << endl;
                    cout << "Absolute error of DMD solution at t = " << tval << " is " << tot_diff_norm << endl;
                    cout << "Relative error of DMD solution at t = " << tval << " is " << rel_error << endl;
                    if (save_csv)
                    {
                        csv_db.putDoubleArray(outputPath + "/" + par_dir + "_" + snap + "_prediction.csv", result->getData(), dim);
                        if (dim < nelements)
                        {
                            csv_db.putDoubleArray(outputPath + "/" + par_dir + "_" + snap + "_state.csv", sample, dim);
                        }
                    }
                }
                delete result;
            }
        }
        if (myid == 0 && t_final <= 0.0)
        {
            csv_db.putDoubleVector(outputPath + "/" + par_dir + "_prediction_time.csv", prediction_time, num_snap);
            csv_db.putDoubleVector(outputPath + "/" + par_dir + "_prediction_error.csv", prediction_error, num_snap);
        }
        prediction_time.clear();
        prediction_error.clear();
        num_tests = (t_final > 0.0) ? num_tests + 1 : num_tests + num_snap;
    }

    CAROM_VERIFY(num_tests > 0);

    if (myid == 0)
    {
        printf("Elapsed time for training DMD: %e second\n", dmd_training_timer.RealTime());
        printf("Elapsed time for preprocessing DMD: %e second\n", dmd_preprocess_timer.RealTime());
        printf("Total elapsed time for predicting DMD: %e second\n", dmd_prediction_timer.RealTime());
        printf("Average elapsed time for predicting DMD: %e second\n", dmd_prediction_timer.RealTime() / num_tests);
    }

    delete[] sample;
    for (int window = 0; window < numWindows; ++window)
    {
        delete dmd[window];
    }

    return 0;
}
