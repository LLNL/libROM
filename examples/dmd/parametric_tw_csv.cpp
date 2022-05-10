/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Parametric time windowing DMD on general CSV datasets.
//
// User specify file locations and names by -list LIST_DIR -data DATA_DIR -var VAR_NAME -o OUT_DIR
//
// File structure:
// 1. LIST_DIR/training_par.csv           -- each row specifies one training DATASET
// 2. LIST_DIR/testing_par.csv            -- each row specifies one testing DATASET
// 3. LIST_DIR/DATASET.csv                -- each row specifies one STATE in DATASET
// 4. DATA_DIR/DATASET/STATE/VAR_NAME.csv -- each row specifies one value of VAR_NAME of STATE
// 5. DATA_DIR/DATASET/STATE/tval.csv     -- specifies the time instance of STATE
// 6. DATA_DIR/dim.csv                    -- specifies the dimension of VAR_NAME
// 7. DATA_DIR/index.csv                  -- (optional) each row specifies one DOF of VAR_NAME
// 8. run/OUT_DIR/indicator_val.csv       -- (optional) each row specifies one indicator endpoint value

#include "mfem.hpp"
#include "algo/DMD.h"
#include "algo/AdaptiveDMD.h"
#include "algo/ParametricDMD.h"
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
    double t_final = -1.0;
    double dtc = 0.0;
    double ddt = 0.0;
    double admd_closest_rbf_val = 0.9;
    double ef = 0.9999;
    int rdim = -1;
    int windowNumSamples = infty;
    int windowOverlapSamples = 0;
    const char *list_dir = "../data/hc_test2";
    const char *data_dir = "../data/hc_data";
    const char *var_name = "sol";
    bool offline = false;
    bool online = false;
    double pdmd_closest_rbf_val = 0.9;
    bool predict = false;
    const char *basename = "";
    bool save_csv = false;

    OptionsParser args(argc, argv);
    args.AddOption(&t_final, "-tf", "--t-final",
                   "Final time.");
    args.AddOption(&dtc, "-dtc", "--dtc",
                   "Fixed (constant) dt.");
    args.AddOption(&ddt, "-ddt", "--dtime-step",
                   "Desired Time step.");
    args.AddOption(&admd_closest_rbf_val, "-admde", "--admde",
                   "Parameter for controlling radial basis function value for the closest two time steps values in adaptive DMD.");
    args.AddOption(&ef, "-ef", "--energy-fraction",
                   "Energy fraction for DMD.");
    args.AddOption(&rdim, "-rdim", "--rdim",
                   "Reduced dimension for DMD.");
    args.AddOption(&windowNumSamples, "-nwinsamp", "--numwindowsamples",
                   "Number of samples in DMD windows.");
    args.AddOption(&windowOverlapSamples, "-nwinover", "--numwindowoverlap",
                   "Number of samples for DMD window overlap.");
    args.AddOption(&list_dir, "-list", "--list-directory",
                   "Location of training and testing data list.");
    args.AddOption(&data_dir, "-data", "--data-directory",
                   "Location of training and testing data.");
    args.AddOption(&var_name, "-var", "--variable-name",
                   "Name of variable.");
    args.AddOption(&offline, "-offline", "--offline", "-no-offline", "--no-offline",
                   "Enable or disable the offline phase.");
    args.AddOption(&online, "-online", "--online", "-no-online", "--no-online",
                   "Enable or disable the online phase.");
    args.AddOption(&pdmd_closest_rbf_val, "-pdmde", "--pdmde",
                   "Parameter for controlling radial basis function value for the closest two parameter points in interpolating DMD.");
    args.AddOption(&predict, "-predict", "--predict", "-no-predict", "--no-predict",
                   "Enable or disable DMD prediction.");
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

    CAROM_VERIFY(!(offline && online) && (offline || online));
    CAROM_VERIFY((dtc > 0.0 || ddt > 0.0) && !(dtc > 0.0 && ddt > 0.0));
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

    string variable = string(var_name);
    int nelements = -1;
    csv_db.getIntegerArray(string(data_dir) + "/dim.csv", &nelements, 1);
    if (myid == 0)
    {
        cout << "Variable " << var_name << " has dimension " << nelements << "." << endl;
    }

    int dim = nelements;
    vector<int> idx_state;
    csv_db.getIntegerVector(string(data_dir) + "/index.csv", idx_state, false);
    if (idx_state.size() > 0)
    {
        dim = idx_state.size();
        if (myid == 0)
        {
            cout << "Restricting on " << dim << " entries out of " << nelements << "." << endl;
        }
    }

    int numWindows = 1;
    vector<double> indicator_val;
    csv_db.getDoubleVector(string(outputPath) + "/indicator_val.csv", indicator_val, false);
    if (indicator_val.size() > 0)
    {
        numWindows = indicator_val.size();
        windowNumSamples = infty;
        if (myid == 0)
        {
            cout << "Read indicator range partition with " << numWindows << " windows." << endl;
        }
    }

    vector<string> training_par_list, testing_par_list; // DATASET info
    vector<string> par_dir_list; // DATASET name
    vector<CAROM::Vector*> par_vectors; // DATASET param
    vector<int> num_train_snap; // DATASET size
    vector<double> indicator_init, indicator_last; // DATASET indicator range

    int npar = csv_db.getLineCount(string(list_dir) + "/training_par.csv");
    if (myid == 0)
    {
        cout << "Loading " << npar << " training datasets." << endl;
    }

    csv_db.getStringVector(string(list_dir) + "/training_par.csv", training_par_list, false);
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
        vector<string> snap_list;
        csv_db.getStringVector(string(list_dir) + "/" + par_dir + ".csv", snap_list, false);
        num_train_snap.push_back(snap_list.size());

        for (int par_order = 0; par_order < dpar; ++par_order)
        {
            curr_par->item(par_order) = stod(par_info[par_order+1]);
        }
        par_vectors.push_back(curr_par);

        if (offline)
        {
            double tval = 0.0;
            csv_db.getDoubleArray(string(data_dir) + "/" + par_dir + "/" + snap_list[0] + "/tval.csv", &tval, 1);
            indicator_init.push_back(tval);

            csv_db.getDoubleArray(string(data_dir) + "/" + par_dir + "/" + snap_list[snap_list.size()-1] + "/tval.csv", &tval, 1);
            indicator_last.push_back(tval);
        }

        CAROM_VERIFY(windowOverlapSamples < windowNumSamples);
        if (windowNumSamples < infty && indicator_val.size() == 0)
        {
            double indicator_min = *min_element(indicator_init.begin(), indicator_init.end());
            double indicator_max = *max_element(indicator_last.begin(), indicator_last.end());
            numWindows = ceil((indicator_max - indicator_min) / (dt_est * windowNumSamples));
            for (int window = 0; window < numWindows; ++window)
            {
                indicator_val.push_back(indicator_min + dt_est * windowNumSamples * window);
            }
            if (myid == 0)
            {
                cout << "Created new indicator range partition." << endl;
                csv_db.putDoubleVector(string(outputPath) + "/indicator_val.csv", indicator_val, numWindows);
            }
        }
    }

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

    StopWatch dmd_training_timer, dmd_prediction_timer;
    vector<vector<CAROM::DMD*>> dmd;
    vector<CAROM::DMD*> dmd_w;
    double* sample = new double[dim];

    if (offline)
    {
        dmd_training_timer.Start();

        for (int window = 0; window < numWindows; ++window)
        {
            dmd_w.assign(npar, nullptr);
            for (int idx_dataset = 0; idx_dataset < npar; ++idx_dataset)
            {
                if (ddt > 0.0)
                {
                    dmd_w[idx_dataset] = new CAROM::AdaptiveDMD(dim, ddt, "G", "LS", admd_closest_rbf_val);
                }
                else
                {
                    dmd_w[idx_dataset] = new CAROM::DMD(dim, dtc);
                }
            }
            dmd.push_back(dmd_w);
            dmd_w.clear();
        }

        for (int idx_dataset = 0; idx_dataset < npar; ++idx_dataset)
        {
            string par_dir = par_dir_list[idx_dataset];
            if (myid == 0)
            {
                cout << "Loading samples for " << par_dir << " to train DMD." << endl;
            }
            vector<string> snap_list;
            csv_db.getStringVector(string(list_dir) + "/" + par_dir + ".csv", snap_list, false);

            int curr_window = 0;
            int overlap_count = 0;
            for (int idx_snap = 0; idx_snap < num_train_snap[idx_dataset]; ++idx_snap)
            {
                string snap = snap_list[idx_snap]; // STATE
                double tval = 0.0;
                csv_db.getDoubleArray(string(data_dir) + "/" + par_dir + "/" + snap + "/tval.csv", &tval, 1);
                string data_filename = string(data_dir) + "/" + par_dir + "/" + snap + "/" + variable + ".csv"; // path to VAR_NAME.csv
                csv_db.getDoubleArray(data_filename, sample, nelements, idx_state);
                dmd[curr_window][idx_dataset]->takeSample(sample, tval);
                if (overlap_count > 0)
                {
                    dmd[curr_window-1][idx_dataset]->takeSample(sample, tval);
                    overlap_count -= 1;
                }
                if (curr_window+1 < numWindows && idx_snap+1 < num_train_snap[idx_dataset] &&
                    tval > indicator_val[curr_window+1] - dt_est / 100.0) // a rough estimate to correct the precision of the indicator range partition
                {
                    overlap_count = windowOverlapSamples;
                    curr_window += 1;
                    if (windowNumSamples < infty)
                    {
                        indicator_val.push_back(tval);
                    }
                    dmd[curr_window][idx_dataset]->takeSample(sample, tval);
                }
            }

            if (myid == 0)
            {
                cout << "Loaded " << num_train_snap[idx_dataset] << " samples for " << par_dir << "." << endl;
            }

            for (int window = 0; window < numWindows; ++window)
            {
                if (rdim != -1)
                {
                    if (myid == 0)
                    {
                        cout << "Creating DMD model #" << window << " with rdim: " << rdim << endl;
                    }
                    dmd[window][idx_dataset]->train(rdim);
                }
                else if (ef != -1)
                {
                    if (myid == 0)
                    {
                        cout << "Creating DMD model #" << window << " with energy fraction: " << ef << endl;
                    }
                    dmd[window][idx_dataset]->train(ef);
                }
                dmd[window][idx_dataset]->save(outputPath + "/window" + to_string(window) + "_par" + to_string(idx_dataset));
                if (myid == 0)
                {
                    dmd[window][idx_dataset]->summary(outputPath + "/window" + to_string(window) + "_par" + to_string(idx_dataset));
                }
            } // escape for-loop over window
        } // escape for-loop over idx_dataset
        dmd_training_timer.Stop();
    } // escape if-statement of offline

    CAROM::Vector* curr_par = new CAROM::Vector(dpar, false);

    if (online)
    {
        par_dir_list.clear();

        dmd_training_timer.Start();
        npar = csv_db.getLineCount(string(list_dir) + "/testing_par.csv");
        if (myid == 0)
        {
            cout << "Loading " << npar << " testing datasets." << endl;
        }

        csv_db.getStringVector(string(list_dir) + "/testing_par.csv", testing_par_list, false);
        dmd_w.assign(npar, nullptr);
        dmd.assign(numWindows, dmd_w);

        int num_tests = 0;
        CAROM::Vector* init_cond = new CAROM::Vector(dim, true);
        vector<double> prediction_time, prediction_error;

        for (int idx_dataset = 0; idx_dataset < npar; ++idx_dataset)
        {
            stringstream par_ss(testing_par_list[idx_dataset]); // testing DATASET
            vector<string> par_info;
            string par_entry;
            while (getline(par_ss, par_entry, ','))
            {
                par_info.push_back(par_entry);
            }

            CAROM_VERIFY(dpar == par_info.size() - 1);

            string par_dir = par_info[0];
            par_dir_list.push_back(par_dir);
            vector<string> snap_list;
            csv_db.getStringVector(string(list_dir) + "/" + par_dir + ".csv", snap_list, false);
            int num_snap = snap_list.size();

            for (int par_order = 0; par_order < dpar; ++par_order)
            {
                curr_par->item(par_order) = stod(par_info[par_order+1]);
            }

            string snap = snap_list[0]; // STATE
            double tval = 0.0;
            csv_db.getDoubleArray(string(data_dir) + "/" + par_dir + "/" + snap + "/tval.csv", &tval, 1);
            string data_filename = string(data_dir) + "/" + par_dir + "/" + snap + "/" + variable + ".csv"; // path to VAR_NAME.csv
            csv_db.getDoubleArray(data_filename, sample, nelements, idx_state);
            for (int window = 0; window < numWindows; ++window)
            {
                std::vector<std::string> dmd_paths;
                for (int idx_trainset = 0; idx_trainset < par_vectors.size(); ++idx_trainset)
                {
                    dmd_paths.push_back(outputPath + "/window" + to_string(window) + "_par" + to_string(idx_trainset));
                }

                dmd[window][idx_dataset] = getParametricDMD(par_vectors, dmd_paths, curr_par, "G", "LS", pdmd_closest_rbf_val);
                if (window == 0)
                {
                    for (int i = 0; i < dim; ++i)
                    {
                        init_cond->item(i) = sample[i];
                    }
                }
                else
                {
                    init_cond = dmd[window-1][idx_dataset]->predict(indicator_val[window]);
                }
                if (myid == 0)
                {
                    cout << "Projecting initial condition at t = " << indicator_val[window] << " for DMD model #" << window << endl;
                }
                dmd[window][idx_dataset]->projectInitialCondition(init_cond);
            } // escape for-loop over window

        } // escape for-loop over idx_dataset
        delete init_cond;
        dmd_training_timer.Stop();
    } // escape if-statement of online

    if (online || predict)
    {
        int num_tests = 0;
        vector<double> prediction_time, prediction_error;

        dmd_prediction_timer.Start();

        for (int idx_dataset = 0; idx_dataset < npar; ++idx_dataset)
        {
            string par_dir = par_dir_list[idx_dataset];
            if (myid == 0)
            {
                cout << "Predicting solution for " << par_dir << " using DMD." << endl;
            }
            vector<string> snap_list;
            csv_db.getStringVector(string(list_dir) + "/" + par_dir + ".csv", snap_list, false);
            int num_snap = snap_list.size();

            int curr_window = 0;
            double tval = 0.0;
            for (int idx_snap = 0; idx_snap < num_snap; ++idx_snap)
            {
                string snap = snap_list[idx_snap]; // STATE
                csv_db.getDoubleArray(string(data_dir) + "/" + par_dir + "/" + snap + "/tval.csv", &tval, 1);
                string data_filename = string(data_dir) + "/" + par_dir + "/" + snap + "/" + variable + ".csv"; // path to VAR_NAME.csv
                csv_db.getDoubleArray(data_filename, sample, nelements, idx_state);
                if (myid == 0)
                {
                    cout << "State " << data_filename << " read." << endl;
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
                    CAROM::Vector* result = dmd[curr_window][idx_dataset]->predict(t_final);
                    if (myid == 0)
                    {
                        csv_db.putDoubleArray(outputPath + "/" + par_dir + "_final_time_prediction.csv", result->getData(), dim);
                    }
                    idx_snap = num_snap; // escape for-loop over idx_snap
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
                    CAROM::Vector* result = dmd[curr_window][idx_dataset]->predict(tval);

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

        dmd_prediction_timer.Stop();
        CAROM_VERIFY(num_tests > 0);

        if (myid == 0)
        {
            printf("Elapsed time for training DMD: %e second\n", dmd_training_timer.RealTime());
            printf("Total elapsed time for predicting DMD: %e second\n", dmd_prediction_timer.RealTime());
            printf("Average elapsed time for predicting DMD: %e second\n", dmd_prediction_timer.RealTime() / num_tests);
        }
    }

    delete[] sample;
    delete curr_par;
    for (int window = 0; window < numWindows; ++window)
    {
        for (int idx_dataset = 0; idx_dataset < npar; ++idx_dataset)
        {
            delete dmd[window][idx_dataset];
        }
    }

    return 0;
}
