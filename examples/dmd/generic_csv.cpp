/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: DMD on general CSV datasets.
//
// User specify file locations and names by -list LIST_DIR -data DATA_DIR -var VAR_NAME
// 
// File structure:
// 1. LIST_DIR/training_par.csv           -- each row specifies one training DATASET 
// 2. LIST_DIR/testing_par.csv            -- each row specifies one testing DATASET
// 3. LIST_DIR/DATASET.csv                -- each row specifies one STATE in DATASET
// 4. DATA_DIR/DATASET/STATE/VAR_NAME.csv -- each row specifies one value of VAR_NAME of STATE
// 5. DATA_DIR/DATASET/STATE/tval.csv     -- specifies the time instance of STATE
// 6. DATA_DIR/dim.csv                    -- specifies the dimension of VAR_NAME
// 7. DATA_DIR/index.csv                  -- (optional) each row specifies one DOF of VAR_NAME

#include "mfem.hpp"
#include "algo/DMD.h"
#include "algo/AdaptiveDMD.h"
#include "linalg/Vector.h"
#include "linalg/Matrix.h"
#include "utils/HDFDatabase.h"
#include "utils/CSVDatabase.h"
#include <cmath>
#include <iostream>

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
    // 1. Initialize MPI.
    MPI_Session mpi;
    int num_procs = mpi.WorldSize();
    int myid = mpi.WorldRank();

    double t_final = -1.0;
    double dtc = 0.0;
    double ddt = 0.0;
    double dmd_epsilon = -1.0;
    double ef = 1.0;
    int rdim = -1;
    const char *list_dir = "";
    const char *data_dir = "";
    const char *var_name = "var";
    const char *basename = "";
    int precision = 16;
    cout.precision(precision);

    OptionsParser args(argc, argv);
    args.AddOption(&t_final, "-tf", "--t-final",
                   "Final time.");
    args.AddOption(&dtc, "-dtc", "--dtc", 
                   "Fixed (constant) dt.");
    args.AddOption(&ddt, "-ddt", "--dtime-step",
                   "Desired Time step.");
    args.AddOption(&dmd_epsilon, "-dmde", "--dmde",
                   "Parameter epsilon for controlling radial basis functions.");
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
    args.AddOption(&basename, "-o", "--outputfile-name",
                   "Name of the sub-folder to dump files within the run directory.");
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

    std::string outputPath = "run";
    if (std::string(basename) != "") {
        outputPath += "/" + std::string(basename);
    }

    if (mpi.Root()) {
        const char path_delim = '/';
        std::string::size_type pos = 0;
        do {
            pos = outputPath.find(path_delim, pos+1);
            std::string subdir = outputPath.substr(0, pos);
            mkdir(subdir.c_str(), 0777);
        }
        while (pos != std::string::npos);
    }

    CAROM::CSVDatabase* csv_db(new CAROM::CSVDatabase);

    std::string variable = std::string(var_name);
    int nelements = -1;
    csv_db->getIntegerArray(std::string(data_dir) + "/dim.csv", &nelements, 1);
    if (myid == 0)
    {
        cout << "Variable " << var_name << " has dimension " << nelements << "." << endl;
    }

    int dim = nelements;
    std::vector<int> idx_state;
    csv_db->getIntegerVector(std::string(data_dir) + "/index.csv", idx_state, false);
    if (idx_state.size() > 0)
    {
        dim = idx_state.size();
        cout << "Restricting on " << dim << " entries out of " << nelements << "." << endl;
    }

    CAROM::DMD* dmd = nullptr;
    CAROM::AdaptiveDMD* admd = nullptr;
    if (dtc > 0.0){
        dmd = new CAROM::DMD(dim, dtc);
    }
    else
    {
        dmd = new CAROM::AdaptiveDMD(dim, ddt, "G", "LS", dmd_epsilon);
        admd = dynamic_cast<CAROM::AdaptiveDMD*> (dmd);
    }

    std::vector<std::string> training_par_list;
    csv_db->getStringVector(std::string(list_dir) + "/training_par.csv", training_par_list, false);
    int npar = training_par_list.size();

    double* sample = new double[dim];
    StopWatch dmd_training_timer, dmd_prediction_timer;
    dmd_training_timer.Start();

    for (int idx_test = 0; idx_test < npar; ++idx_test) 
    {
        std::string par_dir = training_par_list[idx_test]; // training DATASET
        if (myid == 0)
        {
            cout << "Loading samples for " << par_dir << " to train DMD." << endl;
        }
        std::vector<std::string> snap_list;
        csv_db->getStringVector(std::string(list_dir) + "/" + par_dir + ".csv", snap_list, false);
        int num_snap = snap_list.size();
        double tval = 0.0;
        for (int idx_snap = 0; idx_snap < num_snap; ++idx_snap)
        {
            std::string snap = snap_list[idx_snap]; // STATE
            csv_db->getDoubleArray(std::string(data_dir) + "/" + par_dir + "/" + snap + "/tval.csv", &tval, 1);
            std::string data_filename = std::string(data_dir) + "/" + par_dir + "/" + snap + "/" + variable + ".csv"; // path to VAR_NAME.csv
            csv_db->getDoubleArray(data_filename, sample, nelements, idx_state);
            dmd->takeSample(sample, tval);
        }
        if (myid == 0)
        {
            cout << "Loaded " << num_snap << " samples for " << par_dir << "." << endl;
        }
    }

    if (rdim != -1)
    {
        if (myid == 0)
        {
            cout << "Creating DMD with rdim: " << rdim << endl;
        }
        dmd->train(rdim);
    }
    else if (ef != -1)
    {
        if (myid == 0)
        {
            cout << "Creating DMD with energy fraction: " << ef << endl;
        }
        dmd->train(ef);
    }

    dmd_training_timer.Stop();
    dmd->summary(outputPath);

    CAROM::Vector* result = new CAROM::Vector(dim, true);
    if (admd)
    {
        CAROM_VERIFY(npar == 1); // TODO
        double t_init = dmd->getTimeOffset();
        dtc = admd->getTrueDt();
        const CAROM::Matrix* f_snapshots = admd->getInterpolatedSnapshots();
        CAROM::Vector* isnap = new CAROM::Vector(dim, true);
        if (myid == 0)
        {
            cout << "Verifying Adaptive DMD solution." << endl;
        }
        for (int k = 0; k < f_snapshots->numColumns(); ++k)
        {
            f_snapshots->getColumn(k, isnap);
            result = admd->predict(t_init+k*dtc);

            Vector dmd_solution(result->getData(), dim);
            Vector true_solution(isnap->getData(), dim);
            Vector diff(true_solution.Size());
            subtract(dmd_solution, true_solution, diff);

            double tot_diff_norm = sqrt(InnerProduct(MPI_COMM_WORLD, diff, diff));
            double tot_true_solution_norm = sqrt(InnerProduct(MPI_COMM_WORLD, true_solution, true_solution));

            if (myid == 0)
            {
                cout << "Norm of " << k << "-th interpolated snapshot is " << tot_true_solution_norm << endl;
                cout << "Absolute error of DMD prediction for " << k << "-th interpolated snapshot is " << tot_diff_norm << endl;
                cout << "Relative error of DMD prediction for " << k << "-th interpolated snapshot is " << tot_diff_norm / tot_true_solution_norm << endl;
            }
        }
        delete isnap;
    }

    std::vector<std::string> testing_par_list;
    csv_db->getStringVector(std::string(list_dir) + "/testing_par.csv", testing_par_list, false);
    npar = testing_par_list.size();

    int num_tests = 0;
    CAROM::Vector* init_cond = new CAROM::Vector(dim, true);

    for (int idx_test = 0; idx_test < npar; ++idx_test) 
    {
        std::string par_dir = testing_par_list[idx_test]; // testing DATASET
        if (myid == 0)
        {
            cout << "Predicting solution for " << par_dir << " using DMD." << endl;
        }
        std::vector<std::string> snap_list;
        csv_db->getStringVector(std::string(list_dir) + "/" + par_dir + ".csv", snap_list, false);
        int num_snap = snap_list.size();
        double tval = 0.0;
        for (int idx_snap = 0; idx_snap < num_snap; ++idx_snap)
        {
            std::string snap = snap_list[idx_snap]; // STATE
            csv_db->getDoubleArray(std::string(data_dir) + "/" + par_dir + "/" + snap + "/tval.csv", &tval, 1);
            std::string data_filename = std::string(data_dir) + "/" + par_dir + "/" + snap + "/" + variable + ".csv"; // path to VAR_NAME.csv
            csv_db->getDoubleArray(data_filename, sample, nelements, idx_state);
            if (myid == 0)
            {
                cout << "State " << data_filename << " read." << endl;
            }

            if (idx_snap == 0)
            {
                for (int i = 0; i < dim; ++i)
                {
                    init_cond->item(i) = sample[i]; 
                }
                dmd->projectInitialCondition(init_cond);
                if (t_final > 0.0) // Actual prediction
                {
                    num_tests = 1;
                    if (myid == 0)
                    {
                        cout << "Predicting DMD solution at t = " << t_final << "." << endl;
                    }
                    dmd_prediction_timer.Start();
                    result = dmd->predict(t_final);
                    dmd_prediction_timer.Stop();
                    // TODO: store result
                    return 0;
                }
            }
            else // Verify DMD prediction results against dataset
            {
                if (myid == 0)
                {
                    cout << "Predicting DMD solution #" << idx_snap << " at t = " << tval << "." << endl;
                }
                dmd_prediction_timer.Start();
                result = dmd->predict(tval);
                dmd_prediction_timer.Stop();

                // Calculate the relative error between the DMD final solution and the true solution.
                Vector dmd_solution(result->getData(), result->dim());
                Vector true_solution(sample, dim);
                Vector diff(true_solution.Size());
                subtract(dmd_solution, true_solution, diff);

                double tot_diff_norm = sqrt(InnerProduct(MPI_COMM_WORLD, diff, diff));
                double tot_true_solution_norm = sqrt(InnerProduct(MPI_COMM_WORLD, true_solution, true_solution));

                if (myid == 0)
                {
                    cout << "Norm of true solution at t = " << tval << " is " << tot_true_solution_norm << endl;
                    cout << "Absolute error of DMD solution at t = " << tval << " is " << tot_diff_norm << endl;
                    cout << "Relative error of DMD solution at t = " << tval << " is " << tot_diff_norm / tot_true_solution_norm << endl;
                }

            }
        }
        num_tests += num_snap;
    }

    CAROM_VERIFY(num_tests > 0);

    if (myid == 0)
    {
        printf("Elapsed time for training DMD: %e second\n", dmd_training_timer.RealTime());
        printf("Total elapsed time for predicting DMD: %e second\n", dmd_prediction_timer.RealTime());
        printf("Average elapsed time for predicting DMD: %e second\n", dmd_prediction_timer.RealTime() / num_tests);
    }

    delete[] sample;
    delete result;
    delete init_cond;
    delete dmd;
    return 0;
}
