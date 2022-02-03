/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: DMD on poreCollapse.

#include "mfem.hpp"
#include "algo/DMD.h"
#include "algo/AdaptiveDMD.h"
#include "linalg/Vector.h"
#include "linalg/Matrix.h"
#include "utils/HDFDatabase.h"
#include "utils/CSVDatabase.h"
#include <cmath>
#include <iostream>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
    // 1. Initialize MPI.
    MPI_Session mpi;
    int num_procs = mpi.WorldSize();
    int myid = mpi.WorldRank();

    double t_init = 0.0;
    double t_final = -1.0;
    double dtc = 0.0;
    double ddt = 0.00005;
    double dmd_epsilon = -1.0;
    double ef = 1.0;
    int rdim = -1;
    const char *list_dir = "/usr/workspace/nlrom/poreCollapse/libROM_data/pore_collapse_list";
    const char *data_dir = "/usr/workspace/nlrom/poreCollapse/libROM_data/pore_collapse_data";
    const char *var_name = "tkelv";
    int precision = 16;
    cout.precision(precision);

    OptionsParser args(argc, argv);
    args.AddOption(&t_init, "-ti", "--t-init",
                   "Start time.");
    args.AddOption(&t_final, "-tf", "--t-final",
                   "Final time.");
    args.AddOption(&dtc, "-dtc", "--dtc", 
                   "Fixed (constant) dt.");
    args.AddOption(&ddt, "-ddt", "--dtime-step",
                   "Desired Time step.");
    args.AddOption(&dmd_epsilon, "-dmde", "--dmde",
                   "DMD epsilon.");
    args.AddOption(&ef, "-ef", "--energy_fraction",
                   "Energy fraction for DMD.");
    args.AddOption(&rdim, "-rdim", "--rdim",
                   "Reduced dimension for DMD.");
    args.AddOption(&list_dir, "-list", "--list_dir",
                   "Location of training and testing data list.");
    args.AddOption(&data_dir, "-data", "--data_dir",
                   "Location of training and testing data.");
    args.AddOption(&var_name, "-var", "--var_name",
                   "Name of variable.");
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

    CAROM::CSVDatabase* csv_db(new CAROM::CSVDatabase);

    std::string variable = std::string(var_name);
    int nelements = -1;
    csv_db->getIntegerArray(std::string(data_dir) + "/dim.txt", &nelements, 1);
    if (myid == 0)
    {
        cout << "Variable " << var_name << " has dimension " << nelements << "." << endl;
    }

    int dim = nelements;
    std::vector<int> idx_state;
    csv_db->getIntegerArray(std::string(data_dir) + "/index.csv", idx_state, false);
    if (idx_state.size() > 0)
    {
        cout << "Restricting on " << idx_state.size() << " entries out of " << dim << "." << endl;
        dim = idx_state.size();
    }

    CAROM::DMD* dmd = nullptr;
    CAROM::AdaptiveDMD* admd = nullptr;
    if (dtc > 0.0){
        dmd = new CAROM::DMD(dim);
    }
    else
    {
        dmd = new CAROM::AdaptiveDMD(dim, ddt, "LS", "G", dmd_epsilon);
        admd = dynamic_cast<CAROM::AdaptiveDMD*> (dmd);
    }

    std::vector<std::string> training_par_list;
    csv_db->getStringList(std::string(list_dir) + "/training_gpa", training_par_list, false);
    int npar = training_par_list.size();

    double* sample = new double[dim];
    StopWatch dmd_training_timer, dmd_prediction_timer;
    dmd_training_timer.Start();

    for (int idx_test = 0; idx_test < npar; ++idx_test) 
    {
        std::string par_dir = training_par_list[idx_test]; // *gpa
        if (myid == 0)
        {
            cout << "Loading samples for " << par_dir << " to train DMD." << endl;
        }
        std::vector<std::string> snap_list;
        csv_db->getStringList(std::string(list_dir) + "/" + par_dir, snap_list, false);
        int num_snap = snap_list.size();
        double tval = 0.0;
        for (int idx_snap = 0; idx_snap < num_snap; ++idx_snap)
        {
            std::string snap = snap_list[idx_snap]; // run_036.*
            csv_db->getDoubleArray(std::string(data_dir) + "/" + par_dir + "/" + snap + "/tval.txt", &tval, 1);
            if (idx_snap == 0)
            {
                t_init = tval;
            }
            std::string data_filename = std::string(data_dir) + "/" + par_dir + "/" + snap + "/" + variable + ".csv"; // {zone_*,tkelv}.csv
            csv_db->getDoubleArray(data_filename, sample, nelements, idx_state);
            dmd->takeSample(sample, tval - t_init);
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

    CAROM::Vector* result = new CAROM::Vector(dim, true);
    if (admd)
    {
        MFEM_VERIFY(npar == 1, "Adaptive DMD only works with 1 training parameter.");
        dtc = admd->getTruedt();
        const CAROM::Matrix* f_snapshots = admd->getInterpolatedSnapshots();
        CAROM::Vector* isnap = new CAROM::Vector(dim, true);
        if (myid == 0)
        {
            cout << "Verifying Adaptive DMD solution." << endl;
        }
        for (int k = 0; k < f_snapshots->numColumns(); ++k)
        {
            double tval = t_init + k * dtc;
            f_snapshots->getColumn(k, isnap);
            result = admd->predict(k);

            Vector dmd_solution(result->getData(), dim);
            Vector true_solution(isnap->getData(), dim);
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
        delete isnap;
    }

    std::vector<std::string> testing_par_list;
    csv_db->getStringList(std::string(list_dir) + "/testing_gpa", testing_par_list, false);
    npar = testing_par_list.size();

    int num_tests = 0;
    CAROM::Vector* init_cond = new CAROM::Vector(dim, true);

    for (int idx_test = 0; idx_test < npar; ++idx_test) 
    {
        std::string par_dir = testing_par_list[idx_test]; // *gpa
        if (myid == 0)
        {
            cout << "Predicting solution for " << par_dir << " using DMD." << endl;
        }
        std::vector<std::string> snap_list;
        csv_db->getStringList(std::string(list_dir) + "/" + par_dir, snap_list, false);
        int num_snap = snap_list.size();
        double tval = 0.0;
        for (int idx_snap = 0; idx_snap < num_snap; ++idx_snap)
        {
            std::string snap = snap_list[idx_snap]; // run_036.*
            csv_db->getDoubleArray(std::string(data_dir) + "/" + par_dir + "/" + snap + "/tval.txt", &tval, 1);
            std::string data_filename = std::string(data_dir) + "/" + par_dir + "/" + snap + "/" + variable + ".csv"; // {zone_*,tkelv}.csv
            csv_db->getDoubleArray(data_filename, sample, nelements, idx_state);
            if (myid == 0)
            {
                cout << "State " << data_filename << " read." << endl;
            }

            if (idx_snap == 0)
            {
                t_init = tval;
                for (int i = 0; i < dim; ++i)
                {
                    init_cond->item(i) = sample[i]; 
                }
                dmd->projectInitialCondition(init_cond);
                if (t_final > 0.0) // Practical prediction
                {
                    num_tests = 1;
                    if (myid == 0)
                    {
                        cout << "Predicting DMD solution at t = " << t_final << "." << endl;
                    }
                    dmd_prediction_timer.Start();
                    result = dmd->predict((t_final-t_init)/dtc);
                    dmd_prediction_timer.Stop();
                    // TODO: store result
                    return 0;
                }
            }
            else // Verify DMD prediction results against dataset
            {
                double dmd_power = (tval-t_init)/dtc;
                if (myid == 0)
                {
                    cout << "Predicting DMD solution #" << idx_snap << " at t = " << tval << "." << endl;
                    cout << "DMD power = " << dmd_power << "." << endl;
                }
                dmd_prediction_timer.Start();
                result = dmd->predict(dmd_power);
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

    MFEM_VERIFY(num_tests > 0, "No prediction is made.");

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
