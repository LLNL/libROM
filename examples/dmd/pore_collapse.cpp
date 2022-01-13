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
#include "utils/HDFDatabase.h"
#include <cmath>
#include <fstream>
#include <iostream>

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
    double ddt = 0.00005;
    double dmd_epsilon = 1.0;
    double ef = 1.0;
    int rdim = -1;
    const char *list_dir = "/usr/workspace/nlrom/poreCollapse/libROM_data/pore_collapse_list";
    const char *data_dir = "/usr/workspace/nlrom/poreCollapse/libROM_data/pore_collapse_data";
    const char *var_name = "tkelv";
    int precision = 16;
    cout.precision(precision);

    OptionsParser args(argc, argv);
    args.AddOption(&t_final, "-tf", "--t-final",
                   "Final time; start time is 0.");
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

    if (dtc > 0.0) ddt = dtc;

    std::string variable = std::string(var_name);
    int dim = -1;
    std::ifstream var_dim((std::string(data_dir) + "/dim.txt").c_str());
    var_dim >> dim;
    if (myid == 0)
    {
        cout << "Variable " << var_name << " has dimension " << dim << "." << endl;
    }
    double* sample = new double[dim];

    StopWatch dmd_training_timer, dmd_prediction_timer;


    std::string par_dir; // *gpa 
    std::string snap; // run_036.*

    //CAROM::DMD dmd(dim); // TODO: unify?
    CAROM::AdaptiveDMD dmd(dim, ddt, "LS", "G", dmd_epsilon);

    dmd_training_timer.Start();
    std::ifstream training_par_list((std::string(list_dir) + "/training_gpa").c_str());  
    while (std::getline(training_par_list, par_dir))
    {
        if (myid == 0)
        {
            cout << "Loading samples for " << par_dir << " to train DMD." << endl;
        }
        std::ifstream snap_ifs((std::string(list_dir) + "/" + par_dir).c_str());
        int num_samp = 0;
        double tval = 0.0;
        while (std::getline(snap_ifs, snap))
        {
            std::ifstream tval_ifs((std::string(data_dir) + "/" + par_dir + "/" + snap + "/tval.txt").c_str());
            tval_ifs >> tval;
            std::ifstream data_ifs((std::string(data_dir) + "/" + par_dir + "/" + snap + "/" + variable + ".csv").c_str()); // {zone*, tkelv}.csv
            int row = 0;
            double data = 0.0; 
            while (data_ifs >> data)
            {
                sample[row++] = data;
            }
            MFEM_VERIFY(row == dim, "Dimension disagree.");
            dmd.takeSample(sample, tval);
            num_samp++;
        }
        if (myid == 0)
        {
            cout << "Loaded " << num_samp << " samples for " << par_dir << "." << endl;
        }
    }
    training_par_list.close();

    if (rdim != -1)
    {
        if (myid == 0)
        {
            cout << "Creating DMD with rdim: " << rdim << endl;
        }
        dmd.train(rdim);
    }
    else if (ef != -1)
    {
        if (myid == 0)
        {
            cout << "Creating DMD with energy fraction: " << ef << endl;
        }
        dmd.train(ef);
    }

    dmd_training_timer.Stop();

    int num_tests = 0;
    std::ifstream testing_par_list((std::string(list_dir) + "/testing_gpa").c_str()); 
    CAROM::Vector* init_cond = new CAROM::Vector(dim, true);
    CAROM::Vector* result = new CAROM::Vector(dim, true);
    while (std::getline(testing_par_list, par_dir))
    {
        if (myid == 0)
        {
            cout << "Predicting solution for " << par_dir << " using DMD." << endl;
        }
        std::ifstream snap_list((std::string(list_dir) + "/" + par_dir).c_str());
        std::string data_filename; // {zone_*,tkelv}.csv
        int num_steps = 0;
        double tval = 0.0;
        while (std::getline(snap_list, snap))
        {
            std::ifstream tval_ifs((std::string(data_dir) + "/" + par_dir + "/" + snap + "/tval.txt").c_str());
            tval_ifs >> tval;
            data_filename = std::string(data_dir) + "/" + par_dir + "/" + snap + "/" + variable + ".csv"; 
            std::ifstream data_ifs(data_filename.c_str()); 
            int row = 0;
            double data = 0.0; 
            while (data_ifs >> data)
            {
                if (num_steps == 0)
                {
                    init_cond->item(row++) = data; 
                }
                else 
                {
                    sample[row++] = data; 
                }
            }
            if (myid == 0)
            {
                cout << "State " << data_filename << " read." << endl;
            }
            MFEM_VERIFY(row == dim, "Dimension disagree.");

            if (num_steps == 0)
            {
                dmd.projectInitialCondition(init_cond);
                if (t_final > 0.0)
                {
                    num_tests = 1;
                    if (myid == 0)
                    {
                        cout << "Predicting DMD solution at t = " << t_final << "." << endl;
                    }
                    dmd_prediction_timer.Start();
                    result = dmd.predict(t_final/ddt);
                    dmd_prediction_timer.Stop();
                    // TODO: store result
                    break;
                }
            }
            else // Verify DMD prediction results
            {
                double dmd_power = (dtc > 0) ? num_steps : tval/ddt;
                if (myid == 0)
                {
                    cout << "Predicting DMD solution #" << num_steps << " at t = " << tval << "." << endl;
                    cout << "DMD power = " << dmd_power << "." << endl;
                }
                dmd_prediction_timer.Start();
                result = dmd.predict(dmd_power);
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
                    printf("Total elapsed time for predicting DMD: %e second\n", dmd_prediction_timer.RealTime());
                    cout << "Norm of true solution at t = " << tval << " is " << tot_true_solution_norm << endl;
                    cout << "Absolute error of DMD solution at t = " << tval << " is " << tot_diff_norm << endl;
                    cout << "Relative error of DMD solution at t = " << tval << " is " << tot_diff_norm / tot_true_solution_norm << endl;
                }

            }
            num_steps++;
        }
        num_tests += num_steps;
    }
    testing_par_list.close();
    delete result;

    MFEM_VERIFY(num_tests > 0, "No prediction is made.");

    if (myid == 0)
    {
        printf("Elapsed time for training DMD: %e second\n", dmd_training_timer.RealTime());
        printf("Total elapsed time for predicting DMD: %e second\n", dmd_prediction_timer.RealTime());
        printf("Average elapsed time for predicting DMD: %e second\n", dmd_prediction_timer.RealTime() / num_tests);
    }

    delete sample;
    return 0;
}
