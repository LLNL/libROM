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

    double t_final = 0.006;
    double dt = 0.0001;
    double ef = 1.0;
    int rdim = -1;
    const char *data_dir = "/usr/workspace/nlrom/poreCollapse";
    const char *list_dir = "pore_collapse_list";
    const char *var_name = "tkelv";
    int precision = 16;
    cout.precision(precision);

    OptionsParser args(argc, argv);
    args.AddOption(&t_final, "-tf", "--t-final",
                   "Final time; start time is 0.");
    args.AddOption(&dt, "-dt", "--time-step",
                   "Time step size.");
    args.AddOption(&ef, "-ef", "--energy_fraction",
                   "Energy fraction for DMD.");
    args.AddOption(&rdim, "-rdim", "--rdim",
                   "Reduced dimension for DMD.");
    args.AddOption(&data_dir, "-data", "--data_dir",
                   "Location of data.");
    args.AddOption(&list_dir, "-list", "--list_dir",
                   "Location of training and testing data list.");
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

    std::string variable = std::string(var_name);
    int dim = 1; // TODO: read
    CAROM::DMD dmd(dim);
    double* sample = new double[dim];

    std::string par_dir; 
    std::ifstream training_par_list((std::string(list_dir) + "/training_par").c_str()); 
    while std::getline(training_par_list, par_dir)
    {
        std::ifstream snap_list((std::string(list_dir) + "/" + par_dir).c_str());
        std::string snap;
        while std::getline(snap_list, snap)
        {
            CAROM::HDFDatabase database;
            database.open(std::string(data_dir) + "/" + par_dir "/" + snap, "r");
            database.getDoubleArray(variable, sample, 1);
            dmd.takeSample(sample);
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

    std::ifstream testing_par_list((std::string(list_dir) + "/testing_par").c_str()); 
    while std::getline(testing_par_list, par_dir)
    {
        if (myid == 0)
        {
            cout << "Predicting solution for " << par_dir << " using DMD" << endl;
        }
        std::ifstream snap_list((std::string(list_dir) + "/" + par_dir).c_str());
        std::string snap;
        int num_steps = 0;
        dmd.clearSample();
        while std::getline(snap_list, snap)
        {
            CAROM::HDFDatabase database;
            database.open(std::string(data_dir) + "/" + par_dir "/" + snap, "r");
            database.getDoubleArray(variable, sample, 1);
            dmd.takeSample(sample);
            num_steps += 1;
        }

        // 14. Predict the state at t_final using DMD.
        num_steps = std::min(num_steps, round(t_final/dt));
        t_final = dt*num_steps;
        if (myid == 0)
        {
            cout << "Final time t_final = " << t_final << endl;
        }
        CAROM::Vector* result = dmd.predict(num_steps - 1);

        // 15. Calculate the relative error between the DMD final solution and the true solution.
        Vector dmd_solution(result->getData(), result->dim());
        Vector true_solution(sample, num_steps);
        Vector diff(true_solution.Size());
        subtract(dmd_solution, true_solution, diff);

        double tot_diff_norm = sqrt(InnerProduct(MPI_COMM_WORLD, diff, diff));
        double tot_true_solution_norm = sqrt(InnerProduct(MPI_COMM_WORLD, true_solution, true_solution));

        if (myid == 0)
        {
            cout << "Relative error of DMD solution at t_final is " << tot_diff_norm / tot_true_solution_norm << endl;
        }

        delete result;
    }
    testing_par_list.close();

    delete sample;
    return 0;
}
