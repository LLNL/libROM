//                         libROM MFEM Example: DG Euler Equations (adapted from ex18p.cpp)
//
// Compile with: make dg_euler
//
// For DMD:
//   mpirun -n 8 dg_euler
//   mpirun -n 8 dg_euler -p 2 -rs 2 -rp 1 -o 1 -s 3 -visit
//
// Sample runs:
//
//       mpirun -np 4 dg_euler -p 1 -rs 2 -rp 1 -o 1 -s 3
//       mpirun -np 4 dg_euler -p 1 -rs 1 -rp 1 -o 3 -s 4
//       mpirun -np 4 dg_euler -p 1 -rs 1 -rp 1 -o 5 -s 6
//       mpirun -np 4 dg_euler -p 2 -rs 1 -rp 1 -o 1 -s 3
//       mpirun -np 4 dg_euler -p 2 -rs 1 -rp 1 -o 3 -s 3
//
// Description:  This example code solves the compressible Euler system of
//               equations, a model nonlinear hyperbolic PDE, with a
//               discontinuous Galerkin (DG) formulation.
//
//               Specifically, it solves for an exact solution of the equations
//               whereby a vortex is transported by a uniform flow. Since all
//               boundaries are periodic here, the method's accuracy can be
//               assessed by measuring the difference between the solution and
//               the initial condition at a later time when the vortex returns
//               to its initial location.
//
//               Note that as the order of the spatial discretization increases,
//               the timestep must become smaller. This example currently uses a
//               simple estimate derived by Cockburn and Shu for the 1D RKDG
//               method. An additional factor can be tuned by passing the --cfl
//               (or -c shorter) flag.
//
//               The example demonstrates user-defined bilinear and nonlinear
//               form integrators for systems of equations that are defined with
//               block vectors, and how these are used with an operator for
//               explicit time integrators. In this case the system also
//               involves an external approximate Riemann solver for the DG
//               interface flux. It also demonstrates how to use GLVis for
//               in-situ visualization of vector grid functions.
//
//               We recommend viewing examples 9, 14 and 17 before viewing this
//               example.

#include "mfem.hpp"
#include "DMD.h"
#include "Vector.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

// Classes FE_Evolution, RiemannSolver, DomainIntegrator and FaceIntegrator
// shared between the serial and parallel version of the example.
#include "dg_euler.hpp"

// Choice for the problem setup. See InitialCondition in dg_euler.hpp.
int problem;

// Equation constant parameters.
const int num_equation = 4;
const double specific_heat_ratio = 1.4;
const double gas_constant = 1.0;

// Maximum characteristic speed (updated by integrators)
double max_char_speed;

int main(int argc, char *argv[])
{
    // 1. Initialize MPI.
    MPI_Session mpi(argc, argv);

    // 2. Parse command-line options.
    problem = 1;
    const char *mesh_file = "../dependencies/mfem/data/periodic-square.mesh";
    int ser_ref_levels = 0;
    int par_ref_levels = 1;
    int order = 3;
    int ode_solver_type = 4;
    double t_final = 2.0;
    double dt = -0.01;
    double cfl = 0.3;
    double ef = 0.9999;
    int rdim = -1;
    bool visualization = true;
    bool visit = false;
    int vis_steps = 50;

    int precision = 8;
    cout.precision(precision);

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&problem, "-p", "--problem",
                   "Problem setup to use. See options in velocity_function().");
    args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                   "Number of times to refine the mesh uniformly before parallel"
                   " partitioning, -1 for auto.");
    args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                   "Number of times to refine the mesh uniformly after parallel"
                   " partitioning.");
    args.AddOption(&order, "-o", "--order",
                   "Order (degree) of the finite elements.");
    args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                   "ODE solver: 1 - Forward Euler,\n\t"
                   "            2 - RK2 SSP, 3 - RK3 SSP, 4 - RK4, 6 - RK6.");
    args.AddOption(&t_final, "-tf", "--t-final",
                   "Final time; start time is 0.");
    args.AddOption(&dt, "-dt", "--time-step",
                   "Time step. Positive number skips CFL timestep calculation.");
    args.AddOption(&cfl, "-c", "--cfl-number",
                   "CFL number for timestep calculation.");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                   "--no-visualization",
                   "Enable or disable GLVis visualization.");
    args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                   "--no-visit-datafiles",
                   "Save data files for VisIt (visit.llnl.gov) visualization.");
    args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                   "Visualize every n-th timestep.");
    args.AddOption(&ef, "-ef", "--energy_fraction",
                   "Energy fraction for DMD.");
    args.AddOption(&rdim, "-rdim", "--rdim",
                   "Reduced dimension for DMD.");

    args.Parse();
    if (!args.Good())
    {
        if (mpi.Root()) {
            args.PrintUsage(cout);
        }
        return 1;
    }
    if (mpi.Root()) {
        args.PrintOptions(cout);
    }

    // 3. Read the mesh from the given mesh file. This example requires a 2D
    //    periodic mesh, such as ../data/periodic-square.mesh.
    Mesh mesh(mesh_file, 1, 1);
    const int dim = mesh.Dimension();

    MFEM_ASSERT(dim == 2, "Need a two-dimensional mesh for the problem definition");

    // 4. Define the ODE solver used for time integration. Several explicit
    //    Runge-Kutta methods are available.
    ODESolver *ode_solver = NULL;
    switch (ode_solver_type)
    {
    case 1:
        ode_solver = new ForwardEulerSolver;
        break;
    case 2:
        ode_solver = new RK2Solver(1.0);
        break;
    case 3:
        ode_solver = new RK3SSPSolver;
        break;
    case 4:
        ode_solver = new RK4Solver;
        break;
    case 6:
        ode_solver = new RK6Solver;
        break;
    default:
        if (mpi.Root())
        {
            cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
        }
        return 3;
    }

    // 5. Refine the mesh in serial to increase the resolution. In this example
    //    we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
    //    a command-line parameter.
    for (int lev = 0; lev < ser_ref_levels; lev++)
    {
        mesh.UniformRefinement();
    }

    // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
    //    this mesh further in parallel to increase the resolution. Once the
    //    parallel mesh is defined, the serial mesh can be deleted.
    ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();
    for (int lev = 0; lev < par_ref_levels; lev++)
    {
        pmesh.UniformRefinement();
    }

    // 7. Define the discontinuous DG finite element space of the given
    //    polynomial order on the refined mesh.
    DG_FECollection fec(order, dim);
    // Finite element space for a scalar (thermodynamic quantity)
    ParFiniteElementSpace fes(&pmesh, &fec);
    // Finite element space for a mesh-dim vector quantity (momentum)
    ParFiniteElementSpace dfes(&pmesh, &fec, dim, Ordering::byNODES);
    // Finite element space for all variables together (total thermodynamic state)
    ParFiniteElementSpace vfes(&pmesh, &fec, num_equation, Ordering::byNODES);

    // This example depends on this ordering of the space.
    MFEM_ASSERT(fes.GetOrdering() == Ordering::byNODES, "");

    HYPRE_BigInt glob_size = vfes.GlobalTrueVSize();
    if (mpi.Root()) {
        cout << "Number of unknowns: " << glob_size << endl;
    }

    // 8. Define the initial conditions, save the corresponding mesh and grid
    //    functions to a file. This can be opened with GLVis with the -gc option.

    // The solution u has components {density, x-momentum, y-momentum, energy}.
    // These are stored contiguously in the BlockVector u_block.
    Array<int> offsets(num_equation + 1);
    for (int k = 0; k <= num_equation; k++) {
        offsets[k] = k * vfes.GetNDofs();
    }
    BlockVector u_block(offsets);

    // Momentum grid function on dfes for visualization.
    ParGridFunction mom(&dfes, u_block.GetData() + offsets[1]);

    // Initialize the state.
    VectorFunctionCoefficient u0(num_equation, InitialCondition);
    ParGridFunction sol(&vfes, u_block.GetData());
    sol.ProjectCoefficient(u0);

    // Output the initial solution.
    {
        ostringstream mesh_name;
        mesh_name << "vortex-mesh." << setfill('0') << setw(6) << mpi.WorldRank();
        ofstream mesh_ofs(mesh_name.str().c_str());
        mesh_ofs.precision(precision);
        mesh_ofs << pmesh;

        for (int k = 0; k < num_equation; k++)
        {
            ParGridFunction uk(&fes, u_block.GetBlock(k));
            ostringstream sol_name;
            sol_name << "vortex-" << k << "-init."
                     << setfill('0') << setw(6) << mpi.WorldRank();
            ofstream sol_ofs(sol_name.str().c_str());
            sol_ofs.precision(precision);
            sol_ofs << uk;
        }
    }

    // 9. Set up the nonlinear form corresponding to the DG discretization of the
    //    flux divergence, and assemble the corresponding mass matrix.
    MixedBilinearForm Aflux(&dfes, &fes);
    Aflux.AddDomainIntegrator(new DomainIntegrator(dim));
    Aflux.Assemble();

    ParNonlinearForm A(&vfes);
    RiemannSolver rsolver;
    A.AddInteriorFaceIntegrator(new FaceIntegrator(rsolver, dim));

    // 10. Define the time-dependent evolution operator describing the ODE
    //     right-hand side, and perform time-integration (looping over the time
    //     iterations, ti, with a time-step dt).
    FE_Evolution euler(vfes, A, Aflux.SpMat());

    // Visualize the density
    socketstream sout;
    if (visualization)
    {
        char vishost[] = "localhost";
        int  visport   = 19916;

        MPI_Barrier(pmesh.GetComm());
        sout.open(vishost, visport);
        if (!sout)
        {
            if (mpi.Root())
            {
                cout << "Unable to connect to GLVis server at "
                     << vishost << ':' << visport << endl;
            }
            visualization = false;
            if (mpi.Root()) {
                cout << "GLVis visualization disabled.\n";
            }
        }
        else
        {
            sout << "parallel " << mpi.WorldSize() << " " << mpi.WorldRank() << "\n";
            sout.precision(precision);
            sout << "solution\n" << pmesh << mom;
            sout << "pause\n";
            sout << flush;
            if (mpi.Root())
            {
                cout << "GLVis visualization paused."
                     << " Press space (in the GLVis window) to resume it.\n";
            }
        }
    }

    VisItDataCollection visit_dc("DG_Euler", &pmesh);
    visit_dc.RegisterField("solution", &mom);
    if (visit)
    {
        visit_dc.SetCycle(0);
        visit_dc.SetTime(0.0);
        visit_dc.Save();
    }


    // Determine the minimum element size.
    double hmin;
    if (cfl > 0)
    {
        double my_hmin = pmesh.GetElementSize(0, 1);
        for (int i = 1; i < pmesh.GetNE(); i++)
        {
            my_hmin = min(pmesh.GetElementSize(i, 1), my_hmin);
        }
        // Reduce to find the global minimum element size
        MPI_Allreduce(&my_hmin, &hmin, 1, MPI_DOUBLE, MPI_MIN, pmesh.GetComm());
    }

    // Start the timer.
    tic_toc.Clear();
    tic_toc.Start();

    StopWatch fom_timer, dmd_training_timer, dmd_prediction_timer;

    fom_timer.Start();

    double t = 0.0;
    euler.SetTime(t);
    ode_solver->Init(euler);

    fom_timer.Stop();

    dmd_training_timer.Start();

    CAROM::DMD dmd_dens(u_block.GetBlock(0).Size());
    CAROM::DMD dmd_x_mom(u_block.GetBlock(1).Size());
    CAROM::DMD dmd_y_mom(u_block.GetBlock(2).Size());
    CAROM::DMD dmd_e(u_block.GetBlock(3).Size());
    dmd_dens.takeSample(u_block.GetBlock(0).GetData());
    dmd_x_mom.takeSample(u_block.GetBlock(1).GetData());
    dmd_y_mom.takeSample(u_block.GetBlock(2).GetData());
    dmd_e.takeSample(u_block.GetBlock(3).GetData());

    dmd_training_timer.Stop();

    if (cfl > 0)
    {
        // Find a safe dt, using a temporary vector. Calling Mult() computes the
        // maximum char speed at all quadrature points on all faces.
        max_char_speed = 0.;
        Vector z(sol.Size());
        A.Mult(sol, z);
        // Reduce to find the global maximum wave speed
        {
            double all_max_char_speed;
            MPI_Allreduce(&max_char_speed, &all_max_char_speed,
                          1, MPI_DOUBLE, MPI_MAX, pmesh.GetComm());
            max_char_speed = all_max_char_speed;
        }
        dt = cfl * hmin / max_char_speed / (2*order+1);
    }

    // Integrate in time.
    bool done = false;
    for (int ti = 0; !done; )
    {
        fom_timer.Start();

        double dt_real = min(dt, t_final - t);

        ode_solver->Step(sol, t, dt_real);
        if (cfl > 0)
        {
            // Reduce to find the global maximum wave speed
            {
                double all_max_char_speed;
                MPI_Allreduce(&max_char_speed, &all_max_char_speed,
                              1, MPI_DOUBLE, MPI_MAX, pmesh.GetComm());
                max_char_speed = all_max_char_speed;
            }
            dt = cfl * hmin / max_char_speed / (2*order+1);
        }
        ti++;

        done = (t >= t_final - 1e-8*dt);

        fom_timer.Stop();

        dmd_training_timer.Start();

        dmd_dens.takeSample(u_block.GetBlock(0).GetData());
        dmd_x_mom.takeSample(u_block.GetBlock(1).GetData());
        dmd_y_mom.takeSample(u_block.GetBlock(2).GetData());
        dmd_e.takeSample(u_block.GetBlock(3).GetData());

        dmd_training_timer.Stop();

        if (done || ti % vis_steps == 0)
        {
            if (mpi.Root())
            {
                cout << "time step: " << ti << ", time: " << t << endl;
            }
            if (visualization)
            {
                MPI_Barrier(pmesh.GetComm());
                sout << "parallel " << mpi.WorldSize() << " " << mpi.WorldRank() << "\n";
                sout << "solution\n" << pmesh << mom << flush;
            }
            if (visit)
            {
                visit_dc.SetCycle(ti);
                visit_dc.SetTime(t);
                visit_dc.Save();
            }
        }
    }

    tic_toc.Stop();
    if (mpi.Root()) {
        cout << " done, " << tic_toc.RealTime() << "s." << endl;
    }

    // 11. Save the final solution. This output can be viewed later using GLVis:
    //     "glvis -np 4 -m vortex-mesh -g vortex-1-final".
    for (int k = 0; k < num_equation; k++)
    {
        ParGridFunction uk(&fes, u_block.GetBlock(k));
        ostringstream sol_name;
        sol_name << "vortex-" << k << "-final."
                 << setfill('0') << setw(6) << mpi.WorldRank();
        ofstream sol_ofs(sol_name.str().c_str());
        sol_ofs.precision(precision);
        sol_ofs << uk;
    }

    // 12. Compute the L2 solution error summed for all components.
    if (t_final == 2.0)
    {
        const double error = sol.ComputeLpError(2, u0);
        if (mpi.Root()) {
            cout << "Solution error: " << error << endl;
        }
    }

    // 13. Calculate the DMD modes.
    if (mpi.WorldRank() == 0 && rdim != -1 && ef != -1)
    {
        std::cout << "Both rdim and ef are set. ef will be ignored." << std::endl;
    }

    dmd_training_timer.Start();

    if (rdim != -1)
    {
        if (mpi.WorldRank() == 0)
        {
            std::cout << "Creating DMD with rdim: " << rdim << std::endl;
        }
        dmd_dens.train(rdim);
        dmd_x_mom.train(rdim);
        dmd_y_mom.train(rdim);
        dmd_e.train(rdim);
    }
    else if (ef != -1)
    {
        if (mpi.WorldRank() == 0)
        {
            std::cout << "Creating DMD with energy fraction: " << ef << std::endl;
        }
        dmd_dens.train(ef);
        dmd_x_mom.train(ef);
        dmd_y_mom.train(ef);
        dmd_e.train(ef);
    }

    dmd_training_timer.Stop();

    dmd_prediction_timer.Start();

    // 14. Predict the state at t_final using DMD.
    if (mpi.WorldRank() == 0)
    {
        std::cout << "Predicting density, momentum, and energy at t_final using DMD" << std::endl;
    }
    CAROM::Vector* result_dens = dmd_dens.predict(t_final/dt);
    CAROM::Vector* result_x_mom = dmd_x_mom.predict(t_final/dt);
    CAROM::Vector* result_y_mom = dmd_y_mom.predict(t_final/dt);
    CAROM::Vector* result_e = dmd_e.predict(t_final/dt);

    dmd_prediction_timer.Stop();

    // 15. Calculate the relative error between the DMD final solution and the true solution.
    Vector dmd_solution_dens(result_dens->getData(), result_dens->dim());
    Vector true_solution_dens(u_block.GetBlock(0).GetData(), u_block.GetBlock(0).Size());
    Vector diff_dens(true_solution_dens.Size());
    subtract(dmd_solution_dens, true_solution_dens, diff_dens);

    double tot_diff_norm_dens = sqrt(InnerProduct(MPI_COMM_WORLD, diff_dens, diff_dens));
    double tot_true_solution_dens_norm = sqrt(InnerProduct(MPI_COMM_WORLD, true_solution_dens, true_solution_dens));

    Vector dmd_solution_x_mom(result_x_mom->getData(), result_x_mom->dim());
    Vector true_solution_x_mom(u_block.GetBlock(1).GetData(), u_block.GetBlock(1).Size());
    Vector diff_x_mom(true_solution_x_mom.Size());
    subtract(dmd_solution_x_mom, true_solution_x_mom, diff_x_mom);

    double tot_diff_norm_x_mom = sqrt(InnerProduct(MPI_COMM_WORLD, diff_x_mom, diff_x_mom));
    double tot_true_solution_x_mom_norm = sqrt(InnerProduct(MPI_COMM_WORLD, true_solution_x_mom, true_solution_x_mom));

    Vector dmd_solution_y_mom(result_y_mom->getData(), result_y_mom->dim());
    Vector true_solution_y_mom(u_block.GetBlock(2).GetData(), u_block.GetBlock(2).Size());
    Vector diff_y_mom(true_solution_y_mom.Size());
    subtract(dmd_solution_y_mom, true_solution_y_mom, diff_y_mom);

    double tot_diff_norm_y_mom = sqrt(InnerProduct(MPI_COMM_WORLD, diff_y_mom, diff_y_mom));
    double tot_true_solution_y_mom_norm = sqrt(InnerProduct(MPI_COMM_WORLD, true_solution_y_mom, true_solution_y_mom));

    Vector dmd_solution_e(result_e->getData(), result_e->dim());
    Vector true_solution_e(u_block.GetBlock(3).GetData(), u_block.GetBlock(3).Size());
    Vector diff_e(true_solution_e.Size());
    subtract(dmd_solution_e, true_solution_e, diff_e);

    double tot_diff_norm_e = sqrt(InnerProduct(MPI_COMM_WORLD, diff_e, diff_e));
    double tot_true_solution_e_norm = sqrt(InnerProduct(MPI_COMM_WORLD, true_solution_e, true_solution_e));

    if (mpi.WorldRank() == 0)
    {
        std::cout << "Relative error of DMD density (dens) at t_final: " << t_final << " is " << tot_diff_norm_dens / tot_true_solution_dens_norm << std::endl;
        std::cout << "Relative error of DMD x-momentum (x_mom) at t_final: " << t_final << " is " << tot_diff_norm_x_mom / tot_true_solution_x_mom_norm << std::endl;
        std::cout << "Relative error of DMD y-momentum (y_mom) at t_final: " << t_final << " is " << tot_diff_norm_y_mom / tot_true_solution_y_mom_norm << std::endl;
        std::cout << "Relative error of DMD energy (e) at t_final: " << t_final << " is " << tot_diff_norm_e / tot_true_solution_e_norm << std::endl;
        printf("Elapsed time for solving FOM: %e second\n", fom_timer.RealTime());
        printf("Elapsed time for training DMD: %e second\n", dmd_training_timer.RealTime());
        printf("Elapsed time for predicting DMD: %e second\n", dmd_prediction_timer.RealTime());
    }

    // Free the used memory.
    delete ode_solver;
    delete result_dens;
    delete result_x_mom;
    delete result_y_mom;
    delete result_e;

    return 0;
}
