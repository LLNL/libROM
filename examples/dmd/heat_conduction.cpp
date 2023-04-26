/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

//                       libROM MFEM Example: Heat_Conduction (adapted from ex16p.cpp)
//
// Compile with: make heat_conduction
//
// =================================================================================
//
// Sample runs and results for DMD:
//
// Command 1:
//   mpirun -np 8 heat_conduction -s 1 -a 0.0 -k 1.0 -visit
//
// Output 1:
//   Relative error of DMD temperature (u) at t_final: 0.5 is 0.00049906934
//
// Command 2:
//   mpirun -np 8 heat_conduction -s 3 -a 0.5 -k 0.5 -o 4 -tf 0.7 -vs 1 -visit
//
// Output 2:
//   Relative error of DMD temperature (u) at t_final: 0.7 is 0.00082289823
//
// =================================================================================
//
// Sample runs and results for NonuniformDMD:
//
// Command 1:
//   mpirun -np 8 heat_conduction -s 1 -a 0.0 -k 1.0 -nonunif -visit
//
// Output 1:
//   Relative error of DMD temperature (u) at t_final: 0.5 is 0.00071958947
//
// Command 2:
//   mpirun -np 8 heat_conduction -s 3 -a 0.5 -k 0.5 -o 4 -tf 0.7 -vs 1 -nonunif -visit
//
// Output 2:
//   Relative error of DMD temperature (u) at t_final: 0.7 is 0.00013450754
//
// =================================================================================
//
// Description:  This example solves a time dependent nonlinear heat equation
//               problem of the form du/dt = C(u), with a non-linear diffusion
//               operator C(u) = \nabla \cdot (\kappa + \alpha u) \nabla u.
//
//               The example demonstrates the use of nonlinear operators (the
//               class ConductionOperator defining C(u)), as well as their
//               implicit time integration. Note that implementing the method
//               ConductionOperator::ImplicitSolve is the only requirement for
//               high-order implicit (SDIRK) time integration. Optional saving
//               with ADIOS2 (adios2.readthedocs.io) is also illustrated.

#include "mfem.hpp"
#include "algo/DMD.h"
#include "algo/NonuniformDMD.h"
#include "linalg/Vector.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

/** After spatial discretization, the conduction model can be written as:
 *
 *     du/dt = M^{-1}(-Ku)
 *
 *  where u is the vector representing the temperature, M is the mass matrix,
 *  and K is the diffusion operator with diffusivity depending on u:
 *  (\kappa + \alpha u).
 *
 *  Class ConductionOperator represents the right-hand side of the above ODE.
 */
class ConductionOperator : public TimeDependentOperator
{
protected:
    ParFiniteElementSpace &fespace;
    Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.

    ParBilinearForm *M;
    ParBilinearForm *K;

    HypreParMatrix Mmat;
    HypreParMatrix Kmat;
    HypreParMatrix *T; // T = M + dt K
    double current_dt;

    CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
    HypreSmoother M_prec; // Preconditioner for the mass matrix M

    CGSolver T_solver;    // Implicit solver for T = M + dt K
    HypreSmoother T_prec; // Preconditioner for the implicit solver

    double alpha, kappa;

    mutable Vector z; // auxiliary vector

public:
    ConductionOperator(ParFiniteElementSpace &f, double alpha, double kappa,
                       const Vector &u);

    virtual void Mult(const Vector &u, Vector &du_dt) const;
    /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
        This is the only requirement for high-order SDIRK implicit integration.*/
    virtual void ImplicitSolve(const double dt, const Vector &u, Vector &k);

    /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
    void SetParameters(const Vector &u);

    virtual ~ConductionOperator();
};

double InitialTemperature(const Vector &x);

int main(int argc, char *argv[])
{
    // 1. Initialize MPI.
    int num_procs, myid;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // 2. Parse command-line options.
    const char *mesh_file = "../data/star.mesh";
    int ser_ref_levels = 2;
    int par_ref_levels = 1;
    int order = 2;
    int ode_solver_type = 3;
    double t_final = 0.5;
    double dt = 1.0e-2;
    double alpha = 1.0e-2;
    double kappa = 0.5;
    double ef = 0.9999;
    int rdim = -1;
    bool visualization = true;
    bool visit = false;
    int vis_steps = 5;
    bool use_nonuniform = false;
    bool adios2 = false;

    int precision = 8;
    cout.precision(precision);

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                   "Number of times to refine the mesh uniformly in serial.");
    args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                   "Number of times to refine the mesh uniformly in parallel.");
    args.AddOption(&order, "-o", "--order",
                   "Order (degree) of the finite elements.");
    args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                   "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3,\n\t"
                   "\t   11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.");
    args.AddOption(&t_final, "-tf", "--t-final",
                   "Final time; start time is 0.");
    args.AddOption(&dt, "-dt", "--time-step",
                   "Time step.");
    args.AddOption(&alpha, "-a", "--alpha",
                   "Alpha coefficient.");
    args.AddOption(&kappa, "-k", "--kappa",
                   "Kappa coefficient offset.");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                   "--no-visualization",
                   "Enable or disable GLVis visualization.");
    args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                   "--no-visit-datafiles",
                   "Save data files for VisIt (visit.llnl.gov) visualization.");
    args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                   "Visualize every n-th timestep.");
    args.AddOption(&adios2, "-adios2", "--adios2-streams", "-no-adios2",
                   "--no-adios2-streams",
                   "Save data using adios2 streams.");
    args.AddOption(&ef, "-ef", "--energy_fraction",
                   "Energy fraction for DMD.");
    args.AddOption(&rdim, "-rdim", "--rdim",
                   "Reduced dimension for DMD.");
    args.AddOption(&use_nonuniform, "-nonunif", "--nonunif", "-no-nonunif",
                   "--no-nonunif",
                   "Use NonuniformDMD.");
    args.Parse();
    if (!args.Good())
    {
        args.PrintUsage(cout);
        MPI_Finalize();
        return 1;
    }

    if (myid == 0)
    {
        args.PrintOptions(cout);
    }

    // 3. Read the serial mesh from the given mesh file on all processors. We can
    //    handle triangular, quadrilateral, tetrahedral and hexahedral meshes
    //    with the same code.
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    // 4. Define the ODE solver used for time integration. Several implicit
    //    singly diagonal implicit Runge-Kutta (SDIRK) methods, as well as
    //    explicit Runge-Kutta methods are available.
    ODESolver *ode_solver;
    switch (ode_solver_type)
    {
    // Implicit L-stable methods
    case 1:
        ode_solver = new BackwardEulerSolver;
        break;
    case 2:
        ode_solver = new SDIRK23Solver(2);
        break;
    case 3:
        ode_solver = new SDIRK33Solver;
        break;
    // Explicit methods
    case 11:
        ode_solver = new ForwardEulerSolver;
        break;
    case 12:
        ode_solver = new RK2Solver(0.5);
        break; // midpoint method
    case 13:
        ode_solver = new RK3SSPSolver;
        break;
    case 14:
        ode_solver = new RK4Solver;
        break;
    case 15:
        ode_solver = new GeneralizedAlphaSolver(0.5);
        break;
    // Implicit A-stable methods (not L-stable)
    case 22:
        ode_solver = new ImplicitMidpointSolver;
        break;
    case 23:
        ode_solver = new SDIRK23Solver;
        break;
    case 24:
        ode_solver = new SDIRK34Solver;
        break;
    default:
        cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
        delete mesh;
        return 3;
    }

    // 5. Refine the mesh in serial to increase the resolution. In this example
    //    we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
    //    a command-line parameter.
    for (int lev = 0; lev < ser_ref_levels; lev++)
    {
        mesh->UniformRefinement();
    }

    // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
    //    this mesh further in parallel to increase the resolution. Once the
    //    parallel mesh is defined, the serial mesh can be deleted.
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
    for (int lev = 0; lev < par_ref_levels; lev++)
    {
        pmesh->UniformRefinement();
    }

    // 7. Define the vector finite element space representing the current and the
    //    initial temperature, u_ref.
    H1_FECollection fe_coll(order, dim);
    ParFiniteElementSpace fespace(pmesh, &fe_coll);

    int fe_size = fespace.GlobalTrueVSize();
    if (myid == 0)
    {
        cout << "Number of temperature unknowns: " << fe_size << endl;
    }

    ParGridFunction u_gf(&fespace);

    // 8. Set the initial conditions for u. All boundaries are considered
    //    natural.
    FunctionCoefficient u_0(InitialTemperature);
    u_gf.ProjectCoefficient(u_0);
    Vector u;
    u_gf.GetTrueDofs(u);

    // 9. Initialize the conduction operator and the VisIt visualization.
    ConductionOperator oper(fespace, alpha, kappa, u);

    u_gf.SetFromTrueDofs(u);
    {
        ostringstream mesh_name, sol_name;
        mesh_name << "heat_conduction-mesh." << setfill('0') << setw(6) << myid;
        sol_name << "heat_conduction-init." << setfill('0') << setw(6) << myid;
        ofstream omesh(mesh_name.str().c_str());
        omesh.precision(precision);
        pmesh->Print(omesh);
        ofstream osol(sol_name.str().c_str());
        osol.precision(precision);
        u_gf.Save(osol);
    }

    VisItDataCollection visit_dc("Heat_Conduction", pmesh);
    visit_dc.RegisterField("temperature", &u_gf);
    if (visit)
    {
        visit_dc.SetCycle(0);
        visit_dc.SetTime(0.0);
        visit_dc.Save();
    }

    // Optionally output a BP (binary pack) file using ADIOS2. This can be
    // visualized with the ParaView VTX reader.
#ifdef MFEM_USE_ADIOS2
    ADIOS2DataCollection* adios2_dc = NULL;
    if (adios2)
    {
        std::string postfix(mesh_file);
        postfix.erase(0, std::string("../data/").size() );
        postfix += "_o" + std::to_string(order);
        postfix += "_solver" + std::to_string(ode_solver_type);
        const std::string collection_name = "heat_conduction-p-" + postfix + ".bp";

        adios2_dc = new ADIOS2DataCollection(MPI_COMM_WORLD, collection_name, pmesh);
        adios2_dc->SetParameter("SubStreams", std::to_string(num_procs/2) );
        adios2_dc->RegisterField("temperature", &u_gf);
        adios2_dc->SetCycle(0);
        adios2_dc->SetTime(0.0);
        adios2_dc->Save();
    }
#endif

    socketstream sout;
    if (visualization)
    {
        char vishost[] = "localhost";
        int  visport   = 19916;
        sout.open(vishost, visport);
        sout << "parallel " << num_procs << " " << myid << endl;
        int good = sout.good(), all_good;
        MPI_Allreduce(&good, &all_good, 1, MPI_INT, MPI_MIN, pmesh->GetComm());
        if (!all_good)
        {
            sout.close();
            visualization = false;
            if (myid == 0)
            {
                cout << "Unable to connect to GLVis server at "
                     << vishost << ':' << visport << endl;
                cout << "GLVis visualization disabled.\n";
            }
        }
        else
        {
            sout.precision(precision);
            sout << "solution\n" << *pmesh << u_gf;
            sout << flush;
        }
    }

    StopWatch fom_timer, dmd_training_timer, dmd_prediction_timer;

    fom_timer.Start();

    // 10. Perform time-integration (looping over the time iterations, ti, with a
    //     time-step dt).
    ode_solver->Init(oper);
    double t = 0.0;
    vector<double> ts;

    fom_timer.Stop();

    dmd_training_timer.Start();

    // 11. Create DMD object and take initial sample.
    u_gf.SetFromTrueDofs(u);

    CAROM::DMD* dmd_u;
    if (use_nonuniform)
    {
        dmd_u = new CAROM::NonuniformDMD(u.Size());
    }
    else
    {
        dmd_u = new CAROM::DMD(u.Size(), dt);
    }
    dmd_u->takeSample(u.GetData(), t);
    ts.push_back(t);

    dmd_training_timer.Stop();

    bool last_step = false;
    for (int ti = 1; !last_step; ti++)
    {
        fom_timer.Start();

        if (t + dt >= t_final - dt/2)
        {
            last_step = true;
        }

        ode_solver->Step(u, t, dt);

        fom_timer.Stop();

        dmd_training_timer.Start();

        u_gf.SetFromTrueDofs(u);
        dmd_u->takeSample(u.GetData(), t);
        ts.push_back(t);

        dmd_training_timer.Stop();

        if (last_step || (ti % vis_steps) == 0)
        {
            if (myid == 0)
            {
                cout << "step " << ti << ", t = " << t << endl;
            }

            u_gf.SetFromTrueDofs(u);
            if (visualization)
            {
                sout << "parallel " << num_procs << " " << myid << "\n";
                sout << "solution\n" << *pmesh << u_gf << flush;
            }

            if (visit)
            {
                visit_dc.SetCycle(ti);
                visit_dc.SetTime(t);
                visit_dc.Save();
            }

#ifdef MFEM_USE_ADIOS2
            if (adios2)
            {
                adios2_dc->SetCycle(ti);
                adios2_dc->SetTime(t);
                adios2_dc->Save();
            }
#endif
        }
        oper.SetParameters(u);
    }

#ifdef MFEM_USE_ADIOS2
    if (adios2)
    {
        delete adios2_dc;
    }
#endif

    // 12. Save the final solution in parallel. This output can be viewed later
    //     using GLVis: "glvis -np <np> -m heat_conduction-mesh -g heat_conduction-final".
    {
        ostringstream sol_name;
        sol_name << "heat_conduction-final." << setfill('0') << setw(6) << myid;
        ofstream osol(sol_name.str().c_str());
        osol.precision(precision);
        u_gf.Save(osol);
    }

    // 13. Calculate the DMD modes.
    if (myid == 0 && rdim != -1 && ef != -1)
    {
        std::cout << "Both rdim and ef are set. ef will be ignored." << std::endl;
    }

    dmd_training_timer.Start();

    if (rdim != -1)
    {
        if (myid == 0)
        {
            std::cout << "Creating DMD with rdim: " << rdim << std::endl;
        }
        dmd_u->train(rdim);
    }
    else if (ef != -1)
    {
        if (myid == 0)
        {
            std::cout << "Creating DMD with energy fraction: " << ef << std::endl;
        }
        dmd_u->train(ef);
    }

    dmd_training_timer.Stop();

    Vector true_solution_u(u.Size());
    true_solution_u = u.GetData();

    dmd_prediction_timer.Start();

    // 14. Predict the state at t_final using DMD.
    if (myid == 0)
    {
        std::cout << "Predicting temperature using DMD" << std::endl;
    }

    CAROM::Vector* result_u = dmd_u->predict(ts[0]);
    Vector initial_dmd_solution_u(result_u->getData(), result_u->dim());
    u_gf.SetFromTrueDofs(initial_dmd_solution_u);

    VisItDataCollection dmd_visit_dc("DMD_Heat_Conduction", pmesh);
    dmd_visit_dc.RegisterField("temperature", &u_gf);
    if (visit)
    {
        dmd_visit_dc.SetCycle(0);
        dmd_visit_dc.SetTime(0.0);
        dmd_visit_dc.Save();
    }

    delete result_u;

    if (visit)
    {
        for (int i = 1; i < ts.size(); i++)
        {
            if (i == ts.size() - 1 || (i % vis_steps) == 0)
            {
                result_u = dmd_u->predict(ts[i]);
                Vector dmd_solution_u(result_u->getData(), result_u->dim());
                u_gf.SetFromTrueDofs(dmd_solution_u);

                dmd_visit_dc.SetCycle(i);
                dmd_visit_dc.SetTime(ts[i]);
                dmd_visit_dc.Save();

                delete result_u;
            }
        }
    }

    dmd_prediction_timer.Stop();

    result_u = dmd_u->predict(t_final);

    // 15. Calculate the relative error between the DMD final solution and the true solution.
    Vector dmd_solution_u(result_u->getData(), result_u->dim());
    Vector diff_u(true_solution_u.Size());
    subtract(dmd_solution_u, true_solution_u, diff_u);

    double tot_diff_norm_u = sqrt(InnerProduct(MPI_COMM_WORLD, diff_u, diff_u));
    double tot_true_solution_u_norm = sqrt(InnerProduct(MPI_COMM_WORLD,
                                           true_solution_u, true_solution_u));

    if (myid == 0)
    {

        std::cout << "Relative error of DMD temperature (u) at t_final: " << t_final <<
                  " is " << tot_diff_norm_u / tot_true_solution_u_norm << std::endl;
        printf("Elapsed time for solving FOM: %e second\n", fom_timer.RealTime());
        printf("Elapsed time for training DMD: %e second\n",
               dmd_training_timer.RealTime());
        printf("Elapsed time for predicting DMD: %e second\n",
               dmd_prediction_timer.RealTime());
    }

    // 16. Free the used memory.
    delete ode_solver;
    delete pmesh;
    delete result_u;

    MPI_Finalize();

    return 0;
}

ConductionOperator::ConductionOperator(ParFiniteElementSpace &f, double al,
                                       double kap, const Vector &u)
    : TimeDependentOperator(f.GetTrueVSize(), 0.0), fespace(f), M(NULL), K(NULL),
      T(NULL), current_dt(0.0),
      M_solver(f.GetComm()), T_solver(f.GetComm()), z(height)
{
    const double rel_tol = 1e-8;

    M = new ParBilinearForm(&fespace);
    M->AddDomainIntegrator(new MassIntegrator());
    M->Assemble(0); // keep sparsity pattern of M and K the same
    M->FormSystemMatrix(ess_tdof_list, Mmat);

    M_solver.iterative_mode = false;
    M_solver.SetRelTol(rel_tol);
    M_solver.SetAbsTol(0.0);
    M_solver.SetMaxIter(100);
    M_solver.SetPrintLevel(0);
    M_prec.SetType(HypreSmoother::Jacobi);
    M_solver.SetPreconditioner(M_prec);
    M_solver.SetOperator(Mmat);

    alpha = al;
    kappa = kap;

    T_solver.iterative_mode = false;
    T_solver.SetRelTol(rel_tol);
    T_solver.SetAbsTol(0.0);
    T_solver.SetMaxIter(100);
    T_solver.SetPrintLevel(0);
    T_solver.SetPreconditioner(T_prec);

    SetParameters(u);
}

void ConductionOperator::Mult(const Vector &u, Vector &du_dt) const
{
    // Compute:
    //    du_dt = M^{-1}*-K(u)
    // for du_dt
    Kmat.Mult(u, z);
    z.Neg(); // z = -z
    M_solver.Mult(z, du_dt);
}

void ConductionOperator::ImplicitSolve(const double dt,
                                       const Vector &u, Vector &du_dt)
{
    // Solve the equation:
    //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
    // for du_dt
    if (!T)
    {
        T = Add(1.0, Mmat, dt, Kmat);
        current_dt = dt;
        T_solver.SetOperator(*T);
    }
    MFEM_VERIFY(dt == current_dt, ""); // SDIRK methods use the same dt
    Kmat.Mult(u, z);
    z.Neg();
    T_solver.Mult(z, du_dt);
}

void ConductionOperator::SetParameters(const Vector &u)
{
    ParGridFunction u_alpha_gf(&fespace);
    u_alpha_gf.SetFromTrueDofs(u);
    for (int i = 0; i < u_alpha_gf.Size(); i++)
    {
        u_alpha_gf(i) = kappa + alpha*u_alpha_gf(i);
    }

    delete K;
    K = new ParBilinearForm(&fespace);

    GridFunctionCoefficient u_coeff(&u_alpha_gf);

    K->AddDomainIntegrator(new DiffusionIntegrator(u_coeff));
    K->Assemble(0); // keep sparsity pattern of M and K the same
    K->FormSystemMatrix(ess_tdof_list, Kmat);
    delete T;
    T = NULL; // re-compute T on the next ImplicitSolve
}

ConductionOperator::~ConductionOperator()
{
    delete T;
    delete M;
    delete K;
}

double InitialTemperature(const Vector &x)
{
    if (x.Norml2() < 0.5)
    {
        return 2.0;
    }
    else
    {
        return 1.0;
    }
}
