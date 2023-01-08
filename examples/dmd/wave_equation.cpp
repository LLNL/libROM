/******************************************************************************
*
* Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
* and other libROM project developers. See the top-level COPYRIGHT
* file for details.
*
* SPDX-License-Identifier: (Apache-2.0 OR MIT)
*
*****************************************************************************/
//
//                       libROM MFEM Example: Wave_Equation (adapted from ex23.cpp)
//
// Compile with: make wave_equation
//
// =================================================================================
//
// Sample runs and results for Time-Windowing DMD:
//
// Command 1:
//  wave_equation -o 4 -tf 5 -nwinsamp 25
//
// Output 1:
// Relative error of DMD solution (u) at t_final: 5 is 3.0483662e-05
// Elapsed time for solving FOM: 3.122185e+00 second
// Elapsed time for training DMD: 6.904051e-01 second
// Elapsed time for predicting DMD: 2.496171e-03 second
//
// =================================================================================
//
//
// Description:  This example solves the wave equation problem of the form:
//
//                               d^2u/dt^2 = c^2 \Delta u.
//
//               The example demonstrates the use of time dependent operators,
//               implicit solvers and second order time integration.
//
//               We recommend viewing examples on DG Advection and Non-linear Elasticity
//               before viewing this example.
//
#include "mfem.hpp"
#include "algo/DMD.h"
#include "algo/NonuniformDMD.h"
#include "linalg/Vector.h"
#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

/** After spatial discretization, the conduction model can be written as:
 *
 *     d^2u/dt^2 = M^{-1}(-Ku)
 *
 *  where u is the vector representing the temperature, M is the mass matrix,
 *  and K is the diffusion operator with diffusivity depending on u:
 *  (\kappa + \alpha u).
 *
 *  Class WaveOperator represents the right-hand side of the above ODE.
 */
class WaveOperator : public SecondOrderTimeDependentOperator
{
protected:
    FiniteElementSpace &fespace;
    Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.

    BilinearForm *M;
    BilinearForm *K;

    SparseMatrix Mmat, Kmat, Kmat0;
    SparseMatrix *T; // T = M + dt K
    double current_dt;

    CGSolver M_solver; // Krylov solver for inverting the mass matrix M
    DSmoother M_prec;  // Preconditioner for the mass matrix M

    CGSolver T_solver; // Implicit solver for T = M + fac0*K
    DSmoother T_prec;  // Preconditioner for the implicit solver

    Coefficient *c2;
    mutable Vector z; // auxiliary vector

public:
    WaveOperator(FiniteElementSpace &f, Array<int> &ess_bdr,double speed);

    using SecondOrderTimeDependentOperator::Mult;
    virtual void Mult(const Vector &u, const Vector &du_dt,
                      Vector &d2udt2) const;

    /** Solve the Backward-Euler equation:
        d2udt2 = f(u + fac0*d2udt2,dudt + fac1*d2udt2, t),
        for the unknown d2udt2. */
    using SecondOrderTimeDependentOperator::ImplicitSolve;
    virtual void ImplicitSolve(const double fac0, const double fac1,
                               const Vector &u, const Vector &dudt, Vector &d2udt2);

    void SetParameters(const Vector &u);

    virtual ~WaveOperator();
};

WaveOperator::WaveOperator(FiniteElementSpace &f,
                           Array<int> &ess_bdr, double speed)
    : SecondOrderTimeDependentOperator(f.GetTrueVSize(), 0.0), fespace(f), M(NULL),
      K(NULL),
      T(NULL), current_dt(0.0), z(height)
{
    const double rel_tol = 1e-8;

    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    c2 = new ConstantCoefficient(speed*speed);

    K = new BilinearForm(&fespace);
    K->AddDomainIntegrator(new DiffusionIntegrator(*c2));
    K->Assemble();

    Array<int> dummy;
    K->FormSystemMatrix(dummy, Kmat0);
    K->FormSystemMatrix(ess_tdof_list, Kmat);

    M = new BilinearForm(&fespace);
    M->AddDomainIntegrator(new MassIntegrator());
    M->Assemble();
    M->FormSystemMatrix(ess_tdof_list, Mmat);

    M_solver.iterative_mode = false;
    M_solver.SetRelTol(rel_tol);
    M_solver.SetAbsTol(0.0);
    M_solver.SetMaxIter(30);
    M_solver.SetPrintLevel(0);
    M_solver.SetPreconditioner(M_prec);
    M_solver.SetOperator(Mmat);

    T_solver.iterative_mode = false;
    T_solver.SetRelTol(rel_tol);
    T_solver.SetAbsTol(0.0);
    T_solver.SetMaxIter(100);
    T_solver.SetPrintLevel(0);
    T_solver.SetPreconditioner(T_prec);

    T = NULL;
}

void WaveOperator::Mult(const Vector &u, const Vector &du_dt,
                        Vector &d2udt2)  const
{
    // Compute:
    //    d2udt2 = M^{-1}*-K(u)
    // for d2udt2
    Kmat.Mult(u, z);
    z.Neg(); // z = -z
    M_solver.Mult(z, d2udt2);
}

void WaveOperator::ImplicitSolve(const double fac0, const double fac1,
                                 const Vector &u, const Vector &dudt, Vector &d2udt2)
{
    // Solve the equation:
    //    d2udt2 = M^{-1}*[-K(u + fac0*d2udt2)]
    // for d2udt2
    if (!T)
    {
        T = Add(1.0, Mmat, fac0, Kmat);
        T_solver.SetOperator(*T);
    }
    Kmat0.Mult(u, z);
    z.Neg();

    for (int i = 0; i < ess_tdof_list.Size(); i++)
    {
        z[ess_tdof_list[i]] = 0.0;
    }
    T_solver.Mult(z, d2udt2);
}

void WaveOperator::SetParameters(const Vector &u)
{
    delete T;
    T = NULL; // re-compute T on the next ImplicitSolve
}

WaveOperator::~WaveOperator()
{
    delete T;
    delete M;
    delete K;
    delete c2;
}

double InitialSolution(const Vector &x)
{
    return exp(-x.Norml2()*x.Norml2()*30);
}

double InitialRate(const Vector &x)
{
    return 0.0;
}

int main(int argc, char *argv[])
{
    // 1. Parse command-line options.
    const char *mesh_file = "../data/star.mesh";
    const char *ref_dir  = "";
    int ref_levels = 2;
    int order = 2;
    int ode_solver_type = 10;
    double t_final = 0.5;
    double dt = 1.0e-2;
    double speed = 1.0;
    double ef = 0.9999;
    int rdim = -1;
    int windowNumSamples = numeric_limits<int>::max();
    bool visualization = true;
    bool visit = false;
    bool dirichlet = true;
    int vis_steps = 5;
    int precision = 8;
    cout.precision(precision);
    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&ref_levels, "-r", "--refine",
                   "Number of times to refine the mesh uniformly.");
    args.AddOption(&order, "-o", "--order",
                   "Order (degree) of the finite elements.");
    args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                   "ODE solver: [0--10] - GeneralizedAlpha(0.1 * s),\n\t"
                   "\t   11 - Average Acceleration, 12 - Linear Acceleration\n"
                   "\t   13 - CentralDifference, 14 - FoxGoodwin");
    args.AddOption(&t_final, "-tf", "--t-final",
                   "Final time; start time is 0.");
    args.AddOption(&dt, "-dt", "--time-step",
                   "Time step.");
    args.AddOption(&speed, "-c", "--speed",
                   "Wave speed.");
    args.AddOption(&dirichlet, "-dir", "--dirichlet", "-neu",
                   "--neumann",
                   "BC switch.");
    args.AddOption(&ref_dir, "-r", "--ref",
                   "Reference directory for checking final solution.");
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
    args.AddOption(&windowNumSamples, "-nwinsamp",    "--numwindowsamples",
                   "Number of samples in DMD windows.");
    args.Parse();
    if (!args.Good())
    {
        args.PrintUsage(cout);
        return 1;
    }
    args.PrintOptions(cout);

    if (rdim <= 0 && rdim != -1) {
        cout << "rdim is set to " << rdim <<
             " , rdim can only be a positive integer or -1" << endl;
        return 1;
    }

    if (ef <= 0.0)
    {
        cout << "ef must be a positive, it is " << ef << endl;
        return 1;
    }
    else if (rdim != -1)
    {
        cout << "rdim is set to " << rdim << endl;
    }

    // 2. Read the mesh from the given mesh file. We can handle triangular,
    //    quadrilateral, tetrahedral and hexahedral meshes with the same code.
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    // 3. Define the ODE solver used for time integration. Several second order
    //    time integrators are available.
    SecondOrderODESolver *ode_solver;
    switch (ode_solver_type)
    {
    // Implicit methods
    case 0:
        ode_solver = new GeneralizedAlpha2Solver(0.0);
        break;
    case 1:
        ode_solver = new GeneralizedAlpha2Solver(0.1);
        break;
    case 2:
        ode_solver = new GeneralizedAlpha2Solver(0.2);
        break;
    case 3:
        ode_solver = new GeneralizedAlpha2Solver(0.3);
        break;
    case 4:
        ode_solver = new GeneralizedAlpha2Solver(0.4);
        break;
    case 5:
        ode_solver = new GeneralizedAlpha2Solver(0.5);
        break;
    case 6:
        ode_solver = new GeneralizedAlpha2Solver(0.6);
        break;
    case 7:
        ode_solver = new GeneralizedAlpha2Solver(0.7);
        break;
    case 8:
        ode_solver = new GeneralizedAlpha2Solver(0.8);
        break;
    case 9:
        ode_solver = new GeneralizedAlpha2Solver(0.9);
        break;
    case 10:
        ode_solver = new GeneralizedAlpha2Solver(1.0);
        break;
    case 11:
        ode_solver = new AverageAccelerationSolver();
        break;
    case 12:
        ode_solver = new LinearAccelerationSolver();
        break;
    case 13:
        ode_solver = new CentralDifferenceSolver();
        break;
    case 14:
        ode_solver = new FoxGoodwinSolver();
        break;
    default:
        cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
        delete mesh;
        return 3;
    }

    // 4. Refine the mesh to increase the resolution. In this example we do
    //    'ref_levels' of uniform refinement, where 'ref_levels' is a
    //    command-line parameter.
    for (int lev = 0; lev < ref_levels; lev++)
    {
        mesh->UniformRefinement();
    }

    // 5. Define the vector finite element space representing the current and the
    //    initial temperature, u_ref.
    H1_FECollection fe_coll(order, dim);
    FiniteElementSpace fespace(mesh, &fe_coll);

    int true_size = fespace.GetTrueVSize();
    cout << "Number of temperature unknowns: " << true_size << endl;
    Array<int> true_offset(3);
    true_offset[0] = 0;
    true_offset[1] = true_size;
    true_offset[2] = 2*true_size;
    BlockVector vx(true_offset);

    GridFunction u_gf(&fespace);
    GridFunction dudt_gf(&fespace);
    dudt_gf.MakeTRef(&fespace, vx, true_offset[0]);
    u_gf.MakeTRef(&fespace, vx, true_offset[1]);

    // 6. Set the initial conditions for u. All boundaries are considered
    //    natural.
    FunctionCoefficient u_0(InitialSolution);
    u_gf.ProjectCoefficient(u_0);
    Vector u;
    u_gf.GetTrueDofs(u);
    FunctionCoefficient dudt_0(InitialRate);
    dudt_gf.ProjectCoefficient(dudt_0);
    Vector dudt;
    dudt_gf.GetTrueDofs(dudt);
    // 7. Initialize the conduction operator and the visualization.
    Array<int> ess_bdr;
    if (mesh->bdr_attributes.Size())
    {
        ess_bdr.SetSize(mesh->bdr_attributes.Max());

        if (dirichlet)
        {
            ess_bdr = 1;
        }
        else
        {
            ess_bdr = 0;
        }
    }

    WaveOperator oper(fespace, ess_bdr, speed);

    {
        ofstream omesh("wave_equation.mesh");
        omesh.precision(precision);
        mesh->Print(omesh);
        ofstream osol("wave_equation-init.gf");
        osol.precision(precision);
        u_gf.Save(osol);
        dudt_gf.Save(osol);
    }

    VisItDataCollection visit_dc("Wave_Equation", mesh);
    visit_dc.RegisterField("solution", &u_gf);
    visit_dc.RegisterField("rate", &dudt_gf);
    if (visit)
    {
        visit_dc.SetCycle(0);
        visit_dc.SetTime(0.0);
        visit_dc.Save();
    }

    socketstream sout;
    if (visualization)
    {
        char vishost[] = "localhost";
        int  visport   = 19916;
        sout.open(vishost, visport);
        if (!sout)
        {
            cout << "Unable to connect to GLVis server at "
                 << vishost << ':' << visport << endl;
            visualization = false;
            cout << "GLVis visualization disabled.\n";
        }
        else
        {
            sout.precision(precision);
            sout << "solution\n" << *mesh << dudt_gf;
            sout << flush;
        }
    }

    StopWatch fom_timer, dmd_training_timer, dmd_prediction_timer;
    fom_timer.Start();

    // 8. Perform time-integration (looping over the time iterations, ti, with a
    //    time-step dt).
    ode_solver->Init(oper);
    double t = 0.0;
    vector<double> ts;
    fom_timer.Stop();
    dmd_training_timer.Start();
    int curr_window = 0;
    vector<CAROM::DMD*> dmd_u;
    dmd_u.push_back(new CAROM::DMD(u.Size(), dt));
    dmd_u[curr_window]->takeSample(u, t);
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

        ode_solver->Step(u, dudt, t, dt);
        fom_timer.Stop();

        dmd_training_timer.Start();
        dmd_u[curr_window]->takeSample(u, t);

        if (last_step || (ti % windowNumSamples) == 0)
        {
            cout << "step " << ti << ", t = " << t << endl;

            if (rdim != -1) {
                cout << "Creating DMD with rdim: " << rdim << " at window index: " <<
                     curr_window << endl;
                dmd_u[curr_window]->train(rdim);
            }

            else {
                cout << "Creating DMD with energy fraction: " << ef << " at window index: " <<
                     curr_window << endl;
                dmd_u[curr_window]->train(ef);
            }

            if (!last_step) {
                curr_window++;
                dmd_u.push_back(new CAROM::DMD(u.Size(), dt));
                dmd_u[curr_window]->takeSample(u, t);
            }

        }
        ts.push_back(t);
        dmd_training_timer.Stop();

        if (last_step || (ti % vis_steps) == 0)
        {

            u_gf.SetFromTrueDofs(u);
            dudt_gf.SetFromTrueDofs(dudt);

            if (visualization)
            {
                sout << "solution\n" << *mesh << u_gf << flush;
            }

            if (visit)
            {
                visit_dc.SetCycle(ti);
                visit_dc.SetTime(t);
                visit_dc.Save();
            }
        }
        oper.SetParameters(u);
    }

    // 9. Save the final solution. This output can be viewed later using GLVis:
    //    "glvis -m wave_equation.mesh -g wave_equation-final.gf".
    {
        ofstream osol("wave_equation-final.gf");
        osol.precision(precision);
        u_gf.Save(osol);
        dudt_gf.Save(osol);
    }

    // 10. Predict the state at t_final using DMD.
    dmd_prediction_timer.Start();
    cout << "Predicting temperature using DMD" << endl;
    CAROM::Vector* result_u = nullptr;
    VisItDataCollection dmd_visit_dc("DMD_Wave_Equation", mesh);
    dmd_visit_dc.RegisterField("solution", &u_gf);
    curr_window = 0;
    if (visit) {
        result_u = dmd_u[curr_window]->predict(ts[0]);
        Vector initial_dmd_solution_u(result_u->getData(), result_u->dim());
        u_gf.SetFromTrueDofs(initial_dmd_solution_u);
        dmd_visit_dc.SetCycle(0);
        dmd_visit_dc.SetTime(0.0);
        dmd_visit_dc.Save();
        delete result_u;
    }

    for (int i = 1; i < ts.size(); i++)
    {
        if (i == ts.size() - 1 || (i % vis_steps) == 0)
        {
            if (visit)
            {
                result_u = dmd_u[curr_window]->predict(ts[i]);
                Vector dmd_solution_u(result_u->getData(), result_u->dim());
                u_gf.SetFromTrueDofs(dmd_solution_u);
                dmd_visit_dc.SetCycle(i);
                dmd_visit_dc.SetTime(ts[i]);
                dmd_visit_dc.Save();
                delete result_u;
            }

            if (i % windowNumSamples == 0 && i < ts.size()-1)
            {
                delete dmd_u[curr_window];
                curr_window++;
            }
        }
    }
    dmd_prediction_timer.Stop();
    result_u = dmd_u[curr_window]->predict(t_final);

    // 11. Calculate the relative error between the DMD final solution and the true solution.
    Vector dmd_solution_u(result_u->getData(), result_u->dim());
    Vector diff_u(u.Size());
    subtract(dmd_solution_u, u, diff_u);
    double tot_diff_norm_u = sqrt(InnerProduct(diff_u, diff_u));
    double tot_true_solution_u_norm = sqrt(InnerProduct(
            u, u));

    cout << "Relative error of DMD solution (u) at t_final: " << t_final <<
         " is " << tot_diff_norm_u / tot_true_solution_u_norm << endl;
    printf("Elapsed time for solving FOM: %e second\n", fom_timer.RealTime());
    printf("Elapsed time for training DMD: %e second\n",
           dmd_training_timer.RealTime());
    printf("Elapsed time for predicting DMD: %e second\n",
           dmd_prediction_timer.RealTime());

    // 12. Free the used memory.
    delete ode_solver;
    delete mesh;
    delete result_u;
    delete dmd_u[curr_window];
    return 0;
}
