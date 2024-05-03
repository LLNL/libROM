/******************************************************************************
*
* Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
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
//  ./wave_equation_dmdc_1win_rewrite -o 4 -tf 0.10 -lscale 0.15 -center 0.25 -amp 10. -sigma 1.0 --speed 0.5 -online -predict -rdim 4
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
#include "algo/DMDc.h"
#include "linalg/Vector.h"
#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>
#include "utils/CSVDatabase.h"
#include "utils/HDFDatabase.h"

#ifndef _WIN32
#include <sys/stat.h>  // mkdir
#else
#include <direct.h>    // _mkdir
#define mkdir(dir, mode) _mkdir(dir)
#endif

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
    LinearForm *B;

    SparseMatrix Mmat, Kmat, Kmat0;
    SparseMatrix *T; // T = M + dt K
    double current_dt;

    CGSolver M_solver; // Krylov solver for inverting the mass matrix M
    DSmoother M_prec;  // Preconditioner for the mass matrix M

    CGSolver T_solver; // Implicit solver for T = M + fac0*K
    DSmoother T_prec;  // Preconditioner for the implicit solver

    FunctionCoefficient &src_func; // Source function coefficient
    // double sigma; // Nonlinear parameters

    Coefficient *c2;
    mutable Vector z; // auxiliary vector
    // HypreParVector *b; // source vector

public:
    WaveOperator(FiniteElementSpace &f, FunctionCoefficient &s,
                 Array<int> &ess_bdr,double speed);

    // WaveOperator(FiniteElementSpace &f,
    //                 Array<int> &ess_bdr,double speed);

    using SecondOrderTimeDependentOperator::Mult;
    virtual void Mult(const Vector &u, const Vector &du_dt,
                      Vector &d2udt2) const;

    /** Solve the Backward-Euler equation:
        d2udt2 = f(u + fac0*d2udt2,dudt + fac1*d2udt2, t),
        for the unknown d2udt2. */
    using SecondOrderTimeDependentOperator::ImplicitSolve;
    virtual void ImplicitSolve(const double fac0, const double fac1,
                               const Vector &u, const Vector &dudt, Vector &d2udt2);

    void SetForcingTime(const double t);

    void SetParameters(const Vector &u);

    virtual ~WaveOperator();
};

double InitialSolution(const Vector &x);
double ForcingAmp(double t);
double ForcingFunction(const Vector &x, double t);

Vector bb_min, bb_max; // Mesh bounding box
double sigma=1;
double amp = 2;
double lscale = 0.15;
double center = 0.25;
double dt = 1.0e-2;

WaveOperator::WaveOperator(FiniteElementSpace &f,FunctionCoefficient &s,
                           Array<int> &ess_bdr, double speed)
    : SecondOrderTimeDependentOperator(f.GetTrueVSize(), 0.0), fespace(f),
      src_func(s),M(NULL),
      K(NULL),
      T(NULL), current_dt(0.0), z(height)
// WaveOperator::WaveOperator(FiniteElementSpace &f,
//                            Array<int> &ess_bdr, double speed)
//     : SecondOrderTimeDependentOperator(f.GetTrueVSize(), 0.0), fespace(f),M(NULL),
//       K(NULL),
//       T(NULL), current_dt(0.0), z(height)
{
    const double rel_tol = 1e-8;

    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    c2 = new ConstantCoefficient(speed*speed);

    // sigma = sig;

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
    // SetParameters(u);

    B = new LinearForm(&fespace);
    B->AddDomainIntegrator(new DomainLFIntegrator(src_func));
}

void WaveOperator::Mult(const Vector &u, const Vector &du_dt,
                        Vector &d2udt2)  const

{
    // Compute:
    //    d2udt2 = M^{-1}*-K(u)
    // for d2udt2
    Kmat.Mult(u, z);
    z.Neg(); // z = -z
    z += *B;
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
        current_dt = dt;
        T_solver.SetOperator(*T);
    }
    MFEM_VERIFY(dt == current_dt, ""); // SDIRK methods use the same dt
    Kmat0.Mult(u, z);
    z.Neg();
    z += *B;

    for (int i = 0; i < ess_tdof_list.Size(); i++)
    {
        z[ess_tdof_list[i]] = 0.0;
    }
    T_solver.Mult(z, d2udt2);
}

void WaveOperator::SetForcingTime(double t)
{
    src_func.SetTime(t);
    B->Assemble();
    // b = B->Assemble();
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
    delete B;
}

// double InitialSolution(const Vector &x)
// {
//     return exp(-x.Norml2()*x.Norml2()*30);
// }

double InitialSolution(const Vector &x)
{
    Vector X(x.Size());
    for (int i = 0; i < x.Size(); i++)
    {
        double domcenter = (bb_min[i] + bb_max[i]) * 0.5;
        X(i) = 2 * (x(i) - domcenter) / (bb_max[i] - bb_min[i]);
    }
    return 0*exp(-X.Norml2()*X.Norml2()*30);
}

double InitialRate(const Vector &x)
{
    return 0.0;
}

// double ForcingFunction(const Vector &x, double t)
// {
//     Vector X(x.Size());
//     for (int i = 0; i < x.Size(); i++)
//     {
//         X(i) = (x(i) - center);
//     }
//     double r1 = X.Norml2();
//     double rw = amp*2/sqrt(3*sigma)/pow(M_PI,1/4)*(1 - pow(t/sigma,2))*exp(-pow(t/sigma,2)/2);
//     return rw*exp(-r1*r1/lscale/lscale);
// }

double ForcingAmp(double t)
{
    double rw = amp*2/sqrt(3*sigma)/pow(M_PI,1/4)*(1 - pow(t/sigma,
                2))*exp(-pow(t/sigma,2)/2);

    return rw;
}

double ForcingFunction(const Vector &x, double t)
{
    Vector X(x.Size());
    for (int i = 0; i < x.Size(); i++)
    {
        double domcenter = (bb_min[i] + bb_max[i]) * 0.5;
        X(i) = 2 * (x(i) - domcenter) / (bb_max[i] - bb_min[i]);
        X(i) = X(i) - center;
    }
    double r1 = X.Norml2();

    return ForcingAmp(t)*exp(-r1*r1/lscale/lscale);
}

const CAROM::Matrix*
createControlMatrix(std::vector<CAROM::Vector*> snapshots)
{
    CAROM_VERIFY(snapshots.size() > 0);
    CAROM_VERIFY(snapshots[0]->dim() > 0);
    for (int i = 0 ; i < snapshots.size() - 1; i++)
    {
        CAROM_VERIFY(snapshots[i]->dim() == snapshots[i + 1]->dim());
        CAROM_VERIFY(snapshots[i]->distributed() == snapshots[i + 1]->distributed());
    }

    CAROM::Matrix* snapshot_mat = new CAROM::Matrix(snapshots[0]->dim(),
            snapshots.size(),
            snapshots[0]->distributed());

    for (int i = 0; i < snapshots[0]->dim(); i++)
    {
        for (int j = 0; j < snapshots.size(); j++)
        {
            snapshot_mat->item(i, j) = snapshots[j]->item(i);
        }
    }

    return snapshot_mat;
}

int main(int argc, char *argv[])
{
    // 1. Parse command-line options.
    // const char *mesh_file = "../data/star.mesh";
    const char *mesh_file = "";
    const char *ref_dir  = "";
    int ref_levels = 2;
    int order = 2;
    int ode_solver_type = 10;
    double t_final = 0.5;
    // double dt = 1.0e-2;
    double speed = 1.0;
    // double ef = 0.9999;
    int rdim = 15;
    int windowNumSamples = numeric_limits<int>::max();
    bool visualization = true;
    bool visit = false;
    bool dirichlet = true;
    int vis_steps = 5;
    bool offline = false;
    bool online = false;

    double closest_rbf_val = 0.9;
    bool predict = false;
    bool csvFormat = true;
    const char *temp_io_dir = "./outputs_wave_equation_parametric_dmdc";

    std::string io_dir;
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
    args.AddOption(&rdim, "-rdim", "--rdim",
                   "Reduced dimension for DMD.");
    args.AddOption(&windowNumSamples, "-nwinsamp",    "--numwindowsamples",
                   "Number of samples in DMD windows.");
    args.AddOption(&sigma, "-sigma", "--sigma",
                   "Sigma for forcing.");
    args.AddOption(&amp, "-amp", "--amplitude",
                   "Amplitude of forcing.");
    args.AddOption(&lscale, "-lscale", "--lengthscale",
                   "Radius of forcing.");
    args.AddOption(&center, "-center", "--center",
                   "Center of forcing.");
    args.AddOption(&offline, "-offline", "--offline", "-no-offline", "--no-offline",
                   "Enable or disable the offline phase.");
    args.AddOption(&online, "-online", "--online", "-no-online", "--no-online",
                   "Enable or disable the online phase.");
    args.AddOption(&closest_rbf_val, "-crv", "--crv",
                   "DMD Closest RBF Value.");
    args.AddOption(&predict, "-predict", "--predict", "-no-predict", "--no-predict",
                   "Enable or disable DMD prediction.");
    args.AddOption(&csvFormat, "-csv", "--csv", "-hdf", "--hdf",
                   "Use CSV or HDF format for files output by -save option.");
    args.AddOption(&temp_io_dir, "-io", "--io-dir-name",
                   "Name of the sub-folder to load/dump input/output files within the current directory.");

    args.Parse();
    if (!args.Good())
    {
        args.PrintUsage(cout);
        return 1;
    }
    args.PrintOptions(cout);

    io_dir = temp_io_dir;
    mkdir(io_dir.c_str(), 0777);

    if (rdim <= 0 && rdim != -1) {
        cout << "rdim is set to " << rdim <<
             " , rdim can only be a positive integer or -1" << endl;
        return 1;
    }

    else if (rdim != -1)
    {
        cout << "rdim is set to " << rdim << endl;
    }

    // 2. Read the mesh from the given mesh file. We can handle triangular,
    //    quadrilateral, tetrahedral and hexahedral meshes with the same code.
    Mesh *mesh;
    if (mesh_file == "")
    {
        mesh = new Mesh(Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL));
    }
    else
    {
        mesh = new Mesh(mesh_file, 1, 1);
    }
    int dim = mesh->Dimension();

    mesh->GetBoundingBox(bb_min, bb_max, max(order, 1));
    // Mesh *mesh = new Mesh(mesh_file, 1, 1);
    // int dim = mesh->Dimension();

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

    FunctionCoefficient s(ForcingFunction);
    WaveOperator oper(fespace,s, ess_bdr, speed);
    // WaveOperator oper(fespace, ess_bdr, speed);


    {
        stringstream mesh_name, sol_name;
        mesh_name << io_dir << "/wave_equation_parametric_dmdc" << "_"  << to_string(
                      sigma) <<
                  "-mesh." << setfill('0') << setw(6);
        sol_name << io_dir << "/wave_equation_parametric_dmdc" << "_" << to_string(
                     sigma) <<  "-init." << setfill('0') << setw(
                     6);
        ofstream omesh(mesh_name.str().c_str());
        omesh.precision(precision);
        mesh->Print(omesh);
        ofstream osol(sol_name.str().c_str());
        osol.precision(precision);
        u_gf.Save(osol);
        dudt_gf.Save(osol);
    }

    VisItDataCollection visit_dc(io_dir + "/wave_equation_parametric_dmdc_FOM_" +
                                 to_string(sigma), mesh);
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
    double f[1];
    fom_timer.Stop();


    dmd_training_timer.Start();
    CAROM::Vector* init = NULL;
    f[0] = ForcingAmp(t);
    int curr_window = 0;
    // CAROM::DMDc* dmd_u;
    vector<CAROM::DMDc*> dmd_u;
    int dim_c = 1;

    cout << "timestart"<< endl;
    // std::vector<std::vector<CAROM::Vector*>> d_controls; // stores controls.
    std::vector<CAROM::Vector*> d_controls; // stores controls.
    cout << "timestart"<< endl;

    CAROM::Vector* control = new CAROM::Vector(f, dim_c, false);
    cout << "timestart"<< endl;

    // d_controls[curr_window].push_back(control);
    d_controls.push_back(control);
    cout << "timestart"<< endl;


    if (offline)
    {
        dmd_u.push_back(new CAROM::DMDc(u.Size(),dim_c, dt));
        dmd_u[curr_window]->takeSample(u.GetData(), t,f,false);
        // dmd_u = new CAROM::DMDc(u.Size(),dim_c, dt);
        // dmd_u->takeSample(u.GetData(), t,f,false);
    }
    vector<Vector> u_winend;

    if (online)
    {   
        // vector<Vector> u_winend;
        init = new CAROM::Vector(u.GetData(), u.Size(), true);
    }

    ts.push_back(t);

    dmd_training_timer.Stop();
    bool last_step = false;
    cout << "timestart"<< endl;
    for (int ti = 1; !last_step; ti++)
    {
        fom_timer.Start();
        if (t + dt >= t_final - dt/2)
        {
            last_step = true;
        }


        oper.SetForcingTime(t + dt); // how do alpha methods work

        ode_solver->Step(u, dudt, t, dt);
        fom_timer.Stop();

        dmd_training_timer.Start();
        f[0] = ForcingAmp(t);
        if (offline)
        {
            if ((ti % windowNumSamples) != 0)
            {
                dmd_u[curr_window]->takeSample(u.GetData(), t,f,false);
                CAROM::Vector* control = new CAROM::Vector(f, dim_c, false);
                // d_controls[curr_window].push_back(control);
                d_controls.push_back(control);
            }
            else
            {
                dmd_u[curr_window]->takeSample(u.GetData(), t,f,true);
            }
            if (last_step || (ti % windowNumSamples) == 0)
            {
                cout << "step " << ti << ", t = " << t << endl;

                if (rdim != -1) {
                    cout << "Creating DMDc with rdim: " << rdim << " at window index: " <<
                         curr_window << endl;
                    dmd_u[curr_window]->train(rdim);
                }

                dmd_u[curr_window]->save(io_dir + "/" +  to_string(sigma) + "_" + to_string(
                                             curr_window));

                if (curr_window==0) //only write params once per run
                {
                    std::ofstream fout;
                    fout.open(io_dir + "/parameters.txt", std::ios::app);
                    fout << sigma << std::endl;
                    fout.close();
                }

                // const CAROM::Matrix* control_mat = createControlMatrix(d_controls[curr_window]);
                const CAROM::Matrix* control_mat = createControlMatrix(d_controls);

                std::string full_file_name;

                full_file_name = io_dir + "/" + to_string(sigma) + "_" + to_string(
                                     curr_window) + "_control";
                control_mat->write(full_file_name);
                d_controls.clear();

                if (!last_step) {
                    curr_window++;
                    dmd_u.push_back(new CAROM::DMDc(u.Size(),dim_c, dt));
                    dmd_u[curr_window]->takeSample(u.GetData(), t, f, false);
                    CAROM::Vector* control = new CAROM::Vector(f, dim_c,
                            false); //controls at current time
                    // d_controls[curr_window].push_back(control);
                    d_controls.push_back(control);
                }

            }

        }

        if (last_step || (ti % windowNumSamples) == 0)
        {
            u_winend.push_back(u);
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
        ostringstream sol_name;
        sol_name << io_dir << "/wave_equation_parametric_dmdc" << "_" <<  to_string(
                     sigma) <<  "-final." << setfill('0') << setw(
                     6);
        ofstream osol(sol_name.str().c_str());
        osol.precision(precision);
        u_gf.Save(osol);
        dudt_gf.Save(osol);
    }

    // 13. Calculate the DMDc modes.
    if (online || predict)
    {

        std::cout << "Creating DMDc using the rdim of the offline phase" << std::endl;

        if (online)
        {

            std::cout << "Creating DMDc using the rdim of the offline phase" << std::endl;
            // std::cout << ((int)round((t_final/dt)/windowNumSamples) ) << std::endl;

            for ( int curr_win = 0;
                    curr_win < ((int)round((t_final/dt)/windowNumSamples)); curr_win++)
            {

                std::fstream fin(io_dir + "/parameters.txt", std::ios_base::in);
                double curr_param;
                std::vector<std::string> dmdc_paths;
                std::vector<CAROM::Matrix*> controls;
                std::vector<CAROM::Vector*> param_vectors;


                while (fin >> curr_param)
                {
                    double curr_sigma = curr_param;

                    dmdc_paths.push_back(io_dir + "/" + to_string(curr_sigma) +  "_" + to_string(
                                             curr_win));

                    CAROM::Matrix* curr_control = new CAROM::Matrix();
                    curr_control->read(io_dir + "/" + to_string(curr_sigma) + "_" + to_string(
                                           curr_win) + "_control");
                    controls.push_back(curr_control);

                    CAROM::Vector* param_vector = new CAROM::Vector(1, false);
                    param_vector->item(0) = curr_sigma;
                    param_vectors.push_back(param_vector);
                }
                fin.close();

                std::cout << "hereparamvec" << std::endl;

                CAROM::Vector* desired_param = new CAROM::Vector(1, false);
                desired_param->item(0) = sigma;

                std::cout << "heredesireparam" << std::endl;

                dmd_training_timer.Start();


                CAROM::Matrix* controls_interpolated = new CAROM::Matrix();

                std::cout << "heredmdcstart" << std::endl;

                CAROM::DMDc* dmd_currwin = NULL;


                // CAROM::getParametricDMDc(dmd_u[curr_window], param_vectors, dmdc_paths,
                //                          controls,
                //                          controls_interpolated, desired_param, "G", "LS", closest_rbf_val);

                CAROM::getParametricDMDc(dmd_currwin, param_vectors, dmdc_paths,
                                         controls,
                                         controls_interpolated, desired_param, "G", "LS", closest_rbf_val);


                std::cout << "heredmdcend" << std::endl;

                dmd_u.push_back(dmd_currwin);

                // CAROM::DMDc* dmd_currwin = NULL;

                if (curr_win == 0)
                {
                    dmd_u[curr_win]->project(init,controls_interpolated);
                    // dmd_currwin->project(init,controls_interpolated);
                    // dmd_u[curr_window] = dmd_currwin;
                }
                else
                {
                    // init = NULL;
                    cout << "time win" << dt*windowNumSamples  << endl;
                    CAROM::Vector* IC = dmd_u[curr_win-1]->predict(dt*windowNumSamples*
                                                       curr_win );
                    init = new CAROM::Vector(IC->getData(), IC->dim(),true);
                    dmd_u[curr_win]->project(init,controls_interpolated);
                }


                dmd_training_timer.Stop();

                delete desired_param;
                delete controls_interpolated;
                for (auto m : param_vectors)
                    delete m;
                for (auto m : controls)
                    delete m;
            }
        }

        cout << "times" << ts[1] << endl;

        if (predict)
        {
            Vector true_solution_u(u.Size());
            true_solution_u = u.GetData();

            // dmd_prediction_timer.Start();

            // 10. Predict using DMDc.
            cout << "Predicting temperature using DMDc" << endl;
            CAROM::Vector* result_u = dmd_u[curr_window]->predict(ts[0]);
            Vector initial_dmd_solution_u(result_u->getData(), result_u->dim());
            u_gf.SetFromTrueDofs(initial_dmd_solution_u);

            VisItDataCollection dmd_visit_dc(io_dir + "/wave_equation_parametric_dmdc_ROM_"
                                             + to_string(sigma), mesh);

            dmd_visit_dc.RegisterField("solution", &u_gf);
            curr_window = 0;
            if (visit) {
                dmd_visit_dc.SetCycle(0);
                dmd_visit_dc.SetTime(0.0);
                dmd_visit_dc.Save();
            }

            delete result_u;

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
                        true_solution_u = u_winend[curr_window].GetData();
                        result_u = dmd_u[curr_window]->predict(ts[i]);
                        Vector dmd_solution_u(result_u->getData(), result_u->dim());
            Vector diff_u(true_solution_u.Size());
            subtract(dmd_solution_u, true_solution_u, diff_u);

            double tot_diff_norm_u = sqrt(InnerProduct(MPI_COMM_WORLD, diff_u, diff_u));
            double tot_true_solution_u_norm = sqrt(InnerProduct(MPI_COMM_WORLD,
                                                   true_solution_u, true_solution_u));

            std::cout << "Relative error of DMDc temperature (u) at end of window: " << curr_window <<
                      " is " << tot_diff_norm_u / tot_true_solution_u_norm << std::endl;
                        delete dmd_u[curr_window];
                        curr_window++;
                    }
                }
            }
            true_solution_u = u.GetData();
            dmd_prediction_timer.Start();
            result_u = dmd_u[curr_window]->predict(t_final);
            dmd_prediction_timer.Stop();


            // 15. Calculate the relative error between the DMDc final solution and the true solution.
            Vector dmd_solution_u(result_u->getData(), result_u->dim());
            Vector diff_u(true_solution_u.Size());
            subtract(dmd_solution_u, true_solution_u, diff_u);

            double tot_diff_norm_u = sqrt(InnerProduct(MPI_COMM_WORLD, diff_u, diff_u));
            double tot_true_solution_u_norm = sqrt(InnerProduct(MPI_COMM_WORLD,
                                                   true_solution_u, true_solution_u));

            std::cout << "Relative error of DMDc temperature (u) at t_final: " << t_final <<
                      " is " << tot_diff_norm_u / tot_true_solution_u_norm << std::endl;
            printf("Elapsed time for solving FOM: %e second\n", fom_timer.RealTime());
            printf("Elapsed time for training DMDc: %e second\n",
                   dmd_training_timer.RealTime());
            printf("Elapsed time for predicting DMDc: %e second\n",
                   dmd_prediction_timer.RealTime());

            delete result_u;

        }

    }

    // 12. Free the used memory.
    delete ode_solver;
    delete mesh;
    // delete dmd_u;

    MPI_Finalize();

    return 0;
}
