/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

//                       libROM MFEM Example: Parametric_Heat_Conduction (adapted from ex16p.cpp)
//
// Compile with: make heat_conduction_dmdc
//
// =================================================================================
//
// In these examples, the radius of the interface between different source functions, the
// alpha coefficient, and two center location variables are modified.
//
// For Parametric DMD (ex. 1) (radius, interpolation):
//   rm -rf parameters.txt
//   mpirun -np 8 heat_conduction_dmdc -r 0.4 -visit -offline -rdim 16
//   mpirun -np 8 heat_conduction_dmdc -r 0.45 -visit -offline -rdim 16
//   mpirun -np 8 heat_conduction_dmdc -r 0.55 -visit -offline -rdim 16
//   mpirun -np 8 heat_conduction_dmdc -r 0.6 -visit -offline -rdim 16
//   mpirun -np 8 heat_conduction_dmdc -r 0.5 -visit -online -predict
//
// For Parametric DMD (ex. 2) (radius & cx & cy, extrapolation):
//   rm -rf parameters.txt
//   mpirun -np 8 heat_conduction_dmdc -r 0.1 -visit -offline -rdim 16
//   mpirun -np 8 heat_conduction_dmdc -r 0.2 -visit -offline -rdim 16
//   mpirun -np 8 heat_conduction_dmdc -r 0.5 -visit -online -predict
//   (performs well, even though extrapolation)
//   mpirun -np 8 heat_conduction_dmdc -r 0.5 -cx 0.5 -cy 0.5 -visit -online -predict
//   (doesn't perform well)
//   mpirun -np 8 heat_conduction_dmdc -r 0.6 -cx 0.6 -cy 0.6 -visit -offline -rdim 16
//   (let's add another training point)
//   mpirun -np 8 heat_conduction_dmdc -r 0.5 -cx 0.5 -cy 0.5 -visit -online -predict
//   (now performs well)
//
// For Parametric DMD (ex. 3) (alpha, interpolation):
//   rm -rf parameters.txt
//   mpirun -np 8 heat_conduction_dmdc -a 0.1 -visit -offline -rdim 16
//   mpirun -np 8 heat_conduction_dmdc -a 0.15 -visit -offline -rdim 16
//   mpirun -np 8 heat_conduction_dmdc -a 0.25 -visit -offline -rdim 16
//   mpirun -np 8 heat_conduction_dmdc -a 0.3 -visit -offline -rdim 16
//   mpirun -np 8 heat_conduction_dmdc -a 0.2 -visit -online -predict
//
// For Parametric DMD (ex. 4) (alpha, interpolation):
//   rm -rf parameters.txt
//   mpirun -np 8 heat_conduction_dmdc -s 3 -a 0.5 -k 0.5 -o 4 -tf 0.7 -vs 1 -visit -offline -rdim 20
//   mpirun -np 8 heat_conduction_dmdc -s 3 -a 0.55 -k 0.5 -o 4 -tf 0.7 -vs 1 -visit -offline -rdim 20
//   mpirun -np 8 heat_conduction_dmdc -s 3 -a 0.65 -k 0.5 -o 4 -tf 0.7 -vs 1 -visit -offline -rdim 20
//   mpirun -np 8 heat_conduction_dmdc -s 3 -a 0.7 -k 0.5 -o 4 -tf 0.7 -vs 1 -visit -offline -rdim 20
//   mpirun -np 8 heat_conduction_dmdc -s 3 -a 0.6 -k 0.5 -o 4 -tf 0.7 -vs 1 -visit -online -predict
//
// Pointwise snapshots for DMD input:
//   mpirun -np 1 heat_conduction_dmdc -m ../../../examples/data/inline-quad.mesh -pwsnap -pwx 101 -pwy 101
//
// =================================================================================
//
// Description:  This example solves a time dependent nonlinear heat equation
//               problem of the form du/dt = C(u), with a non-linear diffusion
//               operator C(u) = s + \nabla \cdot (\kappa + \alpha u) \nabla u.
//
//               The example demonstrates the use of nonlinear operators (the
//               class ConductionOperator defining C(u)), as well as their
//               implicit time integration. Note that implementing the method
//               ConductionOperator::ImplicitSolve is the only requirement for
//               high-order implicit (SDIRK) time integration. Optional saving
//               with ADIOS2 (adios2.readthedocs.io) is also illustrated.

#include "mfem.hpp"
#include "algo/DMD.h"
#include "linalg/Vector.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include "utils/CSVDatabase.h"
#include "utils/HDFDatabase.h"
#include "mfem/PointwiseSnapshot.hpp"

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
 *     du/dt = M^{-1}(b-Ku)
 *
 *  where u is the vector representing the temperature,
 *  M is the mass matrix, b is the load vector,
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
    mutable Vector b; // auxiliary vector

public:
    ConductionOperator(ParFiniteElementSpace &f, double alpha, double kappa,
                       const Vector &u);

    void GetSource(Vector& s) const;

    virtual void Mult(const Vector &u, Vector &du_dt) const;
    /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
        This is the only requirement for high-order SDIRK implicit integration.*/
    virtual void ImplicitSolve(const double dt, const Vector &u, Vector &k);

    /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
    void SetParameters(const Vector &u);

    virtual ~ConductionOperator();
};

double InitialTemperature(const Vector &x);
double SourceFunction(const Vector &x, const double t);

double radius = 0.5;
double cx = 0.0;
double cy = 0.0;
double t_end = 0.1;
double dt = 1.0e-2;
double amp = 1.0;

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
    double alpha = 1.0e-2;
    double kappa = 0.5;
    double closest_rbf_val = 0.9;
    int rdim = -1;
    bool offline = false;
    bool online = false;
    bool predict = false;
    bool visualization = true;
    bool visit = false;
    int vis_steps = 5;
    bool adios2 = false;
    bool save_dofs = false;
    bool csvFormat = true;
    const char *basename = "";

    bool pointwiseSnapshots = false;
    int pwx = 0;
    int pwy = 0;
    int pwz = 0;

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
    args.AddOption(&radius, "-r", "--radius",
                   "Radius of the radial source location.");
    args.AddOption(&cx, "-cx", "--center_x",
                   "Center offset in the x direction.");
    args.AddOption(&cy, "-cy", "--center_y",
                   "Center offset in the y direction.");
    args.AddOption(&t_end, "-te", "--t_end",
                   "End time of the source.");
    args.AddOption(&amp, "-amp", "--amplitude",
                   "Amplitude of the source.");
    args.AddOption(&closest_rbf_val, "-crv", "--crv",
                   "DMD Closest RBF Value.");
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
    args.AddOption(&offline, "-offline", "--offline", "-no-offline", "--no-offline",
                   "Enable or disable the offline phase.");
    args.AddOption(&online, "-online", "--online", "-no-online", "--no-online",
                   "Enable or disable the online phase.");
    args.AddOption(&predict, "-predict", "--predict", "-no-predict", "--no-predict",
                   "Enable or disable DMD prediction.");
    args.AddOption(&adios2, "-adios2", "--adios2-streams", "-no-adios2",
                   "--no-adios2-streams",
                   "Save data using adios2 streams.");
    args.AddOption(&save_dofs, "-save", "--save", "-no-save", "--no-save",
                   "Enable or disable MFEM DOF solution snapshot files).");
    args.AddOption(&csvFormat, "-csv", "--csv", "-hdf", "--hdf",
                   "Use CSV or HDF format for files output by -save option.");
    args.AddOption(&basename, "-out", "--outputfile-name",
                   "Name of the sub-folder to dump files within the run directory.");
    args.AddOption(&pointwiseSnapshots, "-pwsnap", "--pw-snap", "-no-pwsnap",
                   "--no-pw-snap", "Enable or disable writing pointwise snapshots.");
    args.AddOption(&pwx, "-pwx", "--pwx", "Number of snapshot points in x");
    args.AddOption(&pwy, "-pwy", "--pwy", "Number of snapshot points in y");
    args.AddOption(&pwz, "-pwz", "--pwz", "Number of snapshot points in z");

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

    string outputPath = ".";
    if (save_dofs)
    {
        outputPath = "run";
        if (string(basename) != "") {
            outputPath += "/" + string(basename);
        }
        if (myid == 0)
        {
            const char path_delim = '/';
            string::size_type pos = 0;
            do {
                pos = outputPath.find(path_delim, pos+1);
                string subdir = outputPath.substr(0, pos);
                mkdir(subdir.c_str(), 0777);
            }
            while (pos != string::npos);
        }
    }

    const int check = (int) pointwiseSnapshots + (int) offline + (int) online
                      + (int) save_dofs;
    MFEM_VERIFY(check == 1,
                "Only one of offline, online, save, or pwsnap must be true!");

    if (offline)
    {
        MFEM_VERIFY(rdim != -1, "rdim must be set.");
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

#ifndef MFEM_USE_GSLIB
    if (pointwiseSnapshots) {
        cout << "To use pointwise snapshots, compile with -mg option" << endl;
        MFEM_ABORT("Pointwise snapshots aren't available, since the "
                   "compilation is done without the -mg option");
    }
#else
    CAROM::PointwiseSnapshot *pws = nullptr;
    Vector pwsnap;
    CAROM::Vector *pwsnap_CAROM = nullptr;

    if (pointwiseSnapshots)
    {
        pmesh->EnsureNodes();
        const int dmdDim[3] = {pwx, pwy, pwz};
        pws = new CAROM::PointwiseSnapshot(dim, dmdDim);
        pws->SetMesh(pmesh);

        int snapshotSize = dmdDim[0];
        for (int i=1; i<dim; ++i)
            snapshotSize *= dmdDim[i];

        pwsnap.SetSize(snapshotSize);
        if (myid == 0)
            pwsnap_CAROM = new CAROM::Vector(pwsnap.GetData(), pwsnap.Size(),
                                             true, false);
    }
#endif

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
        mesh_name << outputPath << "/heat_conduction_dmdc_" << to_string(
                      radius) << "_"
                  << to_string(alpha) << "_" << to_string(cx) << "_" << to_string(cy)
                  << "-mesh." << setfill('0') << setw(6) << myid;
        sol_name << outputPath << "/heat_conduction_dmdc_" << to_string(
                     radius) << "_"
                 << to_string(alpha) << "_" << to_string(cx) << "_" << to_string(cy)
                 << "-init." << setfill('0') << setw(6) << myid;
        ofstream omesh(mesh_name.str().c_str());
        omesh.precision(precision);
        pmesh->Print(omesh);
        ofstream osol(sol_name.str().c_str());
        osol.precision(precision);
        u_gf.Save(osol);
    }

    VisItDataCollection visit_dc(outputPath + "/Parametric_Heat_Conduction_" +
                                 to_string(radius) + "_" + to_string(alpha) + "_" + to_string(cx) + "_" +
                                 to_string(cy), pmesh);
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
        const std::string collection_name = outputPath +
                                            "/heat_conduction_dmdc-p-" +
                                            postfix + ".bp";

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

#ifdef MFEM_USE_GSLIB
    if (pointwiseSnapshots)
    {
        pws->GetSnapshot(u_gf, pwsnap);

        ostringstream dmd_filename;
        dmd_filename << "snap_" << to_string(radius) << "_" << to_string(alpha)
                     << "_" << to_string(cx) << "_" << to_string(cy) << "_0";
        if (myid == 0)
        {
            cout << "Writing DMD snapshot at step 0, time 0.0" << endl;
            pwsnap_CAROM->write(dmd_filename.str());
        }
    }
#endif

    StopWatch fom_timer, dmd_training_timer, dmd_prediction_timer;

    fom_timer.Start();

    // 10. Perform time-integration (looping over the time iterations, ti, with a
    //     time-step dt).
    ode_solver->Init(oper);
    double t = 0.0;
    vector<double> ts;
    CAROM::Vector* init = NULL;

    CAROM::Database *db = NULL;
    if (csvFormat)
        db = new CAROM::CSVDatabase();
    else
        db = new CAROM::HDFDatabase();

    vector<int> snap_list;

    fom_timer.Stop();

    CAROM::DMD* dmd_u = NULL;

    if (offline)
    {
        dmd_training_timer.Start();

        // 11. Create DMD object and take initial sample.
        u_gf.SetFromTrueDofs(u);
        dmd_u = new CAROM::DMD(u.Size(), dt);
        dmd_u->takeSample(u.GetData(), t);

        if (myid == 0)
        {
            std::cout << "Taking snapshot at: " << t << std::endl;
        }

        dmd_training_timer.Stop();
    }
    if (online)
    {
        u_gf.SetFromTrueDofs(u);
        init = new CAROM::Vector(u.GetData(), u.Size(), true);
    }

    if (save_dofs && myid == 0)
    {
        if (csvFormat)
        {
            mkdir((outputPath + "/step0").c_str(), 0777);
            db->putDoubleArray(outputPath + "/step0/sol.csv", u.GetData(), u.Size());
        }
        else
        {
            db->create(outputPath + "/dmd_0.hdf");
            db->putDoubleArray("step0sol", u.GetData(), u.Size());
        }
    }

    ts.push_back(t);
    snap_list.push_back(0);

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

        if (offline)
        {
            dmd_training_timer.Start();

            u_gf.SetFromTrueDofs(u);
            dmd_u->takeSample(u.GetData(), t);

            if (myid == 0)
            {
                std::cout << "Taking snapshot at: " << t << std::endl;
            }

            dmd_training_timer.Stop();
        }

        if (save_dofs && myid == 0)
        {
            if (csvFormat)
            {
                mkdir((outputPath + "/step" + to_string(ti)).c_str(), 0777);
                db->putDoubleArray(outputPath + "/step" + to_string(ti) + "/sol.csv",
                                   u.GetData(), u.Size());
            }
            else
                db->putDoubleArray("step" + to_string(ti) + "sol",
                                   u.GetData(), u.Size());
        }

        ts.push_back(t);
        snap_list.push_back(ti);

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

#ifdef MFEM_USE_GSLIB
        if (pointwiseSnapshots)
        {
            pws->GetSnapshot(u_gf, pwsnap);

            ostringstream dmd_filename;
            dmd_filename << "snap_" << to_string(radius) << "_" << to_string(alpha)
                         << "_" << to_string(cx) << "_" << to_string(cy) << "_" << ti;
            if (myid == 0)
            {
                cout << "Writing DMD snapshot at step " << ti << ", time " << t << endl;
                pwsnap_CAROM->write(dmd_filename.str());
            }
        }
#endif

        oper.SetParameters(u);
    }

    if (save_dofs && myid == 0)
    {
        if (csvFormat)
        {
            db->putDoubleVector(outputPath + "/tval.csv", ts, ts.size());
            db->putInteger(outputPath + "/numsnap", snap_list.size());
            db->putIntegerArray(outputPath + "/snap_list.csv", snap_list.data(),
                                snap_list.size());
        }
        else
        {
            db->putDoubleVector("tval", ts, ts.size());
            db->putInteger("numsnap", snap_list.size());
            db->putInteger("snap_bound_size", 0);
            db->putIntegerArray("snap_list", snap_list.data(),
                                snap_list.size());
        }
    }

#ifdef MFEM_USE_ADIOS2
    if (adios2)
    {
        delete adios2_dc;
    }
#endif

    // 12. Save the final solution in parallel. This output can be viewed later
    //     using GLVis: "glvis -np <np> -m heat_conduction_dmdc-mesh -g heat_conduction_dmdc-final".
    {
        ostringstream sol_name;
        sol_name << outputPath << "/heat_conduction_dmdc_" << to_string(
                     radius) << "_"
                 << to_string(alpha) << "_" << to_string(cx) << "_" << to_string(cy)
                 << "-final." << setfill('0') << setw(6) << myid;
        ofstream osol(sol_name.str().c_str());
        osol.precision(precision);
        u_gf.Save(osol);
    }

    // 13. Calculate the DMD modes.
    if (offline || online)
    {
        if (offline)
        {
            if (myid == 0)
            {
                std::cout << "Creating DMD with rdim: " << rdim << std::endl;
            }

            dmd_training_timer.Start();

            dmd_u->train(rdim);

            dmd_training_timer.Stop();

            dmd_u->save(outputPath + "/" + to_string(radius) + "_" + to_string(
                            alpha) + "_" +
                        to_string(cx) + "_" + to_string(cy));

            if (myid == 0)
            {
                std::ofstream fout;
                fout.open("parameters.txt", std::ios::app);
                fout << radius << " " << alpha << " " << cx << " " << cy << std::endl;
                fout.close();
            }
        }

        if (online)
        {
            if (myid == 0)
            {
                std::cout << "Creating DMD using the rdim of the offline phase" << std::endl;
            }

            std::fstream fin("parameters.txt", std::ios_base::in);
            double curr_param;
            std::vector<std::string> dmd_paths;
            std::vector<CAROM::Vector*> param_vectors;

            while (fin >> curr_param)
            {
                double curr_radius = curr_param;
                fin >> curr_param;
                double curr_alpha = curr_param;
                fin >> curr_param;
                double curr_cx = curr_param;
                fin >> curr_param;
                double curr_cy = curr_param;

                dmd_paths.push_back(outputPath + "/" + to_string(curr_radius) + "_" +
                                    to_string(curr_alpha) + "_" + to_string(curr_cx) + "_" +
                                    to_string(curr_cy));
                CAROM::Vector* param_vector = new CAROM::Vector(4, false);
                param_vector->item(0) = curr_radius;
                param_vector->item(1) = curr_alpha;
                param_vector->item(2) = curr_cx;
                param_vector->item(3) = curr_cy;
                param_vectors.push_back(param_vector);
            }
            fin.close();

            CAROM::Vector* desired_param = new CAROM::Vector(4, false);
            desired_param->item(0) = radius;
            desired_param->item(1) = alpha;
            desired_param->item(2) = cx;
            desired_param->item(3) = cy;

            dmd_training_timer.Start();

            CAROM::getParametricDMD(dmd_u, param_vectors, dmd_paths, desired_param,
                                    "G", "LS", closest_rbf_val);

            dmd_u->projectInitialCondition(init);

            dmd_training_timer.Stop();
            delete desired_param;
        }

        if (predict)
        {
            Vector true_solution_u(u.Size());
            true_solution_u = u.GetData();

            dmd_prediction_timer.Start();

            // 14. Predict the state at t_final using DMD.
            if (myid == 0)
            {
                std::cout << "Predicting temperature using DMD at: " << ts[0] << std::endl;
            }

            CAROM::Vector* result_u = dmd_u->predict(ts[0]);
            Vector initial_dmd_solution_u(result_u->getData(), result_u->dim());
            u_gf.SetFromTrueDofs(initial_dmd_solution_u);

            VisItDataCollection dmd_visit_dc(outputPath + "/DMD_Parametric_Heat_Conduction_"
                                             +
                                             to_string(radius) + "_" + to_string(alpha) + "_" +
                                             to_string(cx) + "_" + to_string(cy), pmesh);
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
                        if (myid == 0)
                        {
                            std::cout << "Predicting temperature using DMD at: " << ts[i] << std::endl;
                        }

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
                std::cout << "Relative error of DMD temperature (u) at t_final: "
                          << t_final << " is " << tot_diff_norm_u / tot_true_solution_u_norm << std::endl;
                printf("Elapsed time for predicting DMD: %e second\n",
                       dmd_prediction_timer.RealTime());
            }

            delete result_u;
        }

        if (myid == 0)
        {
            printf("Elapsed time for training DMD: %e second\n",
                   dmd_training_timer.RealTime());
        }
    }

    if (myid == 0)
    {
        printf("Elapsed time for solving FOM: %e second\n", fom_timer.RealTime());
    }

    // 16. Free the used memory.
    delete ode_solver;
    delete pmesh;
    if (offline)
    {
        delete dmd_u;
    }

#ifdef MFEM_USE_GSLIB
    delete pws;
    delete pwsnap_CAROM;
#endif

    MPI_Finalize();

    return 0;
}

ConductionOperator::ConductionOperator(ParFiniteElementSpace &f, double al,
                                       double kap, const Vector &u)
    : TimeDependentOperator(f.GetTrueVSize(), 0.0), fespace(f), M(NULL), K(NULL),
      T(NULL), current_dt(0.0),
      M_solver(f.GetComm()), T_solver(f.GetComm()), z(height), b(height)
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

void ConductionOperator::GetSource(Vector& s_v) const
{
    // Set grid function for f
    ParGridFunction s_gf(&fespace);

    FunctionCoefficient s(SourceFunction);
    s.SetTime(GetTime());
    s_gf.ProjectCoefficient(s);
    s_gf.GetTrueDofs(s_v);
}

void ConductionOperator::Mult(const Vector &u, Vector &du_dt) const
{
    // Compute:
    //    du_dt = M^{-1}*[f-K(u)]
    // for du_dt

    GetSource(b);
    Kmat.Mult(u, z);
    b.Add(-1.0, z);
    M_solver.Mult(b, du_dt);
}

void ConductionOperator::ImplicitSolve(const double dt,
                                       const Vector &u, Vector &du_dt)
{
    // Solve the equation:
    //    du_dt = M^{-1}*[f-K(u + dt*du_dt)]
    // for du_dt
    if (!T)
    {
        T = Add(1.0, Mmat, dt, Kmat);
        current_dt = dt;
        T_solver.SetOperator(*T);
    }
    MFEM_VERIFY(dt == current_dt, ""); // SDIRK methods use the same dt

    GetSource(b);
    Kmat.Mult(u, z);
    b.Add(-1.0, z);
    T_solver.Mult(b, du_dt);
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
    return 1.0;
}

double SourceFunction(const Vector &x, const double t)
{
    Vector y(x);
    Vector c(2);
    c.Elem(0) = cx;
    c.Elem(1) = cy;
    y -= c;
    if (y.Norml2() < radius)
    {
        return amp / (1.0 + exp(0.5*(t-t_end)/dt));
    }
}
