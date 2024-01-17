/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

//      libROM MFEM Example: DMDc for Heat_Conduction (adapted from ex16p.cpp)
//
// Compile with: make heat_conduction_dmdc
//
// =================================================================================
//
// Sample runs and results for DMDc:
//
// Command 1:
//  mpirun -np 8 parametric_dmdc_heat_conduction -s 1 -offline --kappa 0.3 -rdim 6
//  mpirun -np 8 parametric_dmdc_heat_conduction -s 1 -offline --kappa 0.35 -rdim 6
//  mpirun -np 8 parametric_dmdc_heat_conduction -s 1 -offline --kappa 0.45 -rdim 6
//  mpirun -np 8 parametric_dmdc_heat_conduction -s 1 -offline --kappa 0.5 -rdim 6
//
//  mpirun -np 8 parametric_dmdc_heat_conduction -s 1 --kappa 0.4 -online -predict

// Output 1:
//   Relative error of DMDc temperature (u) at t_final: 0.5 is 0.0021705658
//
// Command 2:
//   mpirun -np 8 heat_conduction_dmdc -s 1 -a 0.5 -k 0.5 -o 4 -tf 0.7 -vs 1 -visit
//
// Output 2:
//   Relative error of DMDc temperature (u) at t_final: 0.7 is 0.00099736216
//
// =================================================================================
//
// Description:  This example solves a time dependent nonlinear heat equation
//               problem of the form du/dt = C(u) + s, with a non-linear diffusion
//               operator C(u) = \nabla \cdot (\kappa + \alpha u) \nabla u
//               and time-varying external inlet and outlet source.
//               The inlet and the outlet is located at (0,0) and (0.5,0.5)
//               in the reference domain [-1,1]^d, where the shut down time and
//               the amplitude of the sources are the control variables.
//
//               The example demonstrates the use of nonlinear operators (the
//               class ConductionOperator defining C(u)), as well as their
//               implicit time integration. Note that implementing the method
//               ConductionOperator::ImplicitSolve is the only requirement for
//               high-order implicit (SDIRK) time integration. Optional saving
//               with ADIOS2 (adios2.readthedocs.io) is also illustrated.

#include "mfem.hpp"
#include "algo/DMDc.h"
//#include "algo/DMD.h"
#include "linalg/Vector.h"
//#include "linalg/Matrix.h"
#include <cmath>
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
    ParLinearForm *B;

    HypreParMatrix Mmat;
    HypreParMatrix Kmat;
    HypreParMatrix *T; // T = M + dt K
    double current_dt;

    CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
    HypreSmoother M_prec; // Preconditioner for the mass matrix M

    CGSolver T_solver;    // Implicit solver for T = M + dt K
    HypreSmoother T_prec; // Preconditioner for the implicit solver

    FunctionCoefficient &src_func; // Source function coefficient
    double alpha, kappa; // Nonlinear parameters

    mutable Vector z; // auxiliary vector
    HypreParVector *b; // source vector

public:
    ConductionOperator(ParFiniteElementSpace &f, FunctionCoefficient &s,
                       double alpha, double kappa, const Vector &u);

    virtual void Mult(const Vector &u, Vector &du_dt) const;
    /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
        This is the only requirement for high-order SDIRK implicit integration.*/
    virtual void ImplicitSolve(const double dt, const Vector &u, Vector &k);

    /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
    void SetSourceTime(const double t);

    /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
    void SetParameters(const Vector &u);

    virtual ~ConductionOperator();
};

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

    CAROM::Matrix* snapshot_mat = new CAROM::Matrix(snapshots[0]->dim(), snapshots.size(),
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

double InitialTemperature(const Vector &x);
double TimeWindowFunction(const double t, const double t_begin,
                          const double t_end);
double Amplitude(const double t, const int index);
double SourceFunction(const Vector &x, double t);

Vector bb_min, bb_max; // Mesh bounding box
double amp_in = 0.2;
double t_end_in = 0.1;
double amp_out = 0.1;
double t_end_out = 0.3;
double dt = 1.0e-2;

int main(int argc, char *argv[])
{
    // 1. Initialize MPI.
    int num_procs, myid;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // 2. Parse command-line options.
    const char *mesh_file = "";
    int ser_ref_levels = 2;
    int par_ref_levels = 1;
    int order = 2;
    int ode_solver_type = 1;
    double t_final = 0.5;
    double alpha = 1.0e-2;
    double kappa = 0.5;
//    double ef = 0.9999;
    double ef = -1;
    int rdim = -1;
    bool visualization = true;
    bool visit = false;
    int vis_steps = 5;
    bool adios2 = false;
    bool offline = false;
    bool online = false;

    double closest_rbf_val = 0.9;
    bool predict = false;
    bool save_dofs = false;
    bool csvFormat = true;
    const char *basename = "";

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
    args.AddOption(&amp_in, "-amp-in", "--amplitude-in",
                   "Amplitude of inlet source at (0,0).");
    args.AddOption(&amp_out, "-amp-out", "--amplitude-out",
                   "Amplitude of outlet source at (0.5,0.5).");
    args.AddOption(&t_end_in, "-t-end-in", "--t-end-in",
                   "End time of inlet source at (0,0).");
    args.AddOption(&t_end_out, "-t-end-out", "--t-end-out",
                   "End time of outlet source at (0.5,0.5).");
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
                   "Energy fraction for DMDc.");
    args.AddOption(&rdim, "-rdim", "--rdim",
                   "Reduced dimension for DMDc.");

    args.AddOption(&offline, "-offline", "--offline", "-no-offline", "--no-offline",
                   "Enable or disable the offline phase.");
    args.AddOption(&online, "-online", "--online", "-no-online", "--no-online",
                   "Enable or disable the online phase.");
    args.AddOption(&closest_rbf_val, "-crv", "--crv",
                   "DMD Closest RBF Value.");
    args.AddOption(&predict, "-predict", "--predict", "-no-predict", "--no-predict",
                   "Enable or disable DMD prediction.");
    args.AddOption(&save_dofs, "-save", "--save", "-no-save", "--no-save",
                   "Enable or disable MFEM DOF solution snapshot files).");
    args.AddOption(&csvFormat, "-csv", "--csv", "-hdf", "--hdf",
                   "Use CSV or HDF format for files output by -save option.");
    args.AddOption(&basename, "-out", "--outputfile-name",
                   "Name of the sub-folder to dump files within the run directory.");

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

    // 3. Read the serial mesh from the given mesh file on all processors. We can
    //    handle triangular, quadrilateral, tetrahedral and hexahedral meshes
    //    with the same code.
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

    // 4. Define the ODE solver used for time integration. Several implicit
    //    singly diagonal implicit Runge-Kutta (SDIRK) methods, as well as
    //    explicit Runge-Kutta methods are available.
    ODESolver *ode_solver;
    bool explicit_solver;
    switch (ode_solver_type)
    {
    // Implicit L-stable methods
    case 1:
        ode_solver = new BackwardEulerSolver;
        explicit_solver = false;
        break;
    case 2:
        ode_solver = new SDIRK23Solver(2);
        explicit_solver = false;
        break;
    case 3:
        ode_solver = new SDIRK33Solver;
        explicit_solver = false;
        break;
    // Explicit methods
    case 11:
        ode_solver = new ForwardEulerSolver;
        explicit_solver = true;
        break;
    case 12:
        ode_solver = new RK2Solver(0.5);
        explicit_solver = true;
        break; // midpoint method
    case 13:
        ode_solver = new RK3SSPSolver;
        explicit_solver = true;
        break;
    case 14:
        ode_solver = new RK4Solver;
        explicit_solver = true;
        break;
    case 15:
        ode_solver = new GeneralizedAlphaSolver(0.5);
        explicit_solver = true;
        break;
    // Implicit A-stable methods (not L-stable)
    case 22:
        ode_solver = new ImplicitMidpointSolver;
        explicit_solver = false;
        break;
    case 23:
        ode_solver = new SDIRK23Solver;
        explicit_solver = false;
        break;
    case 24:
        ode_solver = new SDIRK34Solver;
        explicit_solver = false;
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
    FunctionCoefficient s(SourceFunction);
    ConductionOperator oper(fespace, s, alpha, kappa, u);

    u_gf.SetFromTrueDofs(u);
    {
        ostringstream mesh_name, sol_name;
        mesh_name << outputPath << "/parametric_dmdc_heat_conduction_" << to_string(
                      alpha) << "_" << to_string(kappa) << "-mesh." << setfill('0') << setw(
                      6) << myid;
        sol_name << outputPath << "/parametric_dmdc_heat_conduction_" << to_string(
                     alpha) << "_" << to_string(kappa) <<  "-init." << setfill('0') << setw(
                     6) << myid;
        ofstream omesh(mesh_name.str().c_str());
        omesh.precision(precision);
        pmesh->Print(omesh);
        ofstream osol(sol_name.str().c_str());
        osol.precision(precision);
        u_gf.Save(osol);
    }
//    {
//        ostringstream mesh_name, sol_name;
//        mesh_name << "heat_conduction_dmdc-mesh." << setfill('0') << setw(6) << myid;
//        sol_name << "heat_conduction_dmdc-init." << setfill('0') << setw(6) << myid;
//        ofstream omesh(mesh_name.str().c_str());
//        omesh.precision(precision);
//        pmesh->Print(omesh);
//        ofstream osol(sol_name.str().c_str());
//        osol.precision(precision);
//        u_gf.Save(osol);
//    }


//    VisItDataCollection visit_dc("Heat_Conduction", pmesh);
    VisItDataCollection visit_dc(outputPath + "/parametric_dmdc_Heat_Conduction_" +
                                 to_string(alpha) + "_" + to_string(kappa), pmesh);
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
        const std::string collection_name = "heat_conduction_dmdc-p-" + postfix + ".bp";

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
    double f[2];

    CAROM::Vector* init = NULL;

    CAROM::Database *db = NULL;
    if (csvFormat)
        db = new CAROM::CSVDatabase();
    else
        db = new CAROM::HDFDatabase();

    vector<int> snap_list;

    fom_timer.Stop();

//    dmd_training_timer.Start();

    CAROM::DMDc* dmd_u = NULL;
    
    
    f[0] = Amplitude(t, 0);
    f[1] = Amplitude(t, 1);
    int dim_c = 2;
    
    std::vector<CAROM::Vector*> d_controls;
//    CAROM::Vector* d_controls;
//    Vector* control = new Vector(f, dim_c, false);
    CAROM::Vector* control = new CAROM::Vector(f, dim_c, false);
//    const std::vector control = new Vector(f, dim_c);
    d_controls.push_back(control);
    
//    std::cout << "actual control type: " << typeid(control).name() << std::endl;
//    std::cout << "actual d_control type: " << typeid(d_controls).name() << std::endl;

    if (offline)
    {
        dmd_training_timer.Start();

        // 11. Create DMDc object and take initial sample.
        u_gf.SetFromTrueDofs(u);
        f[0] = Amplitude(t, 0);
        f[1] = Amplitude(t, 1);

        std::cout << "t =  " << t << std::endl;

//        CAROM::DMDc* dmd_u; //redefined this on 494
        dmd_u = new CAROM::DMDc(u.Size(), 2, dt);
        dmd_u->takeSample(u.GetData(), t, f, false);
//        ts.push_back(t);

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

        // Assuming Euler method is used
        if (explicit_solver)
        {
            oper.SetSourceTime(t);
        }
        else
        {
            oper.SetSourceTime(t + dt);
        }
        ode_solver->Step(u, t, dt);

        fom_timer.Stop();
        
        
        f[0] = Amplitude(t, 0);
        f[1] = Amplitude(t, 1);
        if (offline)
        {
            dmd_training_timer.Start();

            u_gf.SetFromTrueDofs(u);
            dmd_u->takeSample(u.GetData(), t, f, last_step);
//            ts.push_back(t);

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
        if (!last_step)
        {
            CAROM::Vector* control = new CAROM::Vector(f, dim_c, false);
//            std::cout << "dimcon =  " << control->dim() << std::endl;
            d_controls.push_back(control);
        }

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

    // 12. Save the final solution in parallel. This output can be viewed later
    //     using GLVis: "glvis -np <np> -m heat_conduction_dmdc-mesh -g heat_conduction_dmdc-final".
    {
        ostringstream sol_name;
        sol_name << outputPath << "parametric_dmdc_Heat_Conduction" << to_string(
                     alpha) << "_" << to_string(kappa) << "-final." << setfill('0') << setw(
                     6) << myid;
        ofstream osol(sol_name.str().c_str());
        osol.precision(precision);
        u_gf.Save(osol);
//        ostringstream sol_name;
//        sol_name << "heat_conduction_dmdc-final." << setfill('0') << setw(6) << myid;
//        ofstream osol(sol_name.str().c_str());
//        osol.precision(precision);
//        u_gf.Save(osol);
    }

    // 13. Calculate the DMDc modes.
    // pretty sure it is easier to get rid of EF because each snapshot must have same rdim
    if (offline || online)
    {
        if (offline)
        {
            if (myid == 0 && rdim != -1 && ef != -1)
            {
                std::cout << "Both rdim and ef are set. ef will be ignored." << std::endl;
            }

            dmd_training_timer.Start();

            if (rdim != -1)
            {
                if (myid == 0)
                {
                    std::cout << "Creating DMDc with rdim: " << rdim << std::endl;
                }

                dmd_u->train(rdim);

            }
            else if (ef != -1)
            {
                if (myid == 0)
                {
                    std::cout << "Creating DMDc with energy fraction: " << ef << std::endl;
                }
                dmd_u->train(ef);
            }

            dmd_training_timer.Stop();

            dmd_u->save(outputPath + "/" + to_string(alpha) + "_" + to_string(
                            kappa));

            if (myid == 0)
            {
                std::ofstream fout;
                fout.open("parameters.txt", std::ios::app);
                fout << alpha << " " << kappa << std::endl;
                fout.close();
            }
        }

//         probably will need to change a lot here
        if (online)
        {
            if (myid == 0)
            {
                std::cout << "Creating DMD using the rdim of the offline phase" << std::endl;
            }

            std::fstream fin("parameters.txt", std::ios_base::in);
            double curr_param;
            std::vector<std::string> dmdc_paths;
            std::vector<CAROM::Vector*> param_vectors;

            while (fin >> curr_param)
            {
                double curr_alpha = curr_param;
                fin >> curr_param;
                double curr_kappa = curr_param;

                dmdc_paths.push_back(outputPath + "/" + to_string(curr_alpha) + "_" + to_string(
                                        curr_kappa) );
                CAROM::Vector* param_vector = new CAROM::Vector(2, false);
                param_vector->item(0) = curr_alpha;
                param_vector->item(1) = curr_kappa;
                param_vectors.push_back(param_vector);
            }
            fin.close();
            
//            int np = 5;
//            string temp;
//
//                 scanf("%d", &np);
//                 vector <int> money;
//                 vector <string> names;
//                 cin.ignore();
//                    for(int i = 0; i< np; i++){
//                      getline(cin, temp);
//                      names.push_back(temp);
//                      cout << dmdc_paths[i] << endl;
//                   }
//
//            for (int i = 0; i< 10; i++) {
//                std::cout << "d_control dim =  " << d_controls[i]->dim() << std::endl;
//            }

            CAROM::Vector* desired_param = new CAROM::Vector(2, false);
            desired_param->item(0) = alpha;
            desired_param->item(1) = kappa;

            dmd_training_timer.Start();

            std::cout << "here1" << std::endl;
            CAROM::getParametricDMDc(dmd_u, param_vectors, dmdc_paths, desired_param,
                                    "G", "LS", closest_rbf_val);
            std::cout << "here4" << std::endl;
            
            
            const CAROM::Matrix* f_controls = createControlMatrix(d_controls);
            std::cout << "here5" << std::endl;
            dmd_u->project(init,f_controls);
            std::cout << "here6" << std::endl;
//            dmd_u->projectInitialCondition(init)

            dmd_training_timer.Stop();
            delete desired_param;
        }

        if (predict)
        {
            Vector true_solution_u(u.Size());
            true_solution_u = u.GetData();

            dmd_prediction_timer.Start();

            // 14. Predict the state at t_final using DMDc.
            if (myid == 0)
            {
                std::cout << "Predicting temperature using DMDc" << std::endl;
            }

            CAROM::Vector* result_u = dmd_u->predict(ts[0]);
            Vector initial_dmd_solution_u(result_u->getData(), result_u->dim());
            u_gf.SetFromTrueDofs(initial_dmd_solution_u);

            VisItDataCollection dmd_visit_dc("parametric_dmdc_Heat_Conduction", pmesh);
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

            // 15. Calculate the relative error between the DMDc final solution and the true solution.
            Vector dmd_solution_u(result_u->getData(), result_u->dim());
            Vector diff_u(true_solution_u.Size());
            subtract(dmd_solution_u, true_solution_u, diff_u);

            double tot_diff_norm_u = sqrt(InnerProduct(MPI_COMM_WORLD, diff_u, diff_u));
            double tot_true_solution_u_norm = sqrt(InnerProduct(MPI_COMM_WORLD,
                                                   true_solution_u, true_solution_u));

            if (myid == 0)
            {

                std::cout << "Relative error of DMDc temperature (u) at t_final: " << t_final <<
                          " is " << tot_diff_norm_u / tot_true_solution_u_norm << std::endl;
                printf("Elapsed time for solving FOM: %e second\n", fom_timer.RealTime());
                printf("Elapsed time for training DMDc: %e second\n",
                       dmd_training_timer.RealTime());
                printf("Elapsed time for predicting DMDc: %e second\n",
                       dmd_prediction_timer.RealTime());
            }

            delete result_u;

        }

    }

    // 16. Free the used memory.
    delete ode_solver;
    delete pmesh;
//    delete result_u;
    if (offline)
    {
        delete dmd_u;
    }

    MPI_Finalize();

    return 0;
}

ConductionOperator::ConductionOperator(ParFiniteElementSpace &f,
                                       FunctionCoefficient &s,
                                       double al, double kap, const Vector &u)
    : TimeDependentOperator(f.GetTrueVSize(), 0.0), fespace(f), src_func(s),
      M(NULL), K(NULL), T(NULL), current_dt(0.0),
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

    B = new ParLinearForm(&fespace);
    B->AddDomainIntegrator(new DomainLFIntegrator(src_func));
}

void ConductionOperator::Mult(const Vector &u, Vector &du_dt) const
{
    // Compute:
    //    du_dt = M^{-1}*-K(u)
    // for du_dt
    Kmat.Mult(u, z);
    z.Neg(); // z = -z
    z += *b;
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
    z += *b;
    T_solver.Mult(z, du_dt);
}

void ConductionOperator::SetSourceTime(double t)
{
    src_func.SetTime(t);
    B->Assemble();
    b = B->ParallelAssemble();
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
    delete B;
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

double TimeWindowFunction(const double t, const double t_begin,
                          const double t_end)
{
    return 0.5 * (tanh((t - t_begin) / (5*dt)) - tanh((t - t_end) / (5*dt)));
}

double Amplitude(const double t, const int index)
{
    if (index == 0)
    {
        return amp_in * TimeWindowFunction(t, 0.0, t_end_in);
    }
    else
    {
        return amp_out * TimeWindowFunction(t, 0.0, t_end_out);
    }
}

double SourceFunction(const Vector &x, double t)
{
    // map to the reference [-1,1] domain
    Vector X(x.Size()), Y(x.Size());
    for (int i = 0; i < x.Size(); i++)
    {
        double center = (bb_min[i] + bb_max[i]) * 0.5;
        X(i) = 2 * (x(i) - center) / (bb_max[i] - bb_min[i]);
        Y(i) = X(i) - 0.5;
    }

    double r1 = X.Norml2() / 0.01;
    double r2 = Y.Norml2() / 0.01;
    return Amplitude(t,0) * exp(-0.5*r1*r1) - Amplitude(t,1) * exp(-0.5*r2*r2);
}
