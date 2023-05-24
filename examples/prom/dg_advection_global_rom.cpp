/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

//                       libROM MFEM Example: DG Advection (adapted from ex9p.cpp)
//
// Compile with: make dg_advection_global_rom
//
// For ROM (reproductive case):
//    dg_advection_global_rom -offline
//    dg_advection_global_rom -merge -ns 1
//    dg_advection_global_rom -online
//
// For ROM (global rom):
// Offline phase: dg_advection_global_rom -offline -ff 1.0 -id 0
//                dg_advection_global_rom -offline -ff 1.1 -id 1
//                dg_advection_global_rom -offline -ff 1.2 -id 2
//
// Merge phase:   dg_advection_global_rom -merge -ns 3
//
// FOM solution:  dg_advection_global_rom -fom -ff 1.15
//
// Online phase:  dg_advection_global_rom -online -ff 1.15
//
// Sample runs:
//    mpirun -np 4 dg_advection_global_rom -p 0 -dt 0.005
//    mpirun -np 4 dg_advection_global_rom -p 0 -dt 0.01
//    mpirun -np 4 dg_advection_global_rom -p 1 -dt 0.005 -tf 9
//    mpirun -np 4 dg_advection_global_rom -p 1 -rp 1 -dt 0.002 -tf 9
//    mpirun -np 4 dg_advection_global_rom -p 1 -rp 1 -dt 0.02 -s 13 -tf 9
//    mpirun -np 4 dg_advection_global_rom -p 1 -rp 1 -dt 0.004 -tf 9
//    mpirun -np 4 dg_advection_global_rom -p 1 -rp 1 -dt 0.005 -tf 9
//    mpirun -np 4 dg_advection_global_rom -p 3 -rp 2 -dt 0.0025 -tf 9 -vs 20
//    mpirun -np 4 dg_advection_global_rom -p 0 -o 2 -rp 1 -dt 0.01 -tf 8
//    mpirun -np 4 dg_advection_global_rom -p 0 -rs 2 -dt 0.005 -tf 2
//    mpirun -np 4 dg_advection_global_rom -p 0 -rs 1 -o 2 -tf 2
//    mpirun -np 3 dg_advection_global_rom -p 1 -rs 1 -rp 0 -dt 0.005 -tf 0.5
//
// Description:  This example code solves the time-dependent advection equation
//               du/dt + v.grad(u) = 0, where v is a given fluid velocity, and
//               u0(x)=u(0,x) is a given initial condition.
//
//               The example demonstrates the use of Discontinuous Galerkin (DG)
//               bilinear forms in MFEM (face integrators), the use of implicit
//               and explicit ODE time integrators, the definition of periodic
//               boundary conditions through periodic meshes, as well as the use
//               of GLVis for persistent visualization of a time-evolving
//               solution. Saving of time-dependent data files for visualization
//               with VisIt (visit.llnl.gov) and ParaView (paraview.org), as
//               well as the optional saving with ADIOS2 (adios2.readthedocs.io)
//               are also illustrated.
//
// This example runs in parallel with MPI, by using the same number of MPI ranks
// in all phases (offline, merge, fom, online).

#include "mfem.hpp"
#include "linalg/Vector.h"
#include "linalg/BasisGenerator.h"
#include "linalg/BasisReader.h"
#include "mfem/Utilities.hpp"
#include <cmath>
#include <set>
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// Choice for the problem setup. The fluid velocity, initial condition and
// inflow boundary condition are chosen based on this parameter.
int problem;
double f_factor;

// Velocity coefficient
void velocity_function(const Vector &x, Vector &v);

// Initial condition
double u0_function(const Vector &x);

// Inflow boundary condition
double inflow_function(const Vector &x);

// Mesh bounding box
Vector bb_min, bb_max;

// Type of preconditioner for implicit time integrator
enum class PrecType : int
{
    ILU = 0,
    AIR = 1
};

#if MFEM_HYPRE_VERSION >= 21800
// Algebraic multigrid preconditioner for advective problems based on
// approximate ideal restriction (AIR). Most effective when matrix is
// first scaled by DG block inverse, and AIR applied to scaled matrix.
// See https://doi.org/10.1137/17M1144350.
class AIR_prec : public Solver
{
private:
    const HypreParMatrix *A;
    // Copy of A scaled by block-diagonal inverse
    HypreParMatrix A_s;

    HypreBoomerAMG *AIR_solver;
    int blocksize;

public:
    AIR_prec(int blocksize_) : AIR_solver(NULL), blocksize(blocksize_) { }

    void SetOperator(const Operator &op)
    {
        width = op.Width();
        height = op.Height();

        A = dynamic_cast<const HypreParMatrix *>(&op);
        MFEM_VERIFY(A != NULL, "AIR_prec requires a HypreParMatrix.")

        // Scale A by block-diagonal inverse
        BlockInverseScale(A, &A_s, NULL, NULL, blocksize,
                          BlockInverseScaleJob::MATRIX_ONLY);
        delete AIR_solver;
        AIR_solver = new HypreBoomerAMG(A_s);
        AIR_solver->SetAdvectiveOptions(1, "", "FA");
        AIR_solver->SetPrintLevel(
            0);    // 6. Define the parallel mesh by a partitioning of the serial mesh. Refine
        //    this mesh further in parallel to increase the resolution. Once the
        //    parallel mesh is defined, the serial mesh can be deleted.
        AIR_solver->SetMaxLevels(50);
    }

    virtual void Mult(const Vector &x, Vector &y) const
    {
        // Scale the rhs by block inverse and solve system
        HypreParVector z_s;
        BlockInverseScale(A, NULL, &x, &z_s, blocksize,
                          BlockInverseScaleJob::RHS_ONLY);
        AIR_solver->Mult(z_s, y);
    }

    ~AIR_prec()
    {
        delete AIR_solver;
    }
};
#endif

class DG_Solver : public Solver
{
private:
    HypreParMatrix &M, &K;
    SparseMatrix M_diag;
    HypreParMatrix *A;
    GMRESSolver linear_solver;
    Solver *prec;
    double dt;
public:
    DG_Solver(HypreParMatrix &M_, HypreParMatrix &K_, const FiniteElementSpace &fes,
              PrecType prec_type)
        : M(M_),
          K(K_),
          A(NULL),
          linear_solver(M.GetComm()),
          dt(-1.0)
    {
        int block_size = fes.GetFE(0)->GetDof();
        if (prec_type == PrecType::ILU)
        {
            prec = new BlockILU(block_size,
                                BlockILU::Reordering::MINIMUM_DISCARDED_FILL);
        }
        else if (prec_type == PrecType::AIR)
        {
#if MFEM_HYPRE_VERSION >= 21800
            prec = new AIR_prec(block_size);
#else
            MFEM_ABORT("Must have MFEM_HYPRE_VERSION >= 21800 to use AIR.\n");
#endif
        }
        linear_solver.iterative_mode = false;
        linear_solver.SetRelTol(1e-9);
        linear_solver.SetAbsTol(0.0);
        linear_solver.SetMaxIter(100);
        linear_solver.SetPrintLevel(0);
        linear_solver.SetPreconditioner(*prec);

        M.GetDiag(M_diag);
    }

    void SetTimeStep(double dt_)
    {
        if (dt_ != dt)
        {
            dt = dt_;
            // Form operator A = M - dt*K
            delete A;
            A = Add(-dt, K, 0.0, K);
            SparseMatrix A_diag;
            A->GetDiag(A_diag);
            A_diag.Add(1.0, M_diag);
            // this will also call SetOperator on the preconditioner
            linear_solver.SetOperator(*A);
        }
    }

    void SetOperator(const Operator &op)
    {
        linear_solver.SetOperator(op);
    }

    virtual void Mult(const Vector &x, Vector &y) const
    {
        linear_solver.Mult(x, y);
    }

    ~DG_Solver()
    {
        delete prec;
        delete A;
    }
};

/** A time-dependent operator for the right-hand side of the ODE. The DG weak
    form of du/dt = -v.grad(u) is M du/dt = K u + b, where M and K are the mass
    and advection matrices, and b describes the flow on the boundary. This can
    be written as a general ODE, du/dt = M^{-1} (K u + b), and this class is
    used to evaluate the right-hand side. */
class ROM_FE_Evolution : public TimeDependentOperator
{
private:
    DenseMatrix *M, *K, *M_inv, *A_inv;
    Vector *b, *u_init_hat;
    mutable Vector z;

public:
    ROM_FE_Evolution(DenseMatrix* M_, DenseMatrix* K_, Vector* b_,
                     Vector* u_init_hat_, int num_cols);

    virtual void Mult(const Vector &x, Vector &y) const;
    virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k);

    virtual ~ROM_FE_Evolution();
};

/** A time-dependent operator for the right-hand side of the ODE. The DG weak
    form of du/dt = -v.grad(u) is M du/dt = K u + b, where M and K are the mass
    and advection matrices, and b describes the flow on the boundary. This can
    be written as a general ODE, du/dt = M^{-1} (K u + b), and this class is
    used to evaluate the right-hand side. */
class FE_Evolution : public TimeDependentOperator
{
private:
    OperatorHandle M, K;
    const Vector &b;
    Solver *M_prec;
    CGSolver M_solver;
    DG_Solver *dg_solver;

    mutable Vector z;

public:
    FE_Evolution(ParBilinearForm &M_, ParBilinearForm &K_, const Vector &b_,
                 PrecType prec_type);

    virtual void Mult(const Vector &x, Vector &y) const;
    virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k);

    virtual ~FE_Evolution();
};

int main(int argc, char *argv[])
{
    // 1. Initialize MPI.
    MPI_Session mpi;
    int num_procs = mpi.WorldSize();
    int myid = mpi.WorldRank();

    // 2. Parse command-line options.
    problem = 3;
    const char *mesh_file = "../data/periodic-hexagon.mesh";
    int ser_ref_levels = 2;
    int par_ref_levels = 0;
    int order = 3;
    bool pa = false;
    bool ea = false;
    bool fa = false;
    const char *device_config = "cpu";
    int ode_solver_type = 11;
    double t_final = 10.0;
    double dt = 0.01;
    double ef = 0.9999;
    f_factor = 1.0;
    int rdim = -1;
    int id = 0;
    int nsets = 0;
    bool fom = false;
    bool offline = false;
    bool merge = false;
    bool online = false;
    bool visualization = false;
    bool visit = false;
    bool paraview = false;
    bool adios2 = false;
    bool binary = false;
    int vis_steps = 5;
#if MFEM_HYPRE_VERSION >= 21800
    PrecType prec_type = PrecType::AIR;
#else
    PrecType prec_type = PrecType::ILU;
#endif
    int precision = 16;
    cout.precision(precision);

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&problem, "-p", "--problem",
                   "Problem setup to use. See options in velocity_function().");
    args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                   "Number of times to refine the mesh uniformly in serial.");
    args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                   "Number of times to refine the mesh uniformly in parallel.");
    args.AddOption(&order, "-o", "--order",
                   "Order (degree) of the finite elements.");
    args.AddOption(&id, "-id", "--id", "Parametric id");
    args.AddOption(&nsets, "-ns", "--nset", "Number of parametric snapshot sets");
    args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                   "--no-partial-assembly", "Enable Partial Assembly.");
    args.AddOption(&ea, "-ea", "--element-assembly", "-no-ea",
                   "--no-element-assembly", "Enable Element Assembly.");
    args.AddOption(&fa, "-fa", "--full-assembly", "-no-fa",
                   "--no-full-assembly", "Enable Full Assembly.");
    args.AddOption(&device_config, "-d", "--device",
                   "Device configuration string, see Device::Configure().");
    args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                   "ODE solver: 1 - Forward Euler,\n\t"
                   "            2 - RK2 SSP, 3 - RK3 SSP, 4 - RK4, 6 - RK6,\n\t"
                   "            11 - Backward Euler,\n\t"
                   "            12 - SDIRK23 (L-stable), 13 - SDIRK33,\n\t"
                   "            22 - Implicit Midpoint Method,\n\t"
                   "            23 - SDIRK23 (A-stable), 24 - SDIRK34");
    args.AddOption(&t_final, "-tf", "--t-final",
                   "Final time; start time is 0.");
    args.AddOption(&dt, "-dt", "--time-step",
                   "Time step.");
    args.AddOption(&f_factor, "-ff", "--f-factor",
                   "Frequency scalar factor.");
    args.AddOption((int *)&prec_type, "-pt", "--prec-type", "Preconditioner for "
                   "implicit solves. 0 for ILU, 1 for pAIR-AMG.");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                   "--no-visualization",
                   "Enable or disable GLVis visualization.");
    args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                   "--no-visit-datafiles",
                   "Save data files for VisIt (visit.llnl.gov) visualization.");
    args.AddOption(&paraview, "-paraview", "--paraview-datafiles", "-no-paraview",
                   "--no-paraview-datafiles",
                   "Save data files for ParaView (paraview.org) visualization.");
    args.AddOption(&adios2, "-adios2", "--adios2-streams", "-no-adios2",
                   "--no-adios2-streams",
                   "Save data using adios2 streams.");
    args.AddOption(&binary, "-binary", "--binary-datafiles", "-ascii",
                   "--ascii-datafiles",
                   "Use binary (Sidre) or ascii format for VisIt data files.");
    args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                   "Visualize every n-th timestep.");
    args.AddOption(&fom, "-fom", "--fom", "-no-fom", "--no-fom",
                   "Enable or disable the fom phase.");
    args.AddOption(&offline, "-offline", "--offline", "-no-offline", "--no-offline",
                   "Enable or disable the offline phase.");
    args.AddOption(&online, "-online", "--online", "-no-online", "--no-online",
                   "Enable or disable the online phase.");
    args.AddOption(&merge, "-merge", "--merge", "-no-merge", "--no-merge",
                   "Enable or disable the merge phase.");
    args.AddOption(&ef, "-ef", "--energy_fraction",
                   "Energy fraction.");
    args.AddOption(&rdim, "-rdim", "--rdim",
                   "Reduced dimension.");
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

    if (fom)
    {
        MFEM_VERIFY(fom && !offline && !online
                    && !merge, "everything must be turned off if fom is used.");
    }
    else
    {
        bool check = (offline && !merge && !online) || (!offline && merge && !online)
                     || (!offline && !merge && online);
        MFEM_VERIFY(check, "only one of offline, merge, or online must be true!");
    }

    Device device(device_config);
    if (mpi.Root()) {
        device.Print();
    }

    // 3. Read the serial mesh from the given mesh file on all processors. We can
    //    handle geometrically periodic meshes in this code.
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    // 4. Define the ODE solver used for time integration. Several explicit
    //    Runge-Kutta methods are available.
    ODESolver *ode_solver = NULL;
    switch (ode_solver_type)
    {
    // Explicit methods
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
    // Implicit (L-stable) methods
    case 11:
        ode_solver = new BackwardEulerSolver;
        break;
    case 12:
        ode_solver = new SDIRK23Solver(2);
        break;
    case 13:
        ode_solver = new SDIRK33Solver;
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
        if (mpi.Root())
        {
            cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
        }
        delete mesh;    // 6. Define the parallel mesh by a partitioning of the serial mesh. Refine
        //    this mesh further in parallel to increase the resolution. Once the
        //    parallel mesh is defined, the serial mesh can be deleted.
        return 3;
    }

    // 5. Refine the mesh in serial to increase the resolution. In this example
    //    we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
    //    a command-line parameter. If the mesh is of NURBS type, we convert it
    //    to a (piecewise-polynomial) high-order mesh.
    for (int lev = 0; lev < ser_ref_levels; lev++)
    {
        mesh->UniformRefinement();
    }
    if (mesh->NURBSext)
    {
        mesh->SetCurvature(max(order, 1));
    }
    mesh->GetBoundingBox(bb_min, bb_max, max(order, 1));

    // 6. Define the parallel mesh by a partitioning of the serial mesh. Refine
    //    this mesh further in parallel to increase the resolution. Once the
    //    parallel mesh is defined, the serial mesh can be deleted.
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
    for (int lev = 0; lev < par_ref_levels; lev++)
    {
        pmesh->UniformRefinement();
    }

    // 7. Define the parallel discontinuous DG finite element space on the
    //    parallel refined mesh of the given polynomial order.
    DG_FECollection fec(order, dim, BasisType::GaussLobatto);
    ParFiniteElementSpace *fes = new ParFiniteElementSpace(pmesh, &fec);

    HYPRE_BigInt global_vSize = fes->GlobalTrueVSize();
    if (mpi.Root())
    {
        cout << "Number of unknowns: " << global_vSize << endl;
    }

    // 8. Set up and assemble the parallel bilinear and linear forms (and the
    //    parallel hypre matrices) corresponding to the DG discretization. The
    //    DGTraceIntegrator involves integrals over mesh interior faces.
    VectorFunctionCoefficient velocity(dim, velocity_function);
    FunctionCoefficient inflow(inflow_function);
    FunctionCoefficient u0(u0_function);
    StopWatch solveTimer, assembleTimer, mergeTimer;

    assembleTimer.Start();
    ParBilinearForm *m = new ParBilinearForm(fes);
    ParBilinearForm *k = new ParBilinearForm(fes);
    if (pa)
    {
        m->SetAssemblyLevel(AssemblyLevel::PARTIAL);
        k->SetAssemblyLevel(AssemblyLevel::PARTIAL);
    }
    else if (ea)
    {
        m->SetAssemblyLevel(AssemblyLevel::ELEMENT);
        k->SetAssemblyLevel(AssemblyLevel::ELEMENT);
    }
    else if (fa)
    {
        m->SetAssemblyLevel(AssemblyLevel::FULL);
        k->SetAssemblyLevel(AssemblyLevel::FULL);
    }

    m->AddDomainIntegrator(new MassIntegrator);
    constexpr double alpha = -1.0;
    k->AddDomainIntegrator(new ConvectionIntegrator(velocity, alpha));
    k->AddInteriorFaceIntegrator(
        new NonconservativeDGTraceIntegrator(velocity, alpha));
    k->AddBdrFaceIntegrator(
        new NonconservativeDGTraceIntegrator(velocity, alpha));

    ParLinearForm *b = new ParLinearForm(fes);
    b->AddBdrFaceIntegrator(
        new BoundaryFlowIntegrator(inflow, velocity, alpha));

    int skip_zeros = 0;
    m->Assemble();
    k->Assemble(skip_zeros);
    b->Assemble();
    m->Finalize();
    k->Finalize(skip_zeros);

    HypreParVector *B = b->ParallelAssemble();

    // 9. Define the initial conditions, save the corresponding grid function to
    //    a file and (optionally) save data in the VisIt format and initialize
    //    GLVis visualization.
    ParGridFunction *u = new ParGridFunction(fes);
    u->ProjectCoefficient(u0);
    HypreParVector *U = u->GetTrueDofs();

    assembleTimer.Stop();

    {
        ostringstream mesh_name, sol_name;
        mesh_name << "dg_advection_global_rom-mesh." << setfill('0') << setw(6) << myid;
        sol_name << "dg_advection_global_rom-init." << setfill('0') << setw(6) << myid;
        ofstream omesh(mesh_name.str().c_str());
        omesh.precision(precision);
        pmesh->Print(omesh);
        ofstream osol(sol_name.str().c_str());
        osol.precision(precision);
        u->Save(osol);
    }

    // Create data collection for solution output: either VisItDataCollection for
    // ascii data files, or SidreDataCollection for binary data files.
    DataCollection *dc = NULL;
    if (visit)
    {
        const std::string visitSuffix = online ? "_rom" : "_fom";
        if (binary)
        {
#ifdef MFEM_USE_SIDRE
            dc = new SidreDataCollection("DG_Advection" + visitSuffix, pmesh);
#else
            MFEM_ABORT("Must build with MFEM_USE_SIDRE=YES for binary output.");
#endif
        }
        else
        {
            dc = new VisItDataCollection("DG_Advection" + visitSuffix, pmesh);
            dc->SetPrecision(precision);
            // To save the mesh using MFEM's parallel mesh format:
            // dc->SetFormat(DataCollection::PARALLEL_FORMAT);
        }
        dc->RegisterField("solution", u);
        dc->SetCycle(0);
        dc->SetTime(0.0);
        dc->Save();
    }

    ParaViewDataCollection *pd = NULL;
    if (paraview)
    {
        const std::string suffix = online ? "_rom" : "_fom";
        pd = new ParaViewDataCollection("DG_Advection" + suffix, pmesh);
        pd->SetPrefixPath("ParaView");
        pd->RegisterField("solution", u);
        pd->SetLevelsOfDetail(order);
        pd->SetDataFormat(VTKFormat::BINARY);
        pd->SetHighOrderOutput(true);
        pd->SetCycle(0);
        pd->SetTime(0.0);
        pd->Save();
    }

    // Optionally output a BP (binary pack) file using ADIOS2. This can be
    // visualized with the ParaView VTX reader.
#ifdef MFEM_USE_ADIOS2
    ADIOS2DataCollection *adios2_dc = NULL;
    if (adios2)
    {
        std::string postfix(mesh_file);
        postfix.erase(0, std::string("../data/").size() );
        postfix += "_o" + std::to_string(order);
        const std::string collection_name = "dg_advection_global_rom-p-" + postfix +
                                            ".bp";

        adios2_dc = new ADIOS2DataCollection(MPI_COMM_WORLD, collection_name, pmesh);
        // output data substreams are half the number of mpi processes
        adios2_dc->SetParameter("SubStreams", std::to_string(num_procs/2) );
        // adios2_dc->SetLevelsOfDetail(2);
        adios2_dc->RegisterField("solution", u);
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
        if (!sout)
        {
            if (mpi.Root())
                cout << "Unable to connect to GLVis server at "
                     << vishost << ':' << visport << endl;
            visualization = false;
            if (mpi.Root())
            {
                cout << "GLVis visualization disabled.\n";
            }
        }
        else
        {
            sout << "parallel " << num_procs << " " << myid << "\n";
            sout.precision(precision);
            sout << "solution\n" << *pmesh << *u;
            sout << "pause\n";
            sout << flush;
            if (mpi.Root())
                cout << "GLVis visualization paused."
                     << " Press space (in the GLVis window) to resume it.\n";
        }
    }

    double t = 0.0;

    int max_num_snapshots = 100000;
    bool update_right_SV = false;
    bool isIncremental = false;
    const std::string basisName = "basis";
    const std::string basisFileName = basisName + std::to_string(id);
    const CAROM::Matrix* spatialbasis;
    CAROM::Options* options;
    CAROM::BasisGenerator *generator;
    int numRowRB, numColumnRB;

    CAROM::Matrix *M_hat_carom, *K_hat_carom;
    DenseMatrix *M_hat, *K_hat;
    CAROM::Vector *b_hat_carom, *u_init_hat_carom;
    Vector *b_hat, *u_init_hat;

    Vector u_init(*U);
    Vector *u_hat;

    // 10. Set BasisGenerator if offline
    if (offline)
    {
        options = new CAROM::Options(U->Size(), max_num_snapshots, 1, update_right_SV);
        generator = new CAROM::BasisGenerator(*options, isIncremental, basisFileName);
        Vector u_curr(*U);
        Vector u_centered(U->Size());
        subtract(u_curr, u_init, u_centered);
        bool addSample = generator->takeSample(u_centered.GetData(), t, dt);
    }

    // 11. The merge phase
    if (merge)
    {
        mergeTimer.Start();
        std::unique_ptr<CAROM::BasisGenerator> basis_generator;
        options = new CAROM::Options(U->Size(), max_num_snapshots, 1, update_right_SV);
        generator = new CAROM::BasisGenerator(*options, isIncremental, basisName);
        for (int paramID=0; paramID<nsets; ++paramID)
        {
            std::string snapshot_filename = basisName + std::to_string(
                                                paramID) + "_snapshot";
            generator->loadSamples(snapshot_filename,"snapshot");
        }
        generator->endSamples(); // save the merged basis file
        mergeTimer.Stop();
        if (myid == 0)
        {
            printf("Elapsed time for merging and building ROM basis: %e second\n",
                   mergeTimer.RealTime());
        }
        delete generator;
        delete options;
        MPI_Finalize();
        return 0;
    }

    if (online)
    {
        assembleTimer.Start();
        CAROM::BasisReader reader(basisName);
        if (rdim != -1)
        {
            spatialbasis = reader.getSpatialBasis(0.0, rdim);
        }
        else
        {
            spatialbasis = reader.getSpatialBasis(0.0, ef);
        }
        numRowRB = spatialbasis->numRows();
        numColumnRB = spatialbasis->numColumns();
        if (myid == 0) printf("spatial basis dimension is %d x %d\n", numRowRB,
                                  numColumnRB);

        OperatorHandle M, K;

        if (m->GetAssemblyLevel()==AssemblyLevel::LEGACY)
        {
            M.Reset(m->ParallelAssemble(), true);
            K.Reset(k->ParallelAssemble(), true);
        }
        else
        {
            M.Reset(m, false);
            K.Reset(k, false);
        }

        HypreParMatrix &M_mat = *M.As<HypreParMatrix>();
        HypreParMatrix &K_mat = *K.As<HypreParMatrix>();

        M_hat_carom = new CAROM::Matrix(numRowRB, numColumnRB, false);
        ComputeCtAB(M_mat, *spatialbasis, *spatialbasis, *M_hat_carom);

        // libROM stores the matrix row-wise, so wrapping as a DenseMatrix in MFEM means it is transposed.
        M_hat = new DenseMatrix(numColumnRB, numColumnRB);
        M_hat->Set(1, M_hat_carom->getData());
        M_hat->Transpose();

        K_hat_carom = new CAROM::Matrix(numRowRB, numColumnRB, false);
        ComputeCtAB(K_mat, *spatialbasis, *spatialbasis, *K_hat_carom);

        // libROM stores the matrix row-wise, so wrapping as a DenseMatrix in MFEM means it is transposed.
        K_hat = new DenseMatrix(numColumnRB, numColumnRB);
        K_hat->Set(1, K_hat_carom->getData());
        K_hat->Transpose();

        Vector b_vec = *B;
        CAROM::Vector b_carom(b_vec.GetData(), b_vec.Size(), true);
        b_hat_carom = spatialbasis->transposeMult(&b_carom);
        b_hat = new Vector(b_hat_carom->getData(), b_hat_carom->dim());

        u_init_hat_carom = new CAROM::Vector(numColumnRB, false);
        ComputeCtAB_vec(K_mat, *U, *spatialbasis, *u_init_hat_carom);
        u_init_hat = new Vector(u_init_hat_carom->getData(), u_init_hat_carom->dim());

        u_hat = new Vector(numColumnRB);
        *u_hat = 0.0;
        assembleTimer.Stop();
    }

    TimeDependentOperator* adv;

    // 10. Define the time-dependent evolution operator describing the ODE
    //     right-hand side, and perform time-integration (looping over the time
    //     iterations, ti, with a time-step dt).
    if (online)
    {
        adv = new ROM_FE_Evolution(M_hat, K_hat, b_hat, u_init_hat, numColumnRB);
        adv->SetTime(t);
        ode_solver->Init(*adv);
    }
    else
    {
        adv = new FE_Evolution(*m, *k, *B, prec_type);
        adv->SetTime(t);
        ode_solver->Init(*adv);
    }

    bool done = false;
    for (int ti = 0; !done; )
    {
        double dt_real = min(dt, t_final - t);

        solveTimer.Start();
        if (online)
        {
            ode_solver->Step(*u_hat, t, dt_real);
        }
        else
        {
            ode_solver->Step(*U, t, dt_real);
        }
        solveTimer.Stop();

        ti++;
        done = (t >= t_final - 1e-8*dt);

        // 18. take and write snapshot for ROM
        if (offline)
        {
            Vector u_curr(*U);
            Vector u_centered(U->Size());
            subtract(u_curr, u_init, u_centered);
            bool addSample = generator->takeSample(u_centered.GetData(), t, dt);
        }

        if (done || ti % vis_steps == 0)
        {
            if (mpi.Root())
            {
                cout << "time step: " << ti << ", time: " << t << endl;
            }

            // 11. Extract the parallel grid function corresponding to the finite
            //     element approximation U (the local solution on each processor).
            *u = *U;

            if (visualization)
            {
                sout << "parallel " << num_procs << " " << myid << "\n";
                sout << "solution\n" << *pmesh << *u << flush;
            }

            if (visit)
            {
                if (online)
                {
                    if (mpi.Root())
                        cout << "WARNING: FOM lifting for visualization is slow." << endl;

                    CAROM::Vector u_hat_final_carom(u_hat->GetData(), u_hat->Size(), false);
                    CAROM::Vector* u_final_carom = spatialbasis->mult(u_hat_final_carom);
                    Vector u_final(u_final_carom->getData(), u_final_carom->dim());
                    u_final += u_init;
                    u->SetFromTrueDofs(u_final);
                    delete u_final_carom;
                }

                dc->SetCycle(ti);
                dc->SetTime(t);
                dc->Save();
            }

            if (paraview)
            {
                pd->SetCycle(ti);
                pd->SetTime(t);
                pd->Save();
            }

#ifdef MFEM_USE_ADIOS2
            // transient solutions can be visualized with ParaView
            if (adios2)
            {
                adios2_dc->SetCycle(ti);
                adios2_dc->SetTime(t);
                adios2_dc->Save();
            }
#endif
        }
    }

    if (offline)
    {
        generator->writeSnapshot();
        delete generator;
        delete options;
    }

    if (online)
    {
        CAROM::Vector u_hat_final_carom(u_hat->GetData(), u_hat->Size(), false);
        CAROM::Vector* u_final_carom = spatialbasis->mult(u_hat_final_carom);
        Vector u_final(u_final_carom->getData(), u_final_carom->dim());
        u_final += u_init;

        Vector fom_solution(U->Size());
        ifstream solution_file;
        ostringstream solution_filename;
        solution_filename << "dg_advection_global_rom-final." << f_factor << "." <<
                          setfill('0') << setw(6) << myid;
        solution_file.open(solution_filename.str());
        fom_solution.Load(solution_file, U->Size());
        const double fomNorm = sqrt(InnerProduct(MPI_COMM_WORLD, fom_solution,
                                    fom_solution));
        fom_solution -= u_final;
        const double diffNorm = sqrt(InnerProduct(MPI_COMM_WORLD, fom_solution,
                                     fom_solution));
        if (myid == 0) std::cout << "Relative l2 error of ROM solution " << diffNorm /
                                     fomNorm << std::endl;

        delete spatialbasis;
        delete M_hat_carom;
        delete K_hat_carom;
        delete M_hat;
        delete K_hat;
        delete b_hat_carom;
        delete u_init_hat_carom;
        delete b_hat;
        delete u_init_hat;
        delete u_hat;
        delete u_final_carom;
    }

    // 12. Save the final solution in parallel. This output can be viewed later
    //     using GLVis: "glvis -np <np> -m dg_advection_global_rom-mesh -g dg_advection_global_rom-final".
    if (offline || fom)
    {
        *u = *U;
        ostringstream sol_name;
        sol_name << "dg_advection_global_rom-final." << f_factor << "." << setfill('0')
                 << setw(6) << myid;
        ofstream osol(sol_name.str().c_str());
        osol.precision(precision);
        Vector tv(u->ParFESpace()->GetTrueVSize());
        u->GetTrueDofs(tv);

        for (int i=0; i<tv.Size(); ++i)
            osol << tv[i] << std::endl;

        osol.close();
    }

    // 13. print timing info
    if (myid == 0)
    {
        if (fom || offline)
        {
            printf("Elapsed time for assembling FOM: %e second\n",
                   assembleTimer.RealTime());
            printf("Elapsed time for solving FOM: %e second\n", solveTimer.RealTime());
        }
        if (online)
        {
            printf("Elapsed time for assembling ROM: %e second\n",
                   assembleTimer.RealTime());
            printf("Elapsed time for solving ROM: %e second\n", solveTimer.RealTime());
        }
    }

    // 16. Free the used memory.
    delete U;
    delete u;
    delete B;
    delete b;
    delete k;
    delete m;
    delete fes;
    delete pmesh;
    delete ode_solver;
    delete pd;
#ifdef MFEM_USE_ADIOS2
    if (adios2)
    {
        delete adios2_dc;
    }
#endif
    delete dc;

    return 0;
}

// Implementation of class ROM_FE_Evolution
ROM_FE_Evolution::ROM_FE_Evolution(DenseMatrix* M_, DenseMatrix* K_, Vector* b_,
                                   Vector* u_init_hat_, int num_cols)
    : TimeDependentOperator(num_cols),
      z(num_cols)
{
    M = M_;
    K = K_;
    b = b_;
    A_inv = NULL;
    M_inv = NULL;
    u_init_hat = u_init_hat_;
    M_inv = new DenseMatrix(*M_);
    M_inv->Invert();
}

// Solve the equation:
//    u_t = M_hat^{-1}(K_hatu + b_hat + u_init_hat),
// by solving associated linear system
//    (M_hat - dt*K_hat) d = K_hat*u + b_hat + u_init_hat
void ROM_FE_Evolution::ImplicitSolve(const double dt, const Vector &x,
                                     Vector &k)
{
    K->Mult(x, z);
    z += *b;
    z += *u_init_hat;

    // Assume dt is constant. Pre-compute A_inv.
    if (A_inv == NULL)
    {
        DenseMatrix A(K->NumRows(), K->NumCols());
        A.Set(dt, *K);
        A_inv = new DenseMatrix(*M);
        *A_inv -= A;
        A_inv->Invert();
    }

    A_inv->Mult(z, k);
}

void ROM_FE_Evolution::Mult(const Vector &x, Vector &y) const
{
    // y = M_hat^{-1} (K_hat x + b_hat + u_init_hat )
    K->Mult(x, z);
    z += *b;
    z += *u_init_hat;
    M_inv->Mult(z, y);
}

ROM_FE_Evolution::~ROM_FE_Evolution()
{
    delete M_inv;
    delete A_inv;
}

// Implementation of class FE_Evolution
FE_Evolution::FE_Evolution(ParBilinearForm &M_, ParBilinearForm &K_,
                           const Vector &b_, PrecType prec_type)
    : TimeDependentOperator(M_.Height()), b(b_),
      M_solver(M_.ParFESpace()->GetComm()),
      z(M_.Height())
{
    if (M_.GetAssemblyLevel()==AssemblyLevel::LEGACY)
    {
        M.Reset(M_.ParallelAssemble(), true);
        K.Reset(K_.ParallelAssemble(), true);
    }
    else
    {
        M.Reset(&M_, false);
        K.Reset(&K_, false);
    }

    M_solver.SetOperator(*M);

    Array<int> ess_tdof_list;
    if (M_.GetAssemblyLevel()==AssemblyLevel::LEGACY)
    {
        HypreParMatrix &M_mat = *M.As<HypreParMatrix>();
        HypreParMatrix &K_mat = *K.As<HypreParMatrix>();
        HypreSmoother *hypre_prec = new HypreSmoother(M_mat, HypreSmoother::Jacobi);
        M_prec = hypre_prec;

        dg_solver = new DG_Solver(M_mat, K_mat, *M_.FESpace(), prec_type);
    }
    else
    {
        M_prec = new OperatorJacobiSmoother(M_, ess_tdof_list);
        dg_solver = NULL;
    }

    M_solver.SetPreconditioner(*M_prec);
    M_solver.iterative_mode = false;
    M_solver.SetRelTol(1e-9);
    M_solver.SetAbsTol(0.0);
    M_solver.SetMaxIter(100);
    M_solver.SetPrintLevel(0);
}

// Solve the equation:
//    u_t = M^{-1}(Ku + b),
// by solving associated linear system
//    (M - dt*K) d = K*u + b
void FE_Evolution::ImplicitSolve(const double dt, const Vector &x, Vector &k)
{
    K->Mult(x, z);
    z += b;
    dg_solver->SetTimeStep(dt);
    dg_solver->Mult(z, k);
}

void FE_Evolution::Mult(const Vector &x, Vector &y) const
{
    // y = M^{-1} (K x + b)
    K->Mult(x, z);
    z += b;
    M_solver.Mult(z, y);
}

FE_Evolution::~FE_Evolution()
{
    delete M_prec;
    delete dg_solver;
}

// Velocity coefficient
void velocity_function(const Vector &x, Vector &v)
{
    int dim = x.Size();

    // map to the reference [-1,1] domain
    Vector X(dim);
    for (int i = 0; i < dim; i++)
    {
        double center = (bb_min[i] + bb_max[i]) * 0.5;
        X(i) = 2 * (x(i) - center) / (bb_max[i] - bb_min[i]);
    }

    switch (problem)
    {
    case 0:
    {
        // Translations in 1D, 2D, and 3D
        switch (dim)
        {
        case 1:
            v(0) = 1.0;
            break;
        case 2:
            v(0) = sqrt(2./3.);
            v(1) = sqrt(1./3.);
            break;
        case 3:
            v(0) = sqrt(3./6.);
            v(1) = sqrt(2./6.);
            v(2) = sqrt(1./6.);
            break;
        }
        break;
    }
    case 1:
    case 2:
    {
        // Clockwise rotation in 2D around the origin
        const double w = M_PI/2;
        switch (dim)
        {
        case 1:
            v(0) = 1.0;
            break;
        case 2:
            v(0) = w*X(1);
            v(1) = -w*X(0);
            break;
        case 3:
            v(0) = w*X(1);
            v(1) = -w*X(0);
            v(2) = 0.0;
            break;
        }
        break;
    }
    case 3:
    {
        // Clockwise twisting rotation in 2D around the origin
        const double w = M_PI/2;
        double d = max((X(0)+1.)*(1.-X(0)),0.) * max((X(1)+1.)*(1.-X(1)),0.);
        d = d*d;
        switch (dim)
        {
        case 1:
            v(0) = 1.0;
            break;
        case 2:
            v(0) = d*w*X(1);
            v(1) = -d*w*X(0);
            break;
        case 3:
            v(0) = d*w*X(1);
            v(1) = -d*w*X(0);
            v(2) = 0.0;
            break;
        }
        break;
    }
    }
}

// Initial condition
double u0_function(const Vector &x)
{
    int dim = x.Size();

    // map to the reference [-1,1] domain
    Vector X(dim);
    for (int i = 0; i < dim; i++)
    {
        double center = (bb_min[i] + bb_max[i]) * 0.5;
        X(i) = 2 * (x(i) - center) / (bb_max[i] - bb_min[i]);
    }

    switch (problem)
    {
    case 0:
    case 1:
    {
        switch (dim)
        {
        case 1:
            return exp(-40.*pow(X(0)-0.5,2));
        case 2:
        case 3:
        {
            double rx = 0.45, ry = 0.25, cx = 0., cy = -0.2, w = 10.;
            if (dim == 3)
            {
                const double s = (1. + 0.25*cos(2*M_PI*X(2)));
                rx *= s;
                ry *= s;
            }
            return ( erfc(w*(X(0)-cx-rx))*erfc(-w*(X(0)-cx+rx)) *
                     erfc(w*(X(1)-cy-ry))*erfc(-w*(X(1)-cy+ry)) )/16;
        }
        }
    }
    case 2:
    {
        double x_ = X(0), y_ = X(1), rho, phi;
        rho = hypot(x_, y_);
        phi = atan2(y_, x_);
        return pow(sin(M_PI*rho),2)*sin(3*phi);
    }
    case 3:
    {
        const double f = M_PI * f_factor;
        return sin(f*X(0))*sin(f*X(1));
    }
    }
    return 0.0;
}

// Inflow boundary condition (zero for the problems considered in this example)
double inflow_function(const Vector &x)
{
    switch (problem)
    {
    case 0:
    case 1:
    case 2:
    case 3:
        return 0.0;
    }
    return 0.0;
}
