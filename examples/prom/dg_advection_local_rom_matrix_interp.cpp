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
// Compile with: make dg_advection_local_rom_matrix_interp
//
// For ROM (reproductive case):
//    dg_advection_local_rom_matrix_interp -offline
//    dg_advection_local_rom_matrix_interp -online
//
// For ROM (parametric case using matrix interpolation):
//    rm -rf frequencies.txt
//    dg_advection_local_rom_matrix_interp -offline -ff 1.02
//    dg_advection_local_rom_matrix_interp -interp_prep -ff 1.02 -rdim 40
//    dg_advection_local_rom_matrix_interp -offline -ff 1.03
//    dg_advection_local_rom_matrix_interp -interp_prep -ff 1.03 -rdim 40
//    dg_advection_local_rom_matrix_interp -offline -ff 1.04
//    dg_advection_local_rom_matrix_interp -interp_prep -ff 1.04 -rdim 40
//    dg_advection_local_rom_matrix_interp -offline -ff 1.06
//    dg_advection_local_rom_matrix_interp -interp_prep -ff 1.06 -rdim 40
//    dg_advection_local_rom_matrix_interp -offline -ff 1.07
//    dg_advection_local_rom_matrix_interp -interp_prep -ff 1.07 -rdim 40
//    dg_advection_local_rom_matrix_interp -offline -ff 1.08
//    dg_advection_local_rom_matrix_interp -interp_prep -ff 1.08 -rdim 40
//    dg_advection_local_rom_matrix_interp -fom -ff 1.05
//    dg_advection_local_rom_matrix_interp -online_interp -ff 1.05 -rdim 40 (interpolate using a linear solve)
//    dg_advection_local_rom_matrix_interp -online_interp -ff 1.05 -rdim 40 -im "LP" (interpolate using lagragian polynomials)
//    dg_advection_local_rom_matrix_interp -online_interp -ff 1.05 -rdim 40 -im "IDW" (interpolate using inverse distance weighting)
//
// Sample runs:
//    mpirun -np 4 dg_advection_local_rom_matrix_interp -p 0 -dt 0.005
//    mpirun -np 4 dg_advection_local_rom_matrix_interp -p 0 -dt 0.01
//    mpirun -np 4 dg_advection_local_rom_matrix_interp -p 1 -dt 0.005 -tf 9
//    mpirun -np 4 dg_advection_local_rom_matrix_interp -p 1 -rp 1 -dt 0.002 -tf 9
//    mpirun -np 4 dg_advection_local_rom_matrix_interp -p 1 -rp 1 -dt 0.02 -s 13 -tf 9
//    mpirun -np 4 dg_advection_local_rom_matrix_interp -p 1 -rp 1 -dt 0.004 -tf 9
//    mpirun -np 4 dg_advection_local_rom_matrix_interp -p 1 -rp 1 -dt 0.005 -tf 9
//    mpirun -np 4 dg_advection_local_rom_matrix_interp -p 3 -rp 2 -dt 0.0025 -tf 9 -vs 20
//    mpirun -np 4 dg_advection_local_rom_matrix_interp -p 0 -o 2 -rp 1 -dt 0.01 -tf 8
//    mpirun -np 4 dg_advection_local_rom_matrix_interp -p 0 -rs 2 -dt 0.005 -tf 2
//    mpirun -np 4 dg_advection_local_rom_matrix_interp -p 0 -rs 1 -o 2 -tf 2
//    mpirun -np 3 dg_advection_local_rom_matrix_interp -p 1 -rs 1 -rp 0 -dt 0.005 -tf 0.5
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

#include "mfem.hpp"
#include "linalg/Vector.h"
#include "algo/manifold_interp/MatrixInterpolator.h"
#include "algo/manifold_interp/VectorInterpolator.h"
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
    const char *rbf_type = "G";
    const char *interp_method = "LS";
    double closest_rbf_val = 0.9;
    f_factor = 1.0;
    int rdim = -1;
    bool fom = false;
    bool offline = false;
    bool online = false;
    bool online_interp = false;
    bool interp_prep = false;
    bool visualization = true;
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
    args.AddOption(&online_interp, "-online_interp", "--online_interp",
                   "-no-online_interp", "--no-online_interp",
                   "Enable or disable matrix interpolation during the online phase.");
    args.AddOption(&interp_prep, "-interp_prep", "--interp_prep", "-no-interp_prep",
                   "--no-interp_prep",
                   "Enable or disable matrix interpolation preparation during the online phase.");
    args.AddOption(&rbf_type, "-rt", "--rbf_type",
                   "RBF type ('G' == gaussian, 'IQ' == inverse quadratic, 'IMQ' == inverse multiquadric).");
    args.AddOption(&interp_method, "-im", "--interp_method",
                   "Interpolation method ('LS' == linear solve, 'IDW'== inverse distance weighting, 'LP' == lagrangian polynomials).");
    args.AddOption(&closest_rbf_val, "-crv", "--crv",
                   "RBF value of the two closest points.");
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

    Device device(device_config);
    if (mpi.Root()) {
        device.Print();
    }

    bool check = (!online && !interp_prep && !online_interp) ||
                 (online && !interp_prep && !online_interp) || (!online && interp_prep
                         && !online_interp) ||
                 (!online && !interp_prep && online_interp);
    MFEM_VERIFY(check,
                "only one of online, interp_prep, or online_interp can be true!");

    if (interp_prep || online_interp)
    {
        online = true;
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

    {
        ostringstream mesh_name, sol_name;
        mesh_name << "dg_advection_local_rom_matrix_interp-mesh." << setfill('0') <<
                  setw(6) << myid;
        sol_name << "dg_advection_local_rom_matrix_interp-init." << setfill('0') <<
                 setw(6) << myid;
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
        if (binary)
        {
#ifdef MFEM_USE_SIDRE
            dc = new SidreDataCollection("DG_Advection", pmesh);
#else
            MFEM_ABORT("Must build with MFEM_USE_SIDRE=YES for binary output.");
#endif
        }
        else
        {
            dc = new VisItDataCollection("DG_Advection", pmesh);
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
        pd = new ParaViewDataCollection("DG_Advection", pmesh);
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
        const std::string collection_name = "dg_advection_local_rom_matrix_interp-p-" +
                                            postfix + ".bp";

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

    StopWatch fom_timer;
    double t = 0.0;

    int max_num_snapshots = t_final / dt + 1;
    bool update_right_SV = false;
    bool isIncremental = false;
    const std::string basisName = "basis_" + std::to_string(f_factor);
    const CAROM::Matrix* spatialbasis;
    CAROM::Options* options;
    CAROM::BasisGenerator *generator;
    int numRowRB, numColumnRB;

    CAROM::Matrix *M_hat_carom, *K_hat_carom;
    DenseMatrix *M_hat, *K_hat;
    CAROM::Vector *b_hat_carom, *u_init_hat_carom;
    Vector *b_hat, *u_init_hat;

    Vector u_init(*U);
    Vector *u_in;

    // 10. Set BasisGenerator if offline
    if (offline)
    {
        options = new CAROM::Options(U->Size(), max_num_snapshots, 1, update_right_SV);
        generator = new CAROM::BasisGenerator(*options, isIncremental, basisName);
        Vector u_curr(*U);
        Vector u_centered(U->Size());
        subtract(u_curr, u_init, u_centered);
        bool addSample = generator->takeSample(u_centered.GetData(), t, dt);
    }

    if (online)
    {
        if (!online_interp)
        {
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
            if (interp_prep) M_hat_carom->write("M_hat_" + std::to_string(f_factor));

            // libROM stores the matrix row-wise, so wrapping as a DenseMatrix in MFEM means it is transposed.
            M_hat = new DenseMatrix(numColumnRB, numColumnRB);
            M_hat->Set(1, M_hat_carom->getData());
            M_hat->Transpose();

            K_hat_carom = new CAROM::Matrix(numRowRB, numColumnRB, false);
            ComputeCtAB(K_mat, *spatialbasis, *spatialbasis, *K_hat_carom);
            if (interp_prep) K_hat_carom->write("K_hat_" + std::to_string(f_factor));

            // libROM stores the matrix row-wise, so wrapping as a DenseMatrix in MFEM means it is transposed.
            K_hat = new DenseMatrix(numColumnRB, numColumnRB);
            K_hat->Set(1, K_hat_carom->getData());
            K_hat->Transpose();

            Vector b_vec = *B;
            CAROM::Vector b_carom(b_vec.GetData(), b_vec.Size(), true);
            b_hat_carom = spatialbasis->transposeMult(&b_carom);
            if (interp_prep) b_hat_carom->write("b_hat_" + std::to_string(f_factor));
            b_hat = new Vector(b_hat_carom->getData(), b_hat_carom->dim());

            u_init_hat_carom = new CAROM::Vector(numColumnRB, false);
            ComputeCtAB_vec(K_mat, *U, *spatialbasis, *u_init_hat_carom);
            if (interp_prep) u_init_hat_carom->write("u_init_hat_" + std::to_string(
                            f_factor));
            u_init_hat = new Vector(u_init_hat_carom->getData(), u_init_hat_carom->dim());

            if (interp_prep)
            {
                if (myid == 0)
                {
                    std::ofstream fout;
                    fout.open("frequencies.txt", std::ios::app);
                    fout << f_factor << std::endl;
                    fout.close();
                }
                return 0;
            }
        }
        else
        {
            std::fstream fin("frequencies.txt", std::ios_base::in);
            double freq;
            std::set<double> frequencies;
            while (fin >> freq)
            {
                frequencies.insert(freq);
            }
            fin.close();

            std::vector<CAROM::Vector*> parameter_points;
            std::vector<CAROM::Matrix*> bases;
            std::vector<CAROM::Matrix*> K_hats;
            std::vector<CAROM::Matrix*> M_hats;
            std::vector<CAROM::Vector*> b_hats;
            std::vector<CAROM::Vector*> u_init_hats;
            std::ofstream fout;
            fout.open("frequencies.txt");
            for(auto it = frequencies.begin(); it != frequencies.end(); it++)
            {
                CAROM::Vector* point = new CAROM::Vector(1, false);
                point->item(0) = *it;
                parameter_points.push_back(point);
                fout << *it << std::endl;

                std::string parametricBasisName = "basis_" + std::to_string(*it);
                CAROM::BasisReader reader(parametricBasisName);

                MFEM_VERIFY(rdim != -1, "rdim must be used for interpolation.");
                CAROM::Matrix* parametricSpatialBasis = reader.getSpatialBasis(0.0, rdim);
                numRowRB = parametricSpatialBasis->numRows();
                numColumnRB = parametricSpatialBasis->numColumns();
                bases.push_back(parametricSpatialBasis);

                CAROM::Matrix* parametricMhat = new CAROM::Matrix();
                parametricMhat->read("M_hat_" + std::to_string(*it));
                M_hats.push_back(parametricMhat);

                CAROM::Matrix* parametricKhat = new CAROM::Matrix();
                parametricKhat->read("K_hat_" + std::to_string(*it));
                K_hats.push_back(parametricKhat);

                CAROM::Vector* parametricbhat = new CAROM::Vector();
                parametricbhat->read("b_hat_" + std::to_string(*it));
                b_hats.push_back(parametricbhat);

                CAROM::Vector* parametricuinithat = new CAROM::Vector();
                parametricuinithat->read("u_init_hat_" + std::to_string(*it));
                u_init_hats.push_back(parametricuinithat);
            }
            fout.close();
            if (myid == 0) printf("spatial basis dimension is %d x %d\n", numRowRB,
                                      numColumnRB);

            CAROM::Vector* curr_point = new CAROM::Vector(1, false);
            curr_point->item(0) = f_factor;

            int ref_point = getClosestPoint(parameter_points, curr_point);
            std::vector<CAROM::Matrix*> rotation_matrices = obtainRotationMatrices(
                        parameter_points, bases, ref_point);

            CAROM::MatrixInterpolator basis_interpolator(parameter_points,
                    rotation_matrices, bases, ref_point, "B", rbf_type, interp_method,
                    closest_rbf_val);
            CAROM::MatrixInterpolator M_interpolator(parameter_points, rotation_matrices,
                    M_hats, ref_point, "SPD", rbf_type, interp_method, closest_rbf_val);
            CAROM::MatrixInterpolator K_interpolator(parameter_points, rotation_matrices,
                    K_hats, ref_point, "R", rbf_type, interp_method, closest_rbf_val);
            CAROM::VectorInterpolator b_interpolator(parameter_points, rotation_matrices,
                    b_hats, ref_point, rbf_type, interp_method, closest_rbf_val);
            CAROM::VectorInterpolator u_init_interpolator(parameter_points,
                    rotation_matrices, u_init_hats, ref_point, rbf_type, interp_method,
                    closest_rbf_val);
            spatialbasis = basis_interpolator.interpolate(curr_point);
            M_hat_carom = M_interpolator.interpolate(curr_point);
            K_hat_carom = K_interpolator.interpolate(curr_point);
            b_hat_carom = b_interpolator.interpolate(curr_point);
            u_init_hat_carom = u_init_interpolator.interpolate(curr_point);

            // libROM stores the matrix row-wise, so wrapping as a DenseMatrix in MFEM means it is transposed.
            M_hat = new DenseMatrix(numColumnRB, numColumnRB);
            M_hat->Set(1, M_hat_carom->getData());
            M_hat->Transpose();

            // libROM stores the matrix row-wise, so wrapping as a DenseMatrix in MFEM means it is transposed.
            K_hat = new DenseMatrix(numColumnRB, numColumnRB);
            K_hat->Set(1, K_hat_carom->getData());
            K_hat->Transpose();

            b_hat = new Vector(b_hat_carom->getData(), b_hat_carom->dim());
            u_init_hat = new Vector(u_init_hat_carom->getData(), u_init_hat_carom->dim());
        }

        u_in = new Vector(numColumnRB);
        *u_in = 0.0;
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

        fom_timer.Start();

        double dt_real = min(dt, t_final - t);
        if (online)
        {
            ode_solver->Step(*u_in, t, dt_real);
        }
        else
        {
            ode_solver->Step(*U, t, dt_real);
        }
        ti++;
        done = (t >= t_final - 1e-8*dt);

        fom_timer.Stop();

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
        generator->endSamples();
        delete generator;
        delete options;
    }

    if (online)
    {
        CAROM::Vector u_hat_final_carom(u_in->GetData(), u_in->Size(), false);
        CAROM::Vector* u_final_carom = spatialbasis->mult(u_hat_final_carom);
        Vector u_final(u_final_carom->getData(), u_final_carom->dim());
        u_final += u_init;

        Vector fom_solution(U->Size());
        ifstream solution_file;
        ostringstream solution_filename;
        solution_filename << "dg_advection_local_rom_matrix_interp-final." << f_factor
                          << "." << setfill('0') << setw(6) << myid;
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
        delete u_in;
        delete u_final_carom;
    }

    // 12. Save the final solution in parallel. This output can be viewed later
    //     using GLVis: "glvis -np <np> -m dg_advection_local_rom_matrix_interp-mesh -g dg_advection_local_rom_matrix_interp-final".
    if (offline || fom)
    {
        *u = *U;
        ostringstream sol_name;
        sol_name << "dg_advection_local_rom_matrix_interp-final." << f_factor << "." <<
                 setfill('0') << setw(6) << myid;
        ofstream osol(sol_name.str().c_str());
        osol.precision(precision);
        Vector tv(u->ParFESpace()->GetTrueVSize());
        u->GetTrueDofs(tv);

        for (int i=0; i<tv.Size(); ++i)
            osol << tv[i] << std::endl;

        osol.close();
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
