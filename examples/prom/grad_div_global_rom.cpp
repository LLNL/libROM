/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

//               libROM MFEM Example: parametric ROM for Grad-div problem (adapted from MFEM ex4p)
//
// Description:  This example code demonstrates the use of MFEM and libROM to
//               define a simple projection-based reduced order model of a
//               simple 2D/3D H(div) diffusion problem corresponding to the second
//               order definite equation -grad(alpha div F) + beta F = f with
//               boundary condition F dot n = <given normal field>.
//               Here, we use a given exact solution F and compute the corresponding
//               r.h.s. f. In this example, the ROM parameter is a scalar frequency
//               used in defining the exact solution F and corresponding source function f.
//               We discretize with Raviart-Thomas finite elements.
//
//               The example highlights three distinct ROM processes, i.e.,
//               offline, merge, and online. The offline phase runs the full
//               order model and stores the snapshot data in an HDF file. You
//               can run as many offline phases as you wish to sample the
//               parameter space. The merge phase reads all the snapshot files,
//               builds a global reduced basis, and stores the basis in an HDF
//               file. The online phase reads the basis, builds the ROM
//               operator, solves the reduced order system, and lifts the
//               solution to the full order space.
//
// Offline phase: grad_div_global_rom -offline -f 1.0 -id 0
//                grad_div_global_rom -offline -f 1.1 -id 1
//                grad_div_global_rom -offline -f 1.2 -id 2
//
// Merge phase:   grad_div_global_rom -merge -ns 3
//
// FOM run (for error calculation):
//                grad_div_global_rom -fom -f 1.15
//
// Online phase:  grad_div_global_rom -online -f 1.15
//
// This example runs in parallel with MPI, by using the same number of MPI ranks
// in all phases (offline, merge, fom, online).

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include "linalg/BasisGenerator.h"
#include "linalg/BasisReader.h"
#include "mfem/Utilities.hpp"

using namespace std;
using namespace mfem;

// Exact solution, F, and r.h.s., f. See below for implementation.
void F_exact(const Vector &, Vector &);
void f_exact(const Vector &, Vector &);
double freq = 1.0, kappa;

int main(int argc, char *argv[])
{
    // 1. Initialize MPI and HYPRE.
    int num_procs, myid;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    Hypre::Init();

    // 2. Parse command-line options.
    const char *mesh_file = "../data/star.mesh";
    int order = 1;
    bool set_bc = true;
    bool static_cond = false;
    bool hybridization = false;
    bool pa = false;
    const char *device_config = "cpu";
    bool visualization = true;
    bool visit = false;
    bool fom = false;
    bool offline = false;
    bool merge = false;
    bool online = false;
    int precision = 8;
    int id = 0;
    int nsets = 0;
    double coef = 1.0;

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&order, "-o", "--order",
                   "Finite element order (polynomial degree).");
    args.AddOption(&set_bc, "-bc", "--impose-bc", "-no-bc", "--dont-impose-bc",
                   "Impose or not essential boundary conditions.");
    args.AddOption(&freq, "-f", "--frequency", "Set the frequency for the exact"
                   " solution.");
    args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                   "--no-static-condensation", "Enable static condensation.");
    args.AddOption(&hybridization, "-hb", "--hybridization", "-no-hb",
                   "--no-hybridization", "Enable hybridization.");
    args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                   "--no-partial-assembly", "Enable Partial Assembly.");
    args.AddOption(&id, "-id", "--id", "Parametric id");
    args.AddOption(&nsets, "-ns", "--nset", "Number of parametric snapshot sets");
    args.AddOption(&coef, "-cf", "--coefficient",
                   "Coefficient.");
    args.AddOption(&device_config, "-d", "--device",
                   "Device configuration string, see Device::Configure().");
    args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                   "--no-visit-datafiles",
                   "Save data files for VisIt (visit.llnl.gov) visualization.");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                   "--no-visualization",
                   "Enable or disable GLVis visualization.");
    args.AddOption(&fom, "-fom", "--fom", "-no-fom", "--no-fom",
                   "Enable or disable the fom phase.");
    args.AddOption(&offline, "-offline", "--offline", "-no-offline", "--no-offline",
                   "Enable or disable the offline phase.");
    args.AddOption(&online, "-online", "--online", "-no-online", "--no-online",
                   "Enable or disable the online phase.");
    args.AddOption(&merge, "-merge", "--merge", "-no-merge", "--no-merge",
                   "Enable or disable the merge phase.");
    args.Parse();
    if (!args.Good())
    {
        if (myid == 0)
        {
            args.PrintUsage(cout);
        }
        MPI_Finalize();
        return 1;
    }
    if (myid == 0)
    {
        args.PrintOptions(cout);
    }
    kappa = freq * M_PI;

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

    // 3. Enable hardware devices such as GPUs, and programming models such as
    //    CUDA, OCCA, RAJA and OpenMP based on command line options.
    Device device(device_config);
    if (myid == 0) {
        device.Print();
    }

    // 4. Read the (serial) mesh from the given mesh file on all processors.  We
    //    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
    //    and volume meshes with the same code.
    Mesh mesh(mesh_file, 1, 1);
    int dim = mesh.Dimension();
    int sdim = mesh.SpaceDimension();

    // 5. Refine the serial mesh on all processors to increase the resolution. In
    //    this example we do 'ref_levels' of uniform refinement. We choose
    //    'ref_levels' to be the largest number that gives a final mesh with no
    //    more than 1,000 elements.
    {
        int ref_levels =
            (int)floor(log(1000./mesh.GetNE())/log(2.)/dim);
        for (int l = 0; l < ref_levels; l++)
        {
            mesh.UniformRefinement();
        }
    }

    // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
    //    this mesh further in parallel to increase the resolution. Once the
    //    parallel mesh is defined, the serial mesh can be deleted.
    ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();
    {
        int par_ref_levels = 2;
        for (int l = 0; l < par_ref_levels; l++)
        {
            pmesh.UniformRefinement();
        }
    }

    // 7. Define a parallel finite element space on the parallel mesh. Here we
    //    use the Raviart-Thomas finite elements of the specified order.
    FiniteElementCollection *fec = new RT_FECollection(order-1, dim);
    ParFiniteElementSpace fespace(&pmesh, fec);
    HYPRE_BigInt size = fespace.GlobalTrueVSize();
    if (myid == 0)
    {
        cout << "Number of finite element unknowns: " << size << endl;
    }

    // 8. Determine the list of true (i.e. parallel conforming) essential
    //    boundary dofs. In this example, the boundary conditions are defined
    //    by marking all the boundary attributes from the mesh as essential
    //    (Dirichlet) and converting them to a list of true dofs.
    Array<int> ess_tdof_list;
    if (pmesh.bdr_attributes.Size())
    {
        Array<int> ess_bdr(pmesh.bdr_attributes.Max());
        ess_bdr = set_bc ? 1 : 0;
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    }

    // 9. Initiate ROM related variables
    int max_num_snapshots = 100;
    bool update_right_SV = false;
    bool isIncremental = false;
    const std::string basisName = "basis";
    const std::string basisFileName = basisName + std::to_string(id);
    CAROM::Options* options;
    CAROM::BasisGenerator *generator;
    StopWatch solveTimer, assembleTimer, mergeTimer;

    // 10. Set BasisGenerator if offline
    if (offline)
    {
        options = new CAROM::Options(fespace.GetTrueVSize(), max_num_snapshots,
                                     update_right_SV);
        generator = new CAROM::BasisGenerator(*options, isIncremental, basisFileName);
    }

    // 11. The merge phase
    if (merge)
    {
        mergeTimer.Start();
        options = new CAROM::Options(fespace.GetTrueVSize(), max_num_snapshots,
                                     update_right_SV);
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

    // 12. Set up the parallel linear form b(.) which corresponds to the
    //     right-hand side of the FEM linear system, which in this case is
    //     (f,phi_i) where f is given by the function f_exact and phi_i are the
    //     basis functions in the finite element fespace.
    assembleTimer.Start();
    VectorFunctionCoefficient f(sdim, f_exact);
    ParLinearForm *b = new ParLinearForm(&fespace);
    b->AddDomainIntegrator(new VectorFEDomainLFIntegrator(f));
    b->Assemble();

    // 13. Define the solution vector x as a parallel finite element grid function
    //    corresponding to fespace. Initialize x by projecting the exact
    //    solution. Note that only values from the boundary faces will be used
    //    when eliminating the non-homogeneous boundary condition to modify the
    //    r.h.s. vector b.
    ParGridFunction x(&fespace);
    VectorFunctionCoefficient F(sdim, F_exact);
    x.ProjectCoefficient(F);

    // 14. Set up the parallel bilinear form corresponding to the H(div)
    //     diffusion operator grad alpha div + beta I, by adding the div-div and
    //     the mass domain integrators.
    ConstantCoefficient alpha(1.0);
    ConstantCoefficient beta(1.0);
    ParBilinearForm a(&fespace);
    if (pa) {
        a.SetAssemblyLevel(AssemblyLevel::PARTIAL);
    }
    a.AddDomainIntegrator(new DivDivIntegrator(alpha));
    a.AddDomainIntegrator(new VectorFEMassIntegrator(beta));

    // 15. Assemble the parallel bilinear form and the corresponding linear
    //     system, applying any necessary transformations such as: parallel
    //     assembly, eliminating boundary conditions, applying conforming
    //     constraints for non-conforming AMR, static condensation,
    //     hybridization, etc.
    FiniteElementCollection *hfec = NULL;
    ParFiniteElementSpace *hfes = NULL;
    if (static_cond)
    {
        a.EnableStaticCondensation();
    }
    else if (hybridization)
    {
        hfec = new DG_Interface_FECollection(order-1, dim);
        hfes = new ParFiniteElementSpace(&pmesh, hfec);
        a.EnableHybridization(hfes, new NormalTraceJumpIntegrator(),
                              ess_tdof_list);
    }
    a.Assemble();

    OperatorPtr A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
    assembleTimer.Stop();

    if (myid == 0 && !pa)
    {
        cout << "Size of linear system: "
             << A.As<HypreParMatrix>()->GetGlobalNumRows() << endl;
    }

    // 16. The offline phase
    if(fom || offline)
    {
        // 17. Define and apply a parallel PCG solver for A X = B with the 2D AMS or
        //     the 3D ADS preconditioners from hypre. If using hybridization, the
        //     system is preconditioned with hypre's BoomerAMG. In the partial
        //     assembly case, use Jacobi preconditioning.
        Solver *prec = NULL;
        CGSolver *pcg = new CGSolver(MPI_COMM_WORLD);
        pcg->SetOperator(*A);
        pcg->SetRelTol(1e-12);
        pcg->SetMaxIter(2000);
        pcg->SetPrintLevel(1);
        if (hybridization) {
            prec = new HypreBoomerAMG(*A.As<HypreParMatrix>());
        }
        else if (pa) {
            prec = new OperatorJacobiSmoother(a, ess_tdof_list);
        }
        else
        {
            ParFiniteElementSpace *prec_fespace =
                (a.StaticCondensationIsEnabled() ? a.SCParFESpace() : &fespace);
            if (dim == 2) {
                prec = new HypreAMS(*A.As<HypreParMatrix>(), prec_fespace);
            }
            else {
                prec = new HypreADS(*A.As<HypreParMatrix>(), prec_fespace);
            }
        }
        pcg->SetPreconditioner(*prec);
        solveTimer.Start();
        pcg->Mult(B, X);
        solveTimer.Stop();
        delete pcg;
        delete prec;

        // 18. take and write snapshot for ROM
        if (offline)
        {
            generator->takeSample(X.GetData());
            generator->writeSnapshot();
            delete generator;
            delete options;
        }
    }

    // 19. The online phase
    if (online) {
        // 20. read the reduced basis
        assembleTimer.Start();
        CAROM::BasisReader reader(basisName);
        std::unique_ptr<const CAROM::Matrix> spatialbasis = reader.getSpatialBasis();
        const int numRowRB = spatialbasis->numRows();
        const int numColumnRB = spatialbasis->numColumns();
        if (myid == 0) printf("spatial basis dimension is %d x %d\n", numRowRB,
                                  numColumnRB);

        // 21. form inverse ROM operator
        CAROM::Matrix invReducedA(numColumnRB, numColumnRB, false);
        ComputeCtAB(*A.As<HypreParMatrix>(), *spatialbasis, *spatialbasis, invReducedA);
        invReducedA.inverse();

        CAROM::Vector B_carom(B.GetData(), B.Size(), true, false);
        CAROM::Vector X_carom(X.GetData(), X.Size(), true, false);
        std::unique_ptr<CAROM::Vector> reducedRHS = spatialbasis->transposeMult(
                    B_carom);
        CAROM::Vector reducedSol(numColumnRB, false);
        assembleTimer.Stop();

        // 22. solve ROM
        solveTimer.Start();
        invReducedA.mult(*reducedRHS, reducedSol);
        solveTimer.Stop();

        // 23. reconstruct FOM state
        spatialbasis->mult(reducedSol, X_carom);
    }

    // 24. Recover the parallel grid function corresponding to X. This is the
    //     local finite element solution on each processor.
    a.RecoverFEMSolution(X, *b, x);

    // 25. Compute and print the L^2 norm of the error.
    {
        double error = x.ComputeL2Error(F);
        if (myid == 0)
        {
            cout << "\nAbsolute error in L^2 norm with respect to the exact solution:\n || F_h - F ||_{L^2} = "
                 << error << '\n' << endl;
        }
    }

    // 26. Calculate the relative error of the ROM prediction compared to FOM
    ostringstream sol_dofs_name, sol_dofs_name_fom;
    if (fom || offline)
    {
        sol_dofs_name << "sol_dofs_fom." << setfill('0') << setw(6) << myid;
    }
    if (online)
    {
        sol_dofs_name << "sol_dofs." << setfill('0') << setw(6) << myid;
        sol_dofs_name_fom << "sol_dofs_fom." << setfill('0') << setw(6) << myid;
    }

    if (online)
    {
        // Initialize FOM solution
        Vector x_fom(x.Size());

        ifstream fom_file;

        // Open and load file
        fom_file.open(sol_dofs_name_fom.str().c_str());

        x_fom.Load(fom_file, x_fom.Size());

        fom_file.close();

        Vector diff_x(x.Size());

        subtract(x, x_fom, diff_x);

        // Get norms
        double tot_diff_norm = sqrt(InnerProduct(MPI_COMM_WORLD, diff_x, diff_x));

        double tot_fom_norm = sqrt(InnerProduct(MPI_COMM_WORLD,
                                                x_fom, x_fom));

        if (myid == 0)
        {
            cout << "Relative error of ROM solution = "
                 << tot_diff_norm / tot_fom_norm << endl;
        }
    }

    // 27. Save the refined mesh and the solution in parallel. This output can
    //     be viewed later using GLVis: "glvis -np <np> -m mesh -g sol".
    {
        ostringstream mesh_name, sol_name;
        mesh_name << "mesh." << setfill('0') << setw(6) << myid;
        sol_name << "sol." << setfill('0') << setw(6) << myid;

        ofstream mesh_ofs(mesh_name.str().c_str());
        mesh_ofs.precision(precision);
        pmesh.Print(mesh_ofs);

        ofstream sol_ofs(sol_name.str().c_str());
        sol_ofs.precision(precision);
        x.Save(sol_ofs);

        ofstream sol_dofs_ofs(sol_dofs_name.str().c_str());
        sol_dofs_ofs.precision(16);
        for (int i = 0; i < x.Size(); ++i)
        {
            sol_dofs_ofs << x[i] << std::endl;
        }
    }

    // 28. Save data in the VisIt format.
    DataCollection *dc = NULL;
    if (visit)
    {
        if (offline) dc = new VisItDataCollection("Example4", &pmesh);
        else if (fom) dc = new VisItDataCollection("Example4_fom", &pmesh);
        else if (online) dc = new VisItDataCollection("Example4_rom", &pmesh);
        dc->SetPrecision(precision);
        dc->RegisterField("solution", &x);
        dc->Save();
        delete dc;
    }

    // 29. Send the solution by socket to a GLVis server.
    if (visualization)
    {
        char vishost[] = "localhost";
        int  visport   = 19916;
        socketstream sol_sock(vishost, visport);
        sol_sock << "parallel " << num_procs << " " << myid << "\n";
        sol_sock.precision(precision);
        sol_sock << "solution\n" << pmesh << x << flush;
    }

    // 30. print timing info
    if (myid == 0)
    {
        if(fom || offline)
        {
            printf("Elapsed time for assembling FOM: %e second\n",
                   assembleTimer.RealTime());
            printf("Elapsed time for solving FOM: %e second\n", solveTimer.RealTime());
        }
        if(online)
        {
            printf("Elapsed time for assembling ROM: %e second\n",
                   assembleTimer.RealTime());
            printf("Elapsed time for solving ROM: %e second\n", solveTimer.RealTime());
        }
    }

    // 31. Free the used memory.
    delete b;
    delete fec;
    if (hybridization)
    {
        delete hfes;
        delete hfec;
    }
    MPI_Finalize();

    return 0;
}

// The exact solution (for non-surface meshes)
void F_exact(const Vector &p, Vector &F)
{
    int dim = p.Size();

    double x = p(0);
    double y = p(1);
    // double z = (dim == 3) ? p(2) : 0.0; // Uncomment if F is changed to depend on z

    F(0) = cos(kappa*x)*sin(kappa*y);
    F(1) = cos(kappa*y)*sin(kappa*x);
    if (dim == 3)
    {
        F(2) = 0.0;
    }
}

// The right hand side
void f_exact(const Vector &p, Vector &f)
{
    int dim = p.Size();

    double x = p(0);
    double y = p(1);
    // double z = (dim == 3) ? p(2) : 0.0; // Uncomment if f is changed to depend on z

    double temp = 1 + 2*kappa*kappa;

    f(0) = temp*cos(kappa*x)*sin(kappa*y);
    f(1) = temp*cos(kappa*y)*sin(kappa*x);
    if (dim == 3)
    {
        f(2) = 0;
    }
}
