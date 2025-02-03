/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

//               libROM MFEM Example: parametric ROM for Maxwell equation (adapted from ex3p.cpp)
//
// Compile with: make maxwell_global_rom
//
// Description:  This example code uses MFEM and librom to define a projection-based
//               reduced-order model for a simple electromagnetic diffusion
//               problem corresponding to the second order definite Maxwell
//               equation: curl curl E + E = f. The boundary condition is given by
//               E x n = <given tangential field>. Here, we use a given exact
//               solution E and compute the corresponding r.h.s. f. Each component
//               of the vector field E is assumed to be a sin wave, with frequency
//               controlled by the parameter kappa.
//
//               The example highlights three distinct ROM processes - the
//               offline, merge, and online phases. The offline phase runs the full
//               order model and stores the snapshot data in an HDF file. You
//               can run as many offline phases as you wish to sample the
//               parameter space. The merge phase reads all the snapshot files,
//               builds a global reduced basis, and stores the basis in an HDF
//               file. The online phase reads the basis, builds the ROM
//               operator, solves the reduced order system, and lifts the
//               solution to the full order space.
//
// For ROM (global rom):
// Offline phase: maxwell_global_rom -offline -f 1.0 -id 0
//                maxwell_global_rom -offline -f 1.1 -id 1
//                maxwell_global_rom -offline -f 1.2 -id 2
//
// Merge phase:   maxwell_global_rom -merge -ns 3
//
// FOM solution:  maxwell_global_rom -fom -f 1.15
//
// Online phase:  maxwell_global_rom -online -f 1.15
//
// This example runs in parallel with MPI, by using the same number of MPI ranks
// in all phases (offline, merge, fom, online).

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include "linalg/Vector.h"
#include "linalg/BasisGenerator.h"
#include "linalg/BasisReader.h"
#include "mfem/Utilities.hpp"
#include <cmath>
#include <set>

using namespace std;
using namespace mfem;

// Exact solution, E, and r.h.s., f. See below for implementation.
void E_exact(const Vector &, Vector &);
void f_exact(const Vector &, Vector &);
double freq = 1.0, kappa;
int dim;

int main(int argc, char *argv[])
{
    // 1. Initialize MPI and HYPRE
    Mpi::Init(argc, argv);
    int num_procs = Mpi::WorldSize();
    int myid = Mpi::WorldRank();
    Hypre::Init();

    // 2. Parse command-line options.
    const char *mesh_file = "../data/star.mesh";
    int order = 1;
    bool static_cond = false;
    bool pa = false;
    const char *device_config = "cpu";
    bool visualization = true;
    bool visit = true;
    bool fom = false;
    bool offline = false;
    bool merge = false;
    bool online = false;
    int precision = 8;
    int id = 0;
    int nsets = 0;
    double coef = 1.0;
#ifdef MFEM_USE_AMGX
    bool useAmgX = false;
#endif


    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&order, "-o", "--order",
                   "Finite element order (polynomial degree).");
    args.AddOption(&freq, "-f", "--frequency", "Set the frequency for the exact"
                   " solution.");
    args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                   "--no-static-condensation", "Enable static condensation.");
    args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                   "--no-partial-assembly", "Enable Partial Assembly.");
    args.AddOption(&device_config, "-d", "--device",
                   "Device configuration string, see Device::Configure().");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                   "--no-visualization",
                   "Enable or disable GLVis visualization.");
#ifdef MFEM_USE_AMGX
    args.AddOption(&useAmgX, "-amgx", "--useAmgX", "-no-amgx",
                   "--no-useAmgX",
                   "Enable or disable AmgX in MatrixFreeAMS.");
#endif
    args.AddOption(&fom, "-fom", "--fom", "-no-fom", "--no-fom",
                   "Enable or disable the fom phase.");
    args.AddOption(&offline, "-offline", "--offline", "-no-offline","--no-offline",
                   "Enable or disable the offline phase.");
    args.AddOption(&online, "-online", "--online", "-no-online", "--no-online",
                   "Enable or disable the online phase.");
    args.AddOption(&merge, "-merge", "--merge", "-no-merge", "--no-merge",
                   "Enable or disable the merge phase.");
    args.AddOption(&id, "-id", "--id", "Parametric id");
    args.AddOption(&nsets, "-ns", "--nset", "Number of parametric snapshot sets");

    args.Parse();
    if (!args.Good())
    {
        if (myid == 0)
        {
            args.PrintUsage(cout);
        }
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
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    dim = mesh->Dimension();
    int sdim = mesh->SpaceDimension();

    // 5. Refine the serial mesh on all processors to increase the resolution. In
    //    this example we do 'ref_levels' of uniform refinement. We choose
    //    'ref_levels' to be the largest number that gives a final mesh with no
    //    more than 1,000 elements.
    {
        int ref_levels = (int)floor(log(1000./mesh->GetNE())/log(2.)/dim);
        for (int l = 0; l < ref_levels; l++)
        {
            mesh->UniformRefinement();
        }
    }

    // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
    //    this mesh further in parallel to increase the resolution. Once the
    //    parallel mesh is defined, the serial mesh can be deleted.
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
    {
        int par_ref_levels = 2;
        for (int l = 0; l < par_ref_levels; l++)
        {
            pmesh->UniformRefinement();
        }
    }

    // 7. Define a parallel finite element space on the parallel mesh. Here we
    //    use the Nedelec finite elements of the specified order.
    FiniteElementCollection *fec = new ND_FECollection(order, dim);
    ParFiniteElementSpace *fespace = new ParFiniteElementSpace(pmesh, fec);
    HYPRE_BigInt size = fespace->GlobalTrueVSize();
    if (myid == 0)
    {
        cout << "Number of finite element unknowns: " << size << endl;
    }

    // 8. Determine the list of true (i.e. parallel conforming) essential
    //    boundary dofs. In this example, the boundary conditions are defined
    //    by marking all the boundary attributes from the mesh as essential
    //    (Dirichlet) and converting them to a list of true dofs.
    Array<int> ess_tdof_list;
    Array<int> ess_bdr;
    if (pmesh->bdr_attributes.Size())
    {
        ess_bdr.SetSize(pmesh->bdr_attributes.Max());
        ess_bdr = 1;
        fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    }

    // 9. Initiate ROM related variables
    int max_num_snapshots = 100;
    bool update_right_SV = false;
    bool isIncremental = false;
    const std::string basisName = "basis";
    const std::string basisFileName = basisName + std::to_string(id);
    std::unique_ptr<const CAROM::Matrix> spatialbasis;
    CAROM::Options* options;
    CAROM::BasisGenerator *generator;
    int numRowRB, numColumnRB;
    StopWatch solveTimer, assembleTimer, mergeTimer;

    // 10. Set BasisGenerator if offline
    if (offline)
    {
        options = new CAROM::Options(fespace->GetTrueVSize(), max_num_snapshots,
                                     update_right_SV);
        generator = new CAROM::BasisGenerator(*options, isIncremental, basisFileName);
    }

    // 11. The merge phase
    if (merge)
    {
        mergeTimer.Start();
        options = new CAROM::Options(fespace->GetTrueVSize(), max_num_snapshots,
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
    //    right-hand side of the FEM linear system, which in this case is
    //    (f,phi_i) where f is given by the function f_exact and phi_i are the
    //    basis functions in the finite element fespace.
    VectorFunctionCoefficient f(sdim, f_exact);
    ParLinearForm *b = new ParLinearForm(fespace);
    b->AddDomainIntegrator(new VectorFEDomainLFIntegrator(f));
    b->Assemble();

    // 13. Define the solution vector x as a parallel finite element grid function
    //     corresponding to fespace. Initialize x by projecting the exact
    //     solution. Note that only values from the boundary edges will be used
    //     when eliminating the non-homogeneous boundary condition to modify the
    //     r.h.s. vector b.
    ParGridFunction x(fespace);
    VectorFunctionCoefficient E(sdim, E_exact);
    x.ProjectCoefficient(E);

    // 14. Set up the parallel bilinear form corresponding to the EM diffusion
    //     operator curl muinv curl + sigma I, by adding the curl-curl and the
    //     mass domain integrators.
    Coefficient *muinv = new ConstantCoefficient(1.0);
    Coefficient *sigma = new ConstantCoefficient(1.0);
    ParBilinearForm *a = new ParBilinearForm(fespace);
    if (pa) {
        a->SetAssemblyLevel(AssemblyLevel::PARTIAL);
    }
    a->AddDomainIntegrator(new CurlCurlIntegrator(*muinv));
    a->AddDomainIntegrator(new VectorFEMassIntegrator(*sigma));

    // 15. Assemble the parallel bilinear form and the corresponding linear
    //     system, applying any necessary transformations such as: parallel
    //     assembly, eliminating boundary conditions, applying conforming
    //     constraints for non-conforming AMR, static condensation, etc.
    if (static_cond) {
        a->EnableStaticCondensation();
    }
    a->Assemble();

    OperatorPtr A;
    Vector B, X;
    a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
    assembleTimer.Stop();

    // 16. Offline Phase: Solve the system AX=B using PCG with an AMS preconditioner.
    if(fom || offline)
    {
        // 17. Solve the full order linear system A X = B
        if (pa)
        {
#ifdef MFEM_USE_AMGX
            MatrixFreeAMS ams(*a, *A, *fespace, muinv, sigma, NULL, ess_bdr, useAmgX);
#else
            MatrixFreeAMS ams(*a, *A, *fespace, muinv, sigma, NULL, ess_bdr);
#endif
            CGSolver cg(MPI_COMM_WORLD);
            cg.SetRelTol(1e-12);
            cg.SetMaxIter(1000);
            cg.SetPrintLevel(1);
            cg.SetOperator(*A);
            cg.SetPreconditioner(ams);
            solveTimer.Start();
            cg.Mult(B, X);
            solveTimer.Stop();
        }
        else
        {
            if (myid == 0)
            {
                cout << "Size of linear system: "
                     << A.As<HypreParMatrix>()->GetGlobalNumRows() << endl;
            }

            ParFiniteElementSpace *prec_fespace =
                (a->StaticCondensationIsEnabled() ? a->SCParFESpace() : fespace);
            HypreAMS ams(*A.As<HypreParMatrix>(), prec_fespace);
            HyprePCG pcg(*A.As<HypreParMatrix>());
            pcg.SetTol(1e-12);
            pcg.SetMaxIter(500);
            pcg.SetPrintLevel(2);
            pcg.SetPreconditioner(ams);
            solveTimer.Start();
            pcg.Mult(B, X);
            solveTimer.Stop();
        }
        // 18. take and write snapshot for ROM
        if (offline)
        {
            bool addSample = generator->takeSample(X.GetData());
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
        spatialbasis = reader.getSpatialBasis();
        numRowRB = spatialbasis->numRows();
        numColumnRB = spatialbasis->numColumns();
        if (myid == 0) printf("spatial basis dimension is %d x %d\n", numRowRB,
                                  numColumnRB);

        // 21. form inverse ROM operator
        CAROM::Matrix invReducedA(numColumnRB, numColumnRB, false);
        ComputeCtAB( *A.As<HypreParMatrix>(), *spatialbasis, *spatialbasis,
                     invReducedA);
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
    a->RecoverFEMSolution(X, *b, x);

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
        mesh_ofs.precision(8);
        pmesh->Print(mesh_ofs);

        ofstream sol_ofs(sol_name.str().c_str());
        sol_ofs.precision(8);
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
        if (offline) dc = new VisItDataCollection("maxwell", pmesh);
        else if (fom) dc = new VisItDataCollection("maxwell_fom", pmesh);
        else if (online) dc = new VisItDataCollection("maxwell_rom", pmesh);
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
        sol_sock.precision(8);
        sol_sock << "solution\n" << *pmesh << x << flush;
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
    delete a;
    delete sigma;
    delete muinv;
    delete b;
    delete fespace;
    delete fec;
    delete pmesh;

    return 0;
}

// 32. define Exact E
void E_exact(const Vector &x, Vector &E)
{
    if (dim == 3)
    {
        E(0) = sin(kappa * x(1));
        E(1) = sin(kappa * x(2));
        E(2) = sin(kappa * x(0));
    }
    else
    {
        E(0) = sin(kappa * x(1));
        E(1) = sin(kappa * x(0));
        if (x.Size() == 3) {
            E(2) = 0.0;
        }
    }
}

// 33. define spatially varying righthand side function
void f_exact(const Vector &x, Vector &f)
{
    if (dim == 3)
    {
        f(0) = (1. + kappa * kappa) * sin(kappa * x(1));
        f(1) = (1. + kappa * kappa) * sin(kappa * x(2));
        f(2) = (1. + kappa * kappa) * sin(kappa * x(0));
    }
    else
    {
        f(0) = (1. + kappa * kappa) * sin(kappa * x(1));
        f(1) = (1. + kappa * kappa) * sin(kappa * x(0));
        if (x.Size() == 3) {
            f(2) = 0.0;
        }
    }
}
