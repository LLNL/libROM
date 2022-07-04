//               libROM MFEM Example: parametric ROM for linear elastic problem (adapted from ex2p.cpp)
//
// Compile with: ./scripts/compile.sh -m
//
// Description:  This example code demonstrates the use of MFEM and libROM to
//               define a simple projection-based reduced order model of a 
//				 simple linear elasticity problem describing a multi-material 
//				 cantilever beam.
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
// Offline phase: 
// 
//
// Merge phase:   
//
// Online phase:  

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include "linalg/BasisGenerator.h"
#include "linalg/BasisReader.h"

using namespace std;
using namespace mfem;


int main(int argc, char* argv[])
{
    // 1. Initialize MPI.
    int num_procs, myid;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // 2. Parse command-line options.
    const char* mesh_file = "../data/beam-tri.mesh";
    int order = 1;
    bool static_cond = false;
    bool visualization = 1;
    bool amg_elast = 0;
    bool reorder_space = false;
    const char* device_config = "cpu";
    bool visit = true;
    bool fom = false;
    bool offline = false;
    bool merge = false;
    bool online = false;
    int precision = 8;
    int id = 0;
    int nsets = 0;
    double coef = 1.0;
    double ext_force = -1.0e-2;


    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
        "Mesh file to use.");
    args.AddOption(&order, "-o", "--order",
        "Finite element order (polynomial degree).");
    args.AddOption(&id, "-id", "--id", "Parametric id");
    args.AddOption(&nsets, "-ns", "--nset", "Number of parametric snapshot sets");
    args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
        "--no-static-condensation", "Enable static condensation.");
    args.AddOption(&ext_force, "-f", "--ext-force",
        "External force applied at end.");
    args.AddOption(&amg_elast, "-elast", "--amg-for-elasticity", "-sys",
        "--amg-for-systems",
        "Use the special AMG elasticity solver (GM/LN approaches), "
        "or standard AMG for systems (unknown approach).");
    args.AddOption(&reorder_space, "-nodes", "--by-nodes", "-vdim", "--by-vdim",
        "Use byNODES ordering of vector space instead of byVDIM");
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
        return 1;
    }
    if (myid == 0)
    {
        args.PrintOptions(cout);
    }

    // 3. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
    Device device(device_config);
    if (myid == 0) { device.Print(); }


    // 4. Read the (serial) mesh from the given mesh file on all processors.  We
   //    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
   //    and volume meshes with the same code.
    Mesh* mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    if (mesh->attributes.Max() < 2 || mesh->bdr_attributes.Max() < 2)
    {
        if (myid == 0)
            cerr << "\nInput mesh should have at least two materials and "
            << "two boundary attributes! (See schematic in ex2.cpp)\n"
            << endl;
        return 3;
    }



    // 5. Refine the serial mesh on all processors to increase the resolution. In
    //    this example we do 'ref_levels' of uniform refinement. We choose
    //    'ref_levels' to be the largest number that gives a final mesh with no
    //    more than 1,000 elements.
    {
        int ref_levels =
            (int)floor(log(1000. / mesh->GetNE()) / log(2.) / dim);
        for (int l = 0; l < ref_levels; l++)
        {
            mesh->UniformRefinement();
        }
    }


    // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution. Once the
   //    parallel mesh is defined, the serial mesh can be deleted.
    ParMesh* pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
    {
        int par_ref_levels = 1;
        for (int l = 0; l < par_ref_levels; l++)
        {
            pmesh->UniformRefinement();
        }
    }



    // 7. Define a parallel finite element space on the parallel mesh. Here we
    //    use vector finite elements, i.e. dim copies of a scalar finite element
    //    space. We use the ordering by vector dimension (the last argument of
    //    the FiniteElementSpace constructor) which is expected in the systems
    //    version of BoomerAMG preconditioner. For NURBS meshes, we use the
    //    (degree elevated) NURBS space associated with the mesh nodes.
    FiniteElementCollection* fec;
    ParFiniteElementSpace* fespace;
    const bool use_nodal_fespace = pmesh->NURBSext && !amg_elast;
    if (use_nodal_fespace)
    {
        fec = NULL;
        fespace = (ParFiniteElementSpace*)pmesh->GetNodes()->FESpace();
    }
    else
    {
        fec = new H1_FECollection(order, dim);
        if (reorder_space)
        {
            fespace = new ParFiniteElementSpace(pmesh, fec, dim, Ordering::byNODES);
        }
        else
        {
            fespace = new ParFiniteElementSpace(pmesh, fec, dim, Ordering::byVDIM);
        }
    }
    HYPRE_BigInt size = fespace->GlobalTrueVSize();
    if (myid == 0)
    {
        cout << "Number of finite element unknowns: " << size << endl
            << "Assembling: " << flush;
    }



    // 8. Determine the list of true (i.e. parallel conforming) essential
    //    boundary dofs. In this example, the boundary conditions are defined by
    //    marking only boundary attribute 1 from the mesh as essential and
    //    converting it to a list of true dofs.
    Array<int> ess_tdof_list, ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[0] = 1;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);



    // 9. Initiate ROM related variables
    int max_num_snapshots = 100;
    bool update_right_SV = false;
    bool isIncremental = false;
    const std::string basisName = "basis";
    const std::string basisFileName = basisName + std::to_string(id);
    const CAROM::Matrix* spatialbasis;
    CAROM::Options* options;
    CAROM::BasisGenerator* generator;
    int numRowRB, numColumnRB;
    StopWatch solveTimer, assembleTimer, mergeTimer;


    // 10. Set BasisGenerator if offline
    if (offline)
    {
        options = new CAROM::Options(fespace->GetTrueVSize(), max_num_snapshots, 1,
            update_right_SV);
        generator = new CAROM::BasisGenerator(*options, isIncremental, basisFileName);
    }

    // 11. The merge phase
    if (merge)
    {
        mergeTimer.Start();
        std::unique_ptr<CAROM::BasisGenerator> basis_generator;
        options = new CAROM::Options(fespace->GetTrueVSize(), max_num_snapshots, 1,
            update_right_SV);
        generator = new CAROM::BasisGenerator(*options, isIncremental, basisName);
        for (int paramID = 0; paramID < nsets; ++paramID)
        {
            std::string snapshot_filename = basisName + std::to_string(
                paramID) + "_snapshot";
            generator->loadSamples(snapshot_filename, "snapshot");
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
    //     right-hand side of the FEM linear system. In this case, b_i equals the
    //     boundary integral of f*phi_i where f represents a "pull down" force on
    //     the Neumann part of the boundary and phi_i are the basis functions in
    //     the finite element fespace. The force is defined by the object f, which
    //     is a vector of Coefficient objects. The fact that f is non-zero on
    //     boundary attribute 2 is indicated by the use of piece-wise constants
    //     coefficient for its last component.
    VectorArrayCoefficient f(dim);
    for (int i = 0; i < dim - 1; i++)
    {
        f.Set(i, new ConstantCoefficient(0.0));
    }
    {
        Vector pull_force(pmesh->bdr_attributes.Max());
        pull_force = 0.0;
        pull_force(1) = ext_force;
        f.Set(dim - 1, new PWConstCoefficient(pull_force));
    }

    ParLinearForm* b = new ParLinearForm(fespace);
    b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));
    if (myid == 0)
    {
        cout << "r.h.s. ... " << flush;
    }
    b->Assemble();



    // 13. Define the solution vector x as a parallel finite element grid
    //     function corresponding to fespace. Initialize x with initial guess of
    //     zero, which satisfies the boundary conditions.
    ParGridFunction x(fespace);
    x = 0.0;



    // 14. Set up the parallel bilinear form a(.,.) on the finite element space
    //     corresponding to the linear elasticity integrator with piece-wise
    //     constants coefficient lambda and mu.
    assembleTimer.Start();
    Vector lambda(pmesh->attributes.Max());
    lambda = 1.0;
    lambda(0) = lambda(1) * 50;
    PWConstCoefficient lambda_func(lambda);
    Vector mu(pmesh->attributes.Max());
    mu = 1.0;
    mu(0) = mu(1) * 50;
    PWConstCoefficient mu_func(mu);

    ParBilinearForm* a = new ParBilinearForm(fespace);
    a->AddDomainIntegrator(new ElasticityIntegrator(lambda_func, mu_func));



    // 15. Assemble the parallel bilinear form and the corresponding linear
    //     system, applying any necessary transformations such as: parallel
    //     assembly, eliminating boundary conditions, applying conforming
    //     constraints for non-conforming AMR, static condensation, etc.
    if (myid == 0) { cout << "matrix ... " << flush; }
    if (static_cond) { a->EnableStaticCondensation(); }
    a->Assemble();

    HypreParMatrix A;
    Vector B, X;
    a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
    if (myid == 0)
    {
        cout << "done." << endl;
        cout << "Size of linear system: " << A.GetGlobalNumRows() << endl;
    }
    assembleTimer.Stop();

    // 16. The offline phase
    if (fom || offline)
    {
        // 17. Define and apply a parallel PCG solver for A X = B with the BoomerAMG
        //     preconditioner from hypre.
        HypreBoomerAMG* amg = new HypreBoomerAMG(A);
        if (amg_elast && !a->StaticCondensationIsEnabled())
        {
            amg->SetElasticityOptions(fespace);
        }
        else
        {
            amg->SetSystemsOptions(dim, reorder_space);
        }
        HyprePCG* pcg = new HyprePCG(A);
        pcg->SetTol(1e-8);
        pcg->SetMaxIter(500);
        pcg->SetPrintLevel(2);
        pcg->SetPreconditioner(*amg);
        solveTimer.Start();
        pcg->Mult(B, X);
        solveTimer.Stop();
        delete amg;
        delete pcg;

        // 18. take and write snapshot for ROM
        if (offline)
        {
            bool addSample = generator->takeSample(X.GetData(), 0.0, 0.01);
            generator->writeSnapshot();
            delete generator;
            delete options;
        }

    }

    // 19 - 23. The online phase
    // TODO
    


    // 24. Recover the parallel grid function corresponding to X. This is the
    //     local finite element solution on each processor.
    a->RecoverFEMSolution(X, *b, x);



    // 16. For non-NURBS meshes, make the mesh curved based on the finite element
    //     space. This means that we define the mesh elements through a fespace
    //     based transformation of the reference element.  This allows us to save
    //     the displaced mesh as a curved mesh when using high-order finite
    //     element displacement field. We assume that the initial mesh (read from
    //     the file) is not higher order curved mesh compared to the chosen FE
    //     space.
    if (!use_nodal_fespace)
    {
        pmesh->SetNodalFESpace(fespace);
    }



    // 25. Save in parallel the displaced mesh and the inverted solution (which
    //     gives the backward displacements to the original grid). This output
    //     can be viewed later using GLVis: "glvis -np <np> -m mesh -g sol".
    {
        GridFunction* nodes = pmesh->GetNodes();
        *nodes += x;
        x *= -1;

        ostringstream mesh_name, sol_name;
        mesh_name << "mesh." << setfill('0') << setw(6) << myid;
        sol_name << "sol." << setfill('0') << setw(6) << myid;

        ofstream mesh_ofs(mesh_name.str().c_str());
        mesh_ofs.precision(precision);
        pmesh->Print(mesh_ofs);

        ofstream sol_ofs(sol_name.str().c_str());
        sol_ofs.precision(precision);
        x.Save(sol_ofs);
    }


    // 26. Save data in the VisIt format.
    // TODO



    // 27. Send the above data by socket to a GLVis server.  Use the "n" and "b"
   //     keys in GLVis to visualize the displacements.
    if (visualization)
    {
        char vishost[] = "localhost";
        int  visport = 19916;
        socketstream sol_sock(vishost, visport);
        sol_sock << "parallel " << num_procs << " " << myid << "\n";
        sol_sock.precision(precision);
        sol_sock << "solution\n" << *pmesh << x << flush;
    }


    // 28. print timing info
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

    // 29. Free the used memory.
    delete a;
    delete b;
    if (fec)
    {
        delete fespace;
        delete fec;
    }
    delete pmesh;

cout << "All good" << endl;

return 0;
}


