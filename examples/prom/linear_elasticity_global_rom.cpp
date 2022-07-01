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


    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
        "Mesh file to use.");
    args.AddOption(&order, "-o", "--order",
        "Finite element order (polynomial degree).");
    args.AddOption(&amg_elast, "-elast", "--amg-for-elasticity", "-sys",
        "--amg-for-systems",
        "Use the special AMG elasticity solver (GM/LN approaches), "
        "or standard AMG for systems (unknown approach).");
    args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
        "--no-static-condensation", "Enable static condensation.");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
        "--no-visualization",
        "Enable or disable GLVis visualization.");
    args.AddOption(&reorder_space, "-nodes", "--by-nodes", "-vdim", "--by-vdim",
        "Use byNODES ordering of vector space instead of byVDIM");
    args.AddOption(&device_config, "-d", "--device",
        "Device configuration string, see Device::Configure().");
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


    // 5. Select the order of the finite element discretization space. For NURBS
    //    meshes, we increase the order by degree elevation.
    if (mesh->NURBSext)
    {
        mesh->DegreeElevate(order, order);
    }


    // 6. Refine the serial mesh on all processors to increase the resolution. In
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


    // 7. Define a parallel mesh by a partitioning of the serial mesh. Refine
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



    // 8. Define a parallel finite element space on the parallel mesh. Here we
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



    // 9. Determine the list of true (i.e. parallel conforming) essential
    //    boundary dofs. In this example, the boundary conditions are defined by
    //    marking only boundary attribute 1 from the mesh as essential and
    //    converting it to a list of true dofs.
    Array<int> ess_tdof_list, ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[0] = 1;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);


cout << "All good" << endl;

return 0;
}
