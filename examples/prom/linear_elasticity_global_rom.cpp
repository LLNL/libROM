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


cout << "All good" << endl;

return 0;
}
