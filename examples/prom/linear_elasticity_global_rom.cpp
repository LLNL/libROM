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


cout << "All good" << endl;

return 0;
}