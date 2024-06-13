/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

//                       libROM MFEM Example: Elliptic Eigenproblem (adapted from ex11p.cpp)
//
// =================================================================================
//
// Description:  This example code demonstrates the use of MFEM and libROM to
//               define a simple projection-based reduced order model of the
//               eigenvalue problem -Delta ((1 + alpha v) u) = lambda u with homogeneous
//               Dirichlet boundary conditions where alpha is a scalar ROM parameter
//               controlling the frequency of v.
//
//               We compute a number of the lowest eigenmodes by discretizing
//               the Laplacian and Mass operators using a FE space of the
//               specified order, or an isoparametric/isogeometric space if
//               order < 1 (quadratic for quadratic curvilinear mesh, NURBS for
//               NURBS mesh, etc.)
//
// Offline phase: elliptic_eigenproblem_global_rom -offline -p 2 -id 0 -a 0.40
//                elliptic_eigenproblem_global_rom -offline -p 2 -id 1 -a 0.45
//                elliptic_eigenproblem_global_rom -offline -p 2 -id 2 -a 0.55
//                elliptic_eigenproblem_global_rom -offline -p 2 -id 3 -a 0.60
//
// Merge phase:   elliptic_eigenproblem_global_rom -merge -p 2 -ns 4
//
// FOM run (for error calculation):
//                elliptic_eigenproblem_global_rom -fom -p 2 -a 0.50
// Example output:
//   Number of unknowns: 289
//    Eigenvalue 0: 0.04533314
//    Eigenvalue 1: 0.11332411
//    Eigenvalue 2: 0.11486387
//    Eigenvalue 3: 0.18192
//    Eigenvalue 4: 0.22964377
//   Elapsed time for assembling FOM: 1.471708e-03 second
//   Elapsed time for solving FOM: 3.416310e-01 second
//
// Online phase:  elliptic_eigenproblem_global_rom -online -p 2 -a 0.50
// Example output:
//   Eigenvalue 0: = 0.048430949
//   Eigenvalue 1: = 0.12021157
//   Eigenvalue 2: = 0.12147847
//   Eigenvalue 3: = 0.19456504
//   Eigenvalue 4: = 0.24285855
//   Relative error of ROM solution for eigenvalue 0 = 0.068334316
//   Relative error of ROM solution for eigenvalue 1 = 0.060776586
//   Relative error of ROM solution for eigenvalue 2 = 0.057586388
//   Relative error of ROM solution for eigenvalue 3 = 0.069508795
//   Relative error of ROM solution for eigenvalue 4 = 0.057544708
//   Elapsed time for assembling ROM: 4.289041e-03 second
//   Elapsed time for solving ROM: 6.225410e-04 second

#include "mfem.hpp"
#include "linalg/BasisGenerator.h"
#include "linalg/BasisReader.h"
#include "linalg/Vector.h"
#include "mfem/Utilities.hpp"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

double Conductivity(const Vector &x);
double Potential(const Vector &x);
int problem = 1;
double kappa;
double gaussian_width;
Vector bb_min, bb_max;

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
    bool dirichlet = true;
    int order = 1;
    int nev = 5;
    int seed = 75;
    bool slu_solver = false;
    bool sp_solver = false;
    bool cpardiso_solver = false;

    double alpha = 1.0;
    bool visualization = true;
    bool visit = false;
    int vis_steps = 5;
    bool fom = false;
    bool offline = false;
    bool merge = false;
    bool online = false;
    int id = 0;
    int nsets = 0;
    double ef = 0.9999;
    int rdim = -1;
    int verbose_level = 0;

    int precision = 8;
    cout.precision(precision);

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&problem, "-p", "--problem",
                   "Problem setup to use. See options in Conductivity and Potential functions.");
    args.AddOption(&dirichlet, "-dir", "--dirichlet", "-neu", "--neumann",
                   "BC switch.");
    args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                   "Number of times to refine the mesh uniformly in serial.");
    args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                   "Number of times to refine the mesh uniformly in parallel.");
    args.AddOption(&order, "-o", "--order",
                   "Order (degree) of the finite elements.");
    args.AddOption(&nev, "-n", "--num-eigs",
                   "Number of desired eigenmodes.");
    args.AddOption(&seed, "-s", "--seed",
                   "Random seed used to initialize LOBPCG.");
    args.AddOption(&alpha, "-a", "--alpha",
                   "Alpha coefficient.");
    args.AddOption(&id, "-id", "--id", "Parametric id");
    args.AddOption(&nsets, "-ns", "--nset", "Number of parametric snapshot sets");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                   "--no-visualization",
                   "Enable or disable GLVis visualization.");
    args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                   "--no-visit-datafiles",
                   "Save data files for VisIt (visit.llnl.gov) visualization.");
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
    args.AddOption(&verbose_level, "-v", "--verbose",
                   "Set the verbosity level of the LOBPCG solver and preconditioner. 0 is off.");
#ifdef MFEM_USE_SUPERLU
    args.AddOption(&slu_solver, "-slu", "--superlu", "-no-slu",
                   "--no-superlu", "Use the SuperLU Solver.");
#endif
#ifdef MFEM_USE_STRUMPACK
    args.AddOption(&sp_solver, "-sp", "--strumpack", "-no-sp",
                   "--no-strumpack", "Use the STRUMPACK Solver.");
#endif
#ifdef MFEM_USE_MKL_CPARDISO
    args.AddOption(&cpardiso_solver, "-cpardiso", "--cpardiso", "-no-cpardiso",
                   "--no-cpardiso", "Use the MKL CPardiso Solver.");
#endif
    args.Parse();
    if (slu_solver && sp_solver)
    {
        if (myid == 0)
            std::cout << "WARNING: Both SuperLU and STRUMPACK have been selected,"
                      << " please choose either one." << std::endl
                      << "         Defaulting to SuperLU." << std::endl;
        sp_solver = false;
    }
    // The command line options are also passed to the STRUMPACK
    // solver. So do not exit if some options are not recognized.
    if (!sp_solver)
    {
        if (!args.Good())
        {
            if (myid == 0)
            {
                args.PrintUsage(cout);
            }
            MPI_Finalize();
            return 1;
        }
    }

    if (myid == 0)
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
    double h_min, h_max, k_min, k_max;
    mesh->GetCharacteristics(h_min, h_max, k_min, k_max);
    gaussian_width = h_max * pow(2.0, 1.0 - ser_ref_levels - par_ref_levels);

    // 4. Refine the mesh in serial to increase the resolution. In this example
    //    we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
    //    a command-line parameter.
    for (int lev = 0; lev < ser_ref_levels; lev++)
    {
        mesh->UniformRefinement();
    }

    // 5. Define a parallel mesh by a partitioning of the serial mesh. Refine
    //    this mesh further in parallel to increase the resolution. Once the
    //    parallel mesh is defined, the serial mesh can be deleted.
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
    for (int lev = 0; lev < par_ref_levels; lev++)
    {
        pmesh->UniformRefinement();
    }

    // 6. Define a parallel finite element space on the parallel mesh. Here we
    //    use continuous Lagrange finite elements of the specified order. If
    //    order < 1, we instead use an isoparametric/isogeometric space.
    FiniteElementCollection *fec;
    if (order > 0)
    {
        fec = new H1_FECollection(order, dim);
    }
    else if (pmesh->GetNodes())
    {
        fec = pmesh->GetNodes()->OwnFEC();
    }
    else
    {
        fec = new H1_FECollection(order = 1, dim);
    }
    ParFiniteElementSpace *fespace = new ParFiniteElementSpace(pmesh, fec);
    HYPRE_BigInt size = fespace->GlobalTrueVSize();
    if (myid == 0)
    {
        cout << "Number of unknowns: " << size << endl;
    }

    // 7. Initiate ROM related variables
    int max_num_snapshots = 100;
    bool update_right_SV = false;
    bool isIncremental = false;
    const std::string baseName = "elliptic_eigenproblem_";
    const std::string basis_filename = baseName + "basis";
    CAROM::Options* options;
    CAROM::BasisGenerator *generator;
    StopWatch solveTimer, assembleTimer, mergeTimer;

    // 8. Set BasisGenerator if offline
    if (offline)
    {
        options = new CAROM::Options(fespace->GetTrueVSize(), nev, nev,
                                     update_right_SV);
        std::string snapshot_basename = baseName + "par" + std::to_string(id);
        generator = new CAROM::BasisGenerator(*options, isIncremental,
                                              snapshot_basename);
    }

    // 9. The merge phase
    if (merge)
    {
        mergeTimer.Start();
        options = new CAROM::Options(fespace->GetTrueVSize(), max_num_snapshots, nev,
                                     update_right_SV);
        generator = new CAROM::BasisGenerator(*options, isIncremental, basis_filename);
        for (int paramID=0; paramID<nsets; ++paramID)
        {
            std::string snapshot_filename = baseName + "par" + std::to_string(
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

    // 10. Set up the parallel bilinear forms a(.,.) and m(.,.) on the finite
    //     element space. The first corresponds to the Laplacian operator -Delta,
    //     while the second is a simple mass matrix needed on the right hand side
    //     of the generalized eigenvalue problem below. The boundary conditions
    //     are implemented by elimination with special values on the diagonal to
    //     shift the Dirichlet eigenvalues out of the computational range. After
    //     serial and parallel assembly we extract the corresponding parallel
    //     matrices A and M.
    assembleTimer.Start();
    ConstantCoefficient one(1.0);

    kappa = alpha * M_PI;
    FunctionCoefficient kappa_0(Conductivity);
    FunctionCoefficient v_0(Potential);

    // Project conductivity and potential to visualize initial condition
    ParGridFunction c_gf(fespace);
    c_gf.ProjectCoefficient(kappa_0);

    ParGridFunction p_gf(fespace);
    p_gf.ProjectCoefficient(v_0);

    Array<int> ess_bdr;
    if (pmesh->bdr_attributes.Size())
    {
        ess_bdr.SetSize(pmesh->bdr_attributes.Max());
        ess_bdr = (dirichlet) ? 1 : 0;
    }

    ParBilinearForm *a = new ParBilinearForm(fespace);
    a->AddDomainIntegrator(new DiffusionIntegrator(kappa_0));
    a->AddDomainIntegrator(new MassIntegrator(v_0));
    a->Assemble();
    a->EliminateEssentialBCDiag(ess_bdr, 1.0);
    a->Finalize();

    ParBilinearForm *m = new ParBilinearForm(fespace);
    m->AddDomainIntegrator(new MassIntegrator(one));
    m->Assemble();
    // shift the eigenvalue corresponding to eliminated dofs to a large value
    m->EliminateEssentialBCDiag(ess_bdr, numeric_limits<double>::min());
    m->Finalize();

    HypreParMatrix *A = a->ParallelAssemble();
    HypreParMatrix *M = m->ParallelAssemble();

#if defined(MFEM_USE_SUPERLU) || defined(MFEM_USE_STRUMPACK)
    Operator * Arow = NULL;
#endif
#ifdef MFEM_USE_SUPERLU
    if (slu_solver)
    {
        Arow = new SuperLURowLocMatrix(*A);
    }
#endif
#ifdef MFEM_USE_STRUMPACK
    if (sp_solver)
    {
        Arow = new STRUMPACKRowLocMatrix(*A);
    }
#endif

    assembleTimer.Stop();

    delete a;
    delete m;

    ParGridFunction x(fespace);
    HypreLOBPCG *lobpcg;
    Array<double> eigenvalues;
    DenseMatrix evect;

    // 11. The offline phase
    if(fom || offline)
    {
        // 12. Define and configure the LOBPCG eigensolver and the BoomerAMG
        //     preconditioner for A to be used within the solver. Set the matrices
        //     which define the generalized eigenproblem A x = lambda M x.
        Solver * precond = NULL;
        if (!slu_solver && !sp_solver && !cpardiso_solver)
        {
            HypreBoomerAMG * amg = new HypreBoomerAMG(*A);
            amg->SetPrintLevel(verbose_level);
            precond = amg;
        }
        else
        {
#ifdef MFEM_USE_SUPERLU
            if (slu_solver)
            {
                SuperLUSolver * superlu = new SuperLUSolver(MPI_COMM_WORLD);
                superlu->SetPrintStatistics(verbose_level > 0 ? true : false);
                superlu->SetSymmetricPattern(true);
                superlu->SetColumnPermutation(superlu::PARMETIS);
                superlu->SetOperator(*Arow);
                precond = superlu;
            }
#endif
#ifdef MFEM_USE_STRUMPACK
            if (sp_solver)
            {
                STRUMPACKSolver * strumpack = new STRUMPACKSolver(MPI_COMM_WORLD, argc, argv);
                strumpack->SetPrintFactorStatistics(true);
                strumpack->SetPrintSolveStatistics(verbose_level > 0 ? true : false);
                strumpack->SetKrylovSolver(strumpack::KrylovSolver::DIRECT);
                strumpack->SetReorderingStrategy(strumpack::ReorderingStrategy::METIS);
                strumpack->SetMatching(strumpack::MatchingJob::NONE);
                strumpack->SetCompression(strumpack::CompressionType::NONE);
                strumpack->SetOperator(*Arow);
                strumpack->SetFromCommandLine();
                precond = strumpack;
            }
#endif
#ifdef MFEM_USE_MKL_CPARDISO
            if (cpardiso_solver)
            {
                auto cpardiso = new CPardisoSolver(A->GetComm());
                cpardiso->SetMatrixType(CPardisoSolver::MatType::REAL_STRUCTURE_SYMMETRIC);
                cpardiso->SetPrintLevel(verbose_level);
                cpardiso->SetOperator(*A);
                precond = cpardiso;
            }
#endif
        }

        lobpcg = new HypreLOBPCG(MPI_COMM_WORLD);
        lobpcg->SetNumModes(nev);
        lobpcg->SetRandomSeed(seed);
        lobpcg->SetPreconditioner(*precond);
        lobpcg->SetMaxIter(200);
        lobpcg->SetTol(1e-8);
        lobpcg->SetPrecondUsageMode(1);
        lobpcg->SetPrintLevel(verbose_level);
        lobpcg->SetMassMatrix(*M);
        lobpcg->SetOperator(*A);

        // 13. Compute the eigenmodes and extract the array of eigenvalues. Define a
        //     parallel grid function to represent each of the eigenmodes returned by
        //     the solver.
        solveTimer.Start();
        lobpcg->Solve();
        solveTimer.Stop();
        lobpcg->GetEigenvalues(eigenvalues);

        // 14. take and write snapshots for ROM
        for (int i = 0; i < nev; i++)
        {
            if (myid == 0)
            {
                std::cout << " Eigenvalue " << i << ": " << eigenvalues[i] << "\n";
            }
            if (offline)
            {
                x = lobpcg->GetEigenvector(i);
                generator->takeSample(x.GetData());
            }
        }

        if (offline)
        {
            generator->writeSnapshot();
            delete generator;
            delete options;
        }

        delete precond;
    }

    // 15. The online phase
    if (online) {
        // 16. read the reduced basis
        assembleTimer.Start();
        CAROM::BasisReader reader(basis_filename);

        Vector ev;
        const CAROM::Matrix *spatialbasis;
        if (rdim != -1)
        {
            spatialbasis = reader.getSpatialBasis(rdim);
        }
        else
        {
            spatialbasis = reader.getSpatialBasis(ef);
        }

        const int numRowRB = spatialbasis->numRows();
        const int numColumnRB = spatialbasis->numColumns();
        if (myid == 0) printf("spatial basis dimension is %d x %d\n", numRowRB,
                                  numColumnRB);

        // 17. form ROM operator
        CAROM::Matrix ReducedA;
        ComputeCtAB(*A, *spatialbasis, *spatialbasis, ReducedA);
        DenseMatrix *A_mat = new DenseMatrix(ReducedA.numRows(),
                                             ReducedA.numColumns());
        A_mat->Set(1, ReducedA.getData());

        CAROM::Matrix ReducedM;
        ComputeCtAB(*M, *spatialbasis, *spatialbasis, ReducedM);
        DenseMatrix *M_mat = new DenseMatrix(ReducedM.numRows(),
                                             ReducedM.numColumns());
        M_mat->Set(1, ReducedM.getData());

        assembleTimer.Stop();

        // 18. solve ROM
        solveTimer.Start();
        // (Q^T A Q) c = \lambda (Q^T M Q) c
        A_mat->Eigenvalues(*M_mat, ev, evect);
        solveTimer.Stop();

        if (myid == 0)
        {
            eigenvalues = Array<double>(ev.GetData(), ev.Size());
            for (int i = 0; i < ev.Size() && i < nev; i++)
            {
                std::cout << "Eigenvalue " << i << ": = " << eigenvalues[i] << "\n";
            }
        }

        DenseMatrix tmp(evect);
        evect = DenseMatrix(nev, numRowRB);
        for (int i = 0; i < ev.Size() && i < nev; i++)
        {
            Vector evector;
            tmp.GetRow(i, evector);
            CAROM::Vector evector_carom(evector.GetData(), evector.Size(), false, false);
            CAROM::Vector *ev_carom = spatialbasis->mult(evector_carom);
            evect.SetRow(i, ev_carom->getData());
            delete ev_carom;
        }

        delete spatialbasis;

        delete A_mat;
        delete M_mat;
    }

    ostringstream sol_ev_name, sol_ev_name_fom;
    if (fom || offline)
    {
        sol_ev_name << "sol_eigenvalues_fom." << setfill('0') << setw(6) << myid;
    }
    if (online)
    {
        sol_ev_name << "sol_eigenvalues." << setfill('0') << setw(6) << myid;
        sol_ev_name_fom << "sol_eigenvalues_fom." << setfill('0') << setw(6) << myid;
    }

    // 19. Save the refined mesh and the modes in parallel. This output can be
    //     viewed later using GLVis: "glvis -np <np> -m mesh -g mode".
    Vector sign_ev(nev);
    {
        ostringstream mesh_name, mode_name, mode_ref_name;
        mesh_name << "elliptic_eigenproblem-mesh." << setfill('0') << setw(6) << myid;

        ofstream mesh_ofs(mesh_name.str().c_str());
        mesh_ofs.precision(8);
        pmesh->Print(mesh_ofs);

        std::string mode_prefix = "mode_";
        if (fom || offline)
        {
            mode_prefix += "fom_";
        }
        else if (online)
        {
            mode_prefix += "rom_";
        }

        for (int i=0; i < nev && i < eigenvalues.Size(); i++)
        {
            Vector ev;
            if (fom || offline)
            {
                ev = lobpcg->GetEigenvector(i);
            }
            else
            {
                // for online, eigenvectors are stored in evect matrix
                evect.GetRow(i, ev);
            }
            Vector mode_ref(ev.Size());

            mode_ref_name.str("");
            mode_ref_name << "mode_ref_" << setfill('0') << setw(2) << i << "."
                              << setfill('0') << setw(6) << myid;
            if (offline && id == 0)
            {
                ofstream mode_ref_ofs(mode_ref_name.str().c_str());
                mode_ref_ofs.precision(16);

                // TODO: issue using .Load() function if file written with .Save()?
                //x.Save(mode_ref_ofs);
                for (int j = 0; j < x.Size(); j++)
                {
                    mode_ref_ofs << x[j] << "\n";
                }
                mode_ref_ofs.close();
            }
            else
            {
                ifstream mode_ref_ifs(mode_ref_name.str().c_str());
                mode_ref_ifs.precision(16);
                mode_ref.Load(mode_ref_ifs, evect.NumCols());
                mode_ref_ifs.close();
            }
            sign_ev[i] = (InnerProduct(mode_ref, ev) >= 0) ? 1 : -1;
            ev *= sign_ev[i];
            // convert eigenvector from HypreParVector to ParGridFunction
            x = ev;

            mode_name.str("");
            mode_name << mode_prefix << setfill('0') << setw(2) << i << "."
                      << setfill('0') << setw(6) << myid;

            ofstream mode_ofs(mode_name.str().c_str());
            mode_ofs.precision(16);

            // TODO: issue using .Load() function if file written with .Save()?
            //x.Save(mode_ofs);
            for (int j = 0; j < x.Size(); j++)
            {
                mode_ofs << x[j] << "\n";
            }
            mode_ofs.close();
            mode_name.str("");
        }

        ofstream sol_ev_ofs(sol_ev_name.str().c_str());
        sol_ev_ofs.precision(16);
        for (int i = 0; i < nev && i < eigenvalues.Size(); ++i)
        {
            sol_ev_ofs << eigenvalues[i] << std::endl;
        }
    }

    if (online)
    {
        // Initialize FOM solution
        Vector ev_fom(nev);

        ifstream fom_file;
        fom_file.open(sol_ev_name_fom.str().c_str());
        ev_fom.Load(fom_file, ev_fom.Size());
        fom_file.close();

        Vector diff_ev(nev);
        for (int i = 0; i < eigenvalues.Size() && i < nev; i++)
        {
            diff_ev[i] = ev_fom[i] - eigenvalues[i];
            if (myid == 0)
            {
                std::cout << "FOM solution for eigenvalue " << i << " = " <<
                          ev_fom[i] << std::endl;
                std::cout << "ROM solution for eigenvalue " << i << " = " <<
                          eigenvalues[i] << std::endl;
                std::cout << "Absolute error of ROM solution for eigenvalue " << i << " = " <<
                          abs(diff_ev[i]) << std::endl;
                std::cout << "Relative error of ROM solution for eigenvalue " << i << " = " <<
                          abs(diff_ev[i]) / abs(ev_fom[i]) << std::endl;
            }
        }

        // Calculate errors of eigenvectors
        ostringstream mode_name, mode_name_fom;
        Vector mode_rom(evect.NumCols());
        Vector mode_fom(evect.NumCols());
        for (int i = 0; i < eigenvalues.Size() && i < nev; i++)
        {
            mode_name << "mode_rom_" << setfill('0') << setw(2) << i << "."
                      << setfill('0') << setw(6) << myid;
            mode_name_fom << "mode_fom_" << setfill('0') << setw(2) << i << "."
                          << setfill('0') << setw(6) << myid;

            ifstream mode_rom_ifs(mode_name.str().c_str());
            mode_rom_ifs.precision(16);
            mode_rom.Load(mode_rom_ifs, evect.NumCols());
            mode_rom_ifs.close();

            ifstream mode_fom_ifs(mode_name_fom.str().c_str());
            mode_fom_ifs.precision(16);
            mode_fom.Load(mode_fom_ifs, evect.NumCols());
            mode_fom_ifs.close();

            const double fomNorm = sqrt(InnerProduct(mode_fom, mode_fom));
            mode_fom.Add(-1.0, mode_rom);
            const double diffNorm = sqrt(InnerProduct(mode_fom, mode_fom));
            if (myid == 0) std::cout << "Relative l2 error of ROM eigenvector " << i <<
                                         " = " << diffNorm /
                                         fomNorm << std::endl;

            mode_name.str("");
            mode_name_fom.str("");
        }
    }

    VisItDataCollection *visit_dc = NULL;
    if (visit)
    {
        if (offline) visit_dc = new VisItDataCollection(baseName + "offline_par" +
                    std::to_string(id), pmesh);
        else if (fom) visit_dc = new VisItDataCollection(baseName + "fom", pmesh);
        else if (online) visit_dc = new VisItDataCollection(baseName + "rom", pmesh);
        visit_dc->RegisterField("Conductivity", &c_gf);
        visit_dc->RegisterField("Potential", &p_gf);
        std::vector<ParGridFunction*> visit_evs;
        for (int i = 0; i < nev && i < eigenvalues.Size(); i++)
        {
            Vector ev;
            if (fom || offline) {
                ev = lobpcg->GetEigenvector(i);
            }
            else {
                // for online, eigenvectors are stored in evect matrix
                evect.GetRow(i, ev);
            }
            // convert eigenvector from HypreParVector to ParGridFunction
            ev *= sign_ev[i];
            x = ev;
            // convert eigenvector from HypreParVector to ParGridFunction
            visit_evs.push_back(new ParGridFunction(x));
            visit_dc->RegisterField("Eigenmode_" + std::to_string(i), visit_evs.back());
        }
        visit_dc->SetCycle(0);
        visit_dc->SetTime(0.0);
        visit_dc->Save();
        for (size_t i = 0; i < visit_evs.size(); i++)
        {
            delete visit_evs[i];
        }
    }

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
            for (int i = 0; i < nev && i < eigenvalues.Size(); i++)
            {
                if ( myid == 0 )
                {
                    cout << "Eigenmode " << i+1 << '/' << nev
                         << ", Lambda = " << eigenvalues[i] << endl;
                }

                Vector ev;
                if (fom || offline) {
                    ev = lobpcg->GetEigenvector(i);
                }
                else {
                    evect.GetRow(i, ev);
                }
                // convert eigenvector from HypreParVector to ParGridFunction
                ev *= sign_ev[i];
                x = ev;

                sout << "parallel " << num_procs << " " << myid << "\n"
                     << "solution\n" << *pmesh << x << flush
                     << "window_title 'Eigenmode " << i+1 << '/' << nev
                     << ", Lambda = " << eigenvalues[i] << "'" << endl;

                char c;
                if (myid == 0)
                {
                    cout << "press (q)uit or (c)ontinue --> " << flush;
                    cin >> c;
                }
                MPI_Bcast(&c, 1, MPI_CHAR, 0, MPI_COMM_WORLD);

                if (c != 'c')
                {
                    break;
                }
            }
            sout.close();
        }
    }

    // 20. print timing info
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

    // 21. Free the used memory.
    if (fom || offline)
    {
        delete lobpcg;
    }
    delete M;
    delete A;

    delete fespace;
    if (order > 0)
    {
        delete fec;
    }
    delete pmesh;

#if defined(MFEM_USE_SUPERLU) || defined(MFEM_USE_STRUMPACK)
    delete Arow;
#endif

    MPI_Finalize();

    return 0;
}

double Conductivity(const Vector &x)
{
    int dim = x.Size();

    switch (problem)
    {
    case 1:
        return 1.0;
    case 2:
        return 1.0 + cos(kappa * x(1)) * sin(kappa * x(0));
    case 3:
    case 4:
        return 1.0;
    }
    return 0.0;
}

double Potential(const Vector &x)
{
    double D = 100.0;
    Vector center(x.Size());
    switch (problem)
    {
    case 1:
    case 2:
        return 0.0;
    case 3:
        for (int i = 0; i < x.Size(); i++)
        {
            center(i) = 0.5 * (bb_min[i] + bb_max[i]);
        }
        return D * std::exp(-x.DistanceSquaredTo(center) / pow(gaussian_width, 2.0));
    case 4:
        for (int i = 0; i < x.Size(); i++)
        {
            center(i) = 0.5 * (bb_min[i] + bb_max[i]);
        }
        return D / x.DistanceSquaredTo(center);
    }
    return 0.0;
}
