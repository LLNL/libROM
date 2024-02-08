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
// Offline phase: elliptic_eigenproblem_global_rom -offline -p 2 -id 0 -a 0.8
//                elliptic_eigenproblem_global_rom -offline -p 2 -id 1 -a 0.9
//                elliptic_eigenproblem_global_rom -offline -p 2 -id 2 -a 1.1
//                elliptic_eigenproblem_global_rom -offline -p 2 -id 3 -a 1.2
//
// Merge phase:   elliptic_eigenproblem_global_rom -merge -p 2 -ns 4
//
// FOM run (for error calculation):
//                elliptic_eigenproblem_global_rom -fom -p 2 -f 1.0
//
// Online phase:  elliptic_eigenproblem_global_rom -online -p 2 -f 1.0

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

    int precision = 8;
    cout.precision(precision);

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&problem, "-p", "--problem",
                   "Problem setup to use. See options in Conductivity and Potential functions.");
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

    kappa = alpha * M_PI;

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

    // 4. Refine the mesh in serial to increase the resolution. In this example
    //    we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
    //    a command-line parameter.
    for (int lev = 0; lev < ser_ref_levels; lev++)
    {
        mesh->UniformRefinement();
    }
    mesh->GetBoundingBox(bb_min, bb_max, max(order, 1));

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
    const std::string basisName = "basis";
    const std::string basisFileName = basisName + std::to_string(id);
    CAROM::Options* options;
    CAROM::BasisGenerator *generator;
    StopWatch solveTimer, assembleTimer, mergeTimer;

    // 8. Set BasisGenerator if offline
    if (offline)
    {
        options = new CAROM::Options(fespace->GetTrueVSize(), nev, nev,
                                     update_right_SV);
        generator = new CAROM::BasisGenerator(*options, isIncremental, basisFileName);
    }

    // 9. The merge phase
    if (merge)
    {
        mergeTimer.Start();
        options = new CAROM::Options(fespace->GetTrueVSize(), max_num_snapshots, nev,
                                     update_right_SV);
        generator = new CAROM::BasisGenerator(*options, isIncremental, basisName);
        for (int paramID=0; paramID<nsets; ++paramID)
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

    // 10. Set up the parallel bilinear forms a(.,.) and m(.,.) on the finite
    //     element space. The first corresponds to the Laplacian operator -Delta,
    //     while the second is a simple mass matrix needed on the right hand side
    //     of the generalized eigenvalue problem below. The boundary conditions
    //     are implemented by elimination with special values on the diagonal to
    //     shift the Dirichlet eigenvalues out of the computational range. After
    //     serial and parallel assembly we extract the corresponding parallel
    //     matrices A and M.
    ConstantCoefficient one(1.0);

    FunctionCoefficient kappa_0(Conductivity);
    FunctionCoefficient v_0(Potential);

    // Project initial conductivity and potential to visualize initial condition
    ParGridFunction c_gf(fespace);
    c_gf.ProjectCoefficient(kappa_0);

    ParGridFunction p_gf(fespace);
    p_gf.ProjectCoefficient(v_0);

    Array<int> ess_bdr;
    if (pmesh->bdr_attributes.Size())
    {
        ess_bdr.SetSize(pmesh->bdr_attributes.Max());
        ess_bdr = 1;
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
            amg->SetPrintLevel(1);
            precond = amg;
        }
        else
        {
            // TODO: preconditioners using MFEM_SUPERLU, STRUMPACK, CPARDISO
        }
        lobpcg = new HypreLOBPCG(MPI_COMM_WORLD);
        lobpcg->SetNumModes(nev);
        lobpcg->SetRandomSeed(seed);
        lobpcg->SetPreconditioner(*precond);
        lobpcg->SetMaxIter(200);
        lobpcg->SetTol(1e-8);
        lobpcg->SetPrecondUsageMode(1);
        lobpcg->SetPrintLevel(0);
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
                generator->takeSample(x.GetData(), (double)i, 0.01);
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
        CAROM::BasisReader reader(basisName);

        Vector ev;
        const CAROM::Matrix *spatialbasis;
        if (rdim != -1)
        {
            spatialbasis = reader.getSpatialBasis(0.0, rdim);
        }
        else
        {
            spatialbasis = reader.getSpatialBasis(0.0, ef);
        }

        const int numRowRB = spatialbasis->numRows();
        const int numColumnRB = spatialbasis->numColumns();
        if (myid == 0) printf("spatial basis dimension is %d x %d\n", numRowRB,
                                  numColumnRB);

        // 17. form ROM operator
        CAROM::Matrix invReducedA;
        ComputeCtAB(*A, *spatialbasis, *spatialbasis, invReducedA);
        DenseMatrix *A_mat = new DenseMatrix(invReducedA.numRows(),
                                             invReducedA.numColumns());
        A_mat->Set(1, invReducedA.getData());

        CAROM::Matrix invReducedM;
        ComputeCtAB(*M, *spatialbasis, *spatialbasis, invReducedM);
        DenseMatrix *M_mat = new DenseMatrix(invReducedM.numRows(),
                                             invReducedM.numColumns());
        M_mat->Set(1, invReducedM.getData());

        assembleTimer.Stop();

        // 18. solve ROM
        solveTimer.Start();
        // (Q^T A Q) c = \lamba (Q^T M Q) c
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
            double ev_diff_norm = sqrt(diff_ev[i] * diff_ev[i]);
            double ev_fom_norm = sqrt(ev_fom[i] * ev_fom[i]);
            if (myid == 0)
            {
                std::cout << "Relative error of ROM solution for eigenvalue " << i << " = " <<
                          ev_diff_norm / ev_fom_norm << std::endl;
            }
        }
    }

    // 19. Save the refined mesh and the modes in parallel. This output can be
    //     viewed later using GLVis: "glvis -np <np> -m mesh -g mode".
    {
        ostringstream mesh_name, mode_name;
        mesh_name << "elliptic_eigenproblem-mesh." << setfill('0') << setw(6) << myid;

        ofstream mesh_ofs(mesh_name.str().c_str());
        mesh_ofs.precision(8);
        pmesh->Print(mesh_ofs);

        for (int i=0; i < nev && i < eigenvalues.Size(); i++)
        {
            if (fom || offline) {
                // convert eigenvector from HypreParVector to ParGridFunction
                x = lobpcg->GetEigenvector(i);
            } else {
                // for online, eigenvectors are stored in evect matrix
                Vector ev;
                evect.GetRow(i, ev);
                x = ev;
            }

            mode_name << "mode_" << setfill('0') << setw(2) << i << "."
                      << setfill('0') << setw(6) << myid;

            ofstream mode_ofs(mode_name.str().c_str());
            mode_ofs.precision(8);
            x.Save(mode_ofs);
            mode_name.str("");
        }

        ofstream sol_ev_ofs(sol_ev_name.str().c_str());
        sol_ev_ofs.precision(16);
        for (int i = 0; i < nev && i < eigenvalues.Size(); ++i)
        {
            sol_ev_ofs << eigenvalues[i] << std::endl;
        }
    }

    VisItDataCollection visit_dc("EllipticEigenproblem", pmesh);
    if (visit)
    {
        visit_dc.RegisterField("InitialConductivity", &c_gf);
        visit_dc.RegisterField("InitialPotential", &p_gf);
        std::vector<ParGridFunction*> visit_evs;
        for (int i = 0; i < nev && i < eigenvalues.Size(); i++)
        {
            if (fom || offline) {
                // convert eigenvector from HypreParVector to ParGridFunction
                x = lobpcg->GetEigenvector(i);
            } else {
                // for online, eigenvectors are stored in evect matrix
                Vector ev;
                evect.GetRow(i, ev);
                x = ev;
            }
            visit_evs.push_back(new ParGridFunction(x));
            visit_dc.RegisterField("x" + std::to_string(i), visit_evs.back());
        }
        visit_dc.SetCycle(0);
        visit_dc.SetTime(0.0);
        visit_dc.Save();
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

                if (fom || offline) {
                    // convert eigenvector from HypreParVector to ParGridFunction
                    x = lobpcg->GetEigenvector(i);
                } else {
                    // for online, eigenvectors are stored in evect matrix
                    Vector ev;
                    evect.GetRow(i, ev);
                    x = ev;
                }

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

    // 20. Free the used memory.
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
    case 5:
    case 6:
        return 1.0;
    }
    return 0.0;
}

double Potential(const Vector &x)
{
    int dim = x.Size();

    // map x to the reference [-1,1] domain
    Vector X(dim);
    Vector center(dim); // center of gaussian for problem 4 and 5
    center = kappa / M_PI;
    for (int i = 0; i < dim; i++)
    {
        double mesh_center = (bb_min[i] + bb_max[i]) * 0.5;
        X(i) = 2.0 * (x(i) - mesh_center) / (bb_max[i] - bb_min[i]);
        center(i) = 2.0 * (center(i) - mesh_center) / (bb_max[i] - bb_min[i]);
    }

    switch (problem)
    {
    case 1:
    case 2:
        return 0.0;
    case 3:
        return 1.0;
    case 4:
    case 5:
        return -std::exp(-X.DistanceSquaredTo(center) / 0.01);
    case 6:
        const double D = 1.0;
        Vector neg_center(center);
        neg_center.Neg();
        return -D * (std::exp(-X.DistanceSquaredTo(center) / 0.01) + std::exp(
                         -X.DistanceSquaredTo(neg_center) / 0.01));
    }
    return 0.0;
}
