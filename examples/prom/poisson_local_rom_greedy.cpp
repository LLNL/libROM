/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

//               libROM MFEM Example: parametric ROM for Poisson problem (adapted from ex1p.cpp)
//
// Compile with: ./scripts/compile.sh -m
//
// Description:  This example code demonstrates the use of MFEM and libROM to
//               define a simple projection-based reduced order model of the
//               Poisson problem -Delta u = f(x) with homogeneous Dirichlet
//               boundary conditions and spatially varying right hand side f.
//
//               The example highlights the greedy algorithm. The build_database phase
//               builds a global ROM database using different frequncies and a latin-hypercube
//               sampling procedure. The use_database phase uses the global ROM database,
//               builds the ROM operator, solves thereduced order system, and
//               lifts the solution to the full order space.
//
// build_database phase: poisson_local_rom_greedy -build_database -greedy-param-min 1.0 -greedy-param-max 1.2 -greedy-param-size 5 -greedysubsize 2 -greedyconvsize 3 -greedyrelerrortol 0.01
// use_database phase:   poisson_local_rom_greedy -fom -f 1.15 (create a new solution to compare with)
// use_database phase:   poisson_local_rom_greedy -use_database -online -f 1.15 (use the database to compute at f 1.15 while comparing to the true offline solution at f 1.15)
//
// Larger example:
// build_database phase: poisson_local_rom_greedy -build_database -greedy-param-min 0.5 -greedy-param-max 1.5 -greedy-param-size 15 -greedysubsize 4 -greedyconvsize 6 -greedyrelerrortol 0.01
// use_database phase:   poisson_local_rom_greedy -fom -f X.XX (create a new solution to compare with. Set X.XX to your desired frequency.)
// use_database phase:   poisson_local_rom_greedy -use_database -online -f X.XX (use the database to compute at f X.XX while comparing to the true offline solution at f X.XX)
//
// This example runs in parallel with MPI, by using the same number of MPI ranks
// in all phases (build_database, fom, online).

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include "algo/greedy/GreedyRandomSampler.h"
#include "linalg/BasisGenerator.h"
#include "linalg/BasisReader.h"
#include "mfem/Utilities.hpp"

using namespace std;
using namespace mfem;

// rhs function. See below for implementation.
double rhs(const Vector &);
double freq = 1.0, kappa;
int dim;

int main(int argc, char *argv[])
{
    // 1. Initialize MPI.
    int num_procs, myid;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // 2. Parse command-line options.
    const char *mesh_file = "../data/star.mesh";
    int order = 1;
    bool static_cond = false;
    bool pa = false;
    const char *device_config = "cpu";
    bool visualization = false;
    bool visit = false;
    bool build_database = false;
    bool use_database = false;
    bool fom = false;
    bool offline = false;
    bool online = false;
    double greedy_param_space_min = 1.0;
    double greedy_param_space_max = 1.0;
    int greedy_param_space_size = 0;
    double greedy_relative_error_tol = 0.01;
    int greedy_subset_size = 0;
    int greedy_convergence_subset_size = 0;
    int precision = 16;
    double coef = 1.0;

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&order, "-o", "--order",
                   "Finite element order (polynomial degree) or -1 for"
                   " isoparametric space.");
    args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                   "--no-static-condensation", "Enable static condensation.");
    args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                   "--no-partial-assembly", "Enable Partial Assembly.");
    args.AddOption(&freq, "-f", "--frequency", "Set the frequency for the exact"
                   " solution.");
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
    args.AddOption(&build_database, "-build_database", "--build_database",
                   "-no-build_database", "--no-build_database",
                   "Enable or disable the build_database phase of the greedy algorithm.");
    args.AddOption(&use_database, "-use_database", "--use_database",
                   "-no-use_database", "--no-use_database",
                   "Enable or disable the use_database phase of the greedy algorithm.");
    args.AddOption(&fom, "-fom", "--fom", "-no-fom", "--no-fom",
                   "Enable or disable the fom phase.");
    args.AddOption(&offline, "-offline", "--offline", "-no-offline", "--no-offline",
                   "Enable or disable the offline phase.");
    args.AddOption(&online, "-online", "--online", "-no-online", "--no-online",
                   "Enable or disable the online phase.");
    args.AddOption(&greedy_param_space_min, "-greedy-param-min",
                   "--greedy-param-min", "The minimum value of the parameter point space.");
    args.AddOption(&greedy_param_space_max, "-greedy-param-max",
                   "--greedy-param-max", "The maximum value of the parameter point space.");
    args.AddOption(&greedy_param_space_size, "-greedy-param-size",
                   "--greedy-param-size",
                   "The number of values to search in the parameter point space.");
    args.AddOption(&greedy_relative_error_tol, "-greedyrelerrortol",
                   "--greedyrelerrortol", "The greedy algorithm relative error tolerance.");
    args.AddOption(&greedy_subset_size, "-greedysubsize", "--greedysubsize",
                   "The greedy algorithm subset size.");
    args.AddOption(&greedy_convergence_subset_size, "-greedyconvsize",
                   "--greedyconvsize", "The greedy algorithm convergence subset size.");
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

    if (fom) MFEM_VERIFY(fom && !build_database && !use_database && !offline
                             && !online, "everything must be turned off if fom is used.");

    CAROM::GreedySampler* greedy_sampler = NULL;
    MFEM_VERIFY(!build_database
                || !use_database,
                "both build_database and use_database can not be used at the same time.");

    // 3. Set up the ROM database for the greedy algorithm to run.
    if (build_database)
    {
        MFEM_VERIFY(!offline
                    && !online,
                    "offline and online must be turned off during the build_database phase.");
        MFEM_VERIFY(!visit
                    && !visualization,
                    "visit and visualization must be turned off during the build_database phase.")
        std::ifstream infile("poisson_local_rom_greedy_algorithm_data");
        if (infile.good())
        {
            if (myid == 0) cout << "The database has already been built."
                                    << "Exiting." << endl;
            return 0;
        }
        infile.close();
        greedy_sampler = new CAROM::GreedyRandomSampler(greedy_param_space_min,
                greedy_param_space_max,
                greedy_param_space_size, false, greedy_relative_error_tol, 1.05,
                2.0, greedy_subset_size, greedy_convergence_subset_size,
                true, "poisson_local_rom_greedy_algorithm_log.txt");
    }
    // 3. Or use the database set up by the greedy algorithm.
    else if (use_database)
    {
        MFEM_VERIFY(!offline
                    && online,
                    "offline must be turned off and online must be turned on during the build_database phase.");
        std::ifstream infile("poisson_local_rom_greedy_algorithm_data");
        if (!infile.good())
        {
            if (myid == 0) cout << "The database has not been built. Exiting."
                                    << endl;
            return 0;
        }
        infile.close();
    }

    vector<string> basisIdentifiers;
    // The simulation is wrapped in a do-while statement so that the greedy
    // algorithm (build_database) can run multiple simulations in succession.
    do
    {
        if (build_database)
        {
            offline = false;
            online = false;
        }
        bool calc_rel_error = false;
        bool calc_err_indicator = false;
        std::string curr_basis_identifier = "";

        // 4. Set the correct stage of the greedy algorithm (i.e. sampling new point,
        //    calculating relative error of the last sampled point, or calculating
        //    an error indicator at a new point.)
        if (build_database)
        {
            double local_rom_freq = 0.0;
            double curr_freq = 0.0;
            struct CAROM::GreedyErrorIndicatorPoint pointRequiringRelativeError =
                greedy_sampler->getNextPointRequiringRelativeError();
            CAROM::Vector* relativeErrorPointData = pointRequiringRelativeError.point.get();
            struct CAROM::GreedyErrorIndicatorPoint pointRequiringErrorIndicator =
                greedy_sampler->getNextPointRequiringErrorIndicator();
            CAROM::Vector* errorIndicatorPointData =
                pointRequiringErrorIndicator.point.get();

            if (relativeErrorPointData != NULL)
            {
                if (myid == 0) cout << "Calculating the relative error of "
                                        << "the last sampled point." << endl;

                local_rom_freq = pointRequiringRelativeError.localROM.get()->item(0);
                curr_freq = pointRequiringRelativeError.point.get()->item(0);
                if (myid == 0) cout << "Using the basis obtained at the frequency: " <<
                                        local_rom_freq << endl;
                online = true;
                calc_rel_error = true;
            }
            else if (errorIndicatorPointData != NULL)
            {
                if (myid == 0) cout << "Calculating an error indicator at "
                                        << "a new point." << endl;

                local_rom_freq = pointRequiringErrorIndicator.localROM.get()->item(0);
                curr_freq = pointRequiringErrorIndicator.point.get()->item(0);
                if (myid == 0) cout << "Using the basis obtained at the frequency: "
                                        << local_rom_freq << endl;
                online = true;
                calc_err_indicator = true;
            }
            else
            {
                std::shared_ptr<CAROM::Vector> nextSampleParameterPoint =
                    greedy_sampler->getNextParameterPoint();
                CAROM::Vector* samplePointData = nextSampleParameterPoint.get();
                if (samplePointData != NULL)
                {
                    if (myid == 0) cout << "Sampling a new point." << endl;
                    local_rom_freq = samplePointData->item(0);
                    curr_freq = samplePointData->item(0);
                    offline = true;
                }
                else
                {
                    if (myid == 0) cout << "The greedy algorithm has finished." << endl;
                    greedy_sampler->save("poisson_local_rom_greedy_algorithm_data");
                    build_database = false;
                    continue;
                }
            }
            // 4a. Set the correct frequency as commanded by the greedy algorithm.
            curr_basis_identifier += "_" + to_string(curr_freq);
            freq = curr_freq;
        }

        if (myid == 0) cout << "Running loop at frequency: " << freq << endl;

        kappa = freq * M_PI;

        // 5. Enable hardware devices such as GPUs, and programming models such as
        //    CUDA, OCCA, RAJA and OpenMP based on command line options.
        Device device(device_config);
        if (myid == 0) {
            device.Print();
        }

        // 6. Read the (serial) mesh from the given mesh file on all processors.  We
        //    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
        //    and volume meshes with the same code.
        Mesh mesh(mesh_file, 1, 1);
        dim = mesh.Dimension();

        // 7. Refine the serial mesh on all processors to increase the resolution. In
        //    this example we do 'ref_levels' of uniform refinement. We choose
        //    'ref_levels' to be the largest number that gives a final mesh with no
        //    more than 10,000 elements.
        {
            int ref_levels =
                (int)floor(log(10000./mesh.GetNE())/log(2.)/dim);
            for (int l = 0; l < ref_levels; l++)
            {
                mesh.UniformRefinement();
            }
        }
        // 8. Define a parallel mesh by a partitioning of the serial mesh. Refine
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

        // 9. Define a parallel finite element space on the parallel mesh. Here we
        //    use continuous Lagrange finite elements of the specified order. If
        //    order < 1, we instead use an isoparametric/isogeometric space.
        FiniteElementCollection *fec;
        bool delete_fec;
        if (order > 0)
        {
            fec = new H1_FECollection(order, dim);
            delete_fec = true;
        }
        else if (pmesh.GetNodes())
        {
            fec = pmesh.GetNodes()->OwnFEC();
            delete_fec = false;
            if (myid == 0)
            {
                cout << "Using isoparametric FEs: " << fec->Name() << endl;
            }
        }
        else
        {
            fec = new H1_FECollection(order = 1, dim);
            delete_fec = true;
        }
        ParFiniteElementSpace fespace(&pmesh, fec);
        HYPRE_Int globalSize = fespace.GlobalTrueVSize();
        if (myid == 0)
        {
            cout << "Number of finite element unknowns: " << globalSize << endl;
        }

        // 10. Determine the list of true (i.e. parallel conforming) essential
        //    boundary dofs. In this example, the boundary conditions are defined
        //    by marking all the boundary attributes from the mesh as essential
        //    (Dirichlet) and converting them to a list of true dofs.
        Array<int> ess_tdof_list;
        if (pmesh.bdr_attributes.Size())
        {
            Array<int> ess_bdr(pmesh.bdr_attributes.Max());
            ess_bdr = 1;
            fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
        }

        // 11. Initiate ROM related variables
        int max_num_snapshots = 100;
        bool update_right_SV = false;
        bool isIncremental = false;
        const std::string saveBasisName = "basis" + curr_basis_identifier;
        std::string loadBasisName = "basis";

        // 11a. If using the greedy algorithm, load the global greedy basis.
        if (build_database || use_database)
        {
            loadBasisName += "_greedy";
        }

        const CAROM::Matrix* spatialbasis;
        CAROM::Options* options;
        CAROM::BasisGenerator *generator;
        int numRowRB, numColumnRB;
        StopWatch solveTimer, assembleTimer, mergeTimer;

        // 12. Set BasisGenerator if offline
        if (offline)
        {
            options = new CAROM::Options(fespace.GetTrueVSize(), max_num_snapshots, 1,
                                         update_right_SV);
            if (myid == 0) cout << "Saving basis to: " << saveBasisName << endl;
            generator = new CAROM::BasisGenerator(*options, isIncremental, saveBasisName);
        }

        // 13. Set up the parallel linear form b(.) which corresponds to the
        //     right-hand side of the FEM linear system, which in this case is
        //     (f,phi_i) where f is given by the function f_exact and phi_i are the
        //     basis functions in the finite element fespace.
        assembleTimer.Start();
        ParLinearForm *b = new ParLinearForm(&fespace);
        FunctionCoefficient f(rhs);
        b->AddDomainIntegrator(new DomainLFIntegrator(f));
        b->Assemble();

        // 14. Define the solution vector x as a parallel finite element grid function
        //     corresponding to fespace. Initialize x with initial guess of zero,
        //     which satisfies the boundary conditions.
        ParGridFunction x(&fespace);
        x = 0.0;

        // 15. Set up the parallel bilinear form a(.,.) on the finite element space
        //     corresponding to the Laplacian operator -Delta, by adding the Diffusion
        //     domain integrator.
        ParBilinearForm a(&fespace);
        ConstantCoefficient one(coef);
        if (pa) {
            a.SetAssemblyLevel(AssemblyLevel::PARTIAL);
        }
        a.AddDomainIntegrator(new DiffusionIntegrator(one));

        // 16. Assemble the parallel bilinear form and the corresponding linear
        //     system, applying any necessary transformations such as: parallel
        //     assembly, eliminating boundary conditions, applying conforming
        //     constraints for non-conforming AMR, static condensation, etc.
        if (static_cond) {
            a.EnableStaticCondensation();
        }
        a.Assemble();

        HypreParMatrix A;
        Vector B, X;
        a.FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
        assembleTimer.Stop();

        // 17. The offline phase
        if(fom || offline)
        {
            // 18. Solve the full order linear system A X = B
            Solver *prec = NULL;
            if (pa)
            {
                if (UsesTensorBasis(fespace))
                {
                    prec = new OperatorJacobiSmoother(a, ess_tdof_list);
                }
            }
            else
            {
                prec = new HypreBoomerAMG;
            }
            CGSolver cg(MPI_COMM_WORLD);
            cg.SetRelTol(1e-12);
            cg.SetMaxIter(2000);
            cg.SetPrintLevel(1);
            if (prec) {
                cg.SetPreconditioner(*prec);
            }
            cg.SetOperator(A);
            solveTimer.Start();
            cg.Mult(B, X);
            solveTimer.Stop();
            delete prec;

            // 19. take and write snapshot for ROM
            if (offline)
            {
                bool addSample = generator->takeSample(X.GetData(), 0.0, 0.01);
                generator->writeSnapshot();
                basisIdentifiers.push_back(saveBasisName);
                delete generator;
                delete options;
            }
        }

        // 20. The online phase
        if (online) {
            // 21. read the reduced basis
            assembleTimer.Start();
            CAROM::BasisReader reader(loadBasisName);
            spatialbasis = reader.getSpatialBasis(0.0);
            numRowRB = spatialbasis->numRows();
            numColumnRB = spatialbasis->numColumns();
            if (myid == 0) printf("spatial basis dimension is %d x %d\n", numRowRB,
                                      numColumnRB);

            // 22. form inverse ROM operator
            CAROM::Matrix invReducedA(numColumnRB, numColumnRB, false);
            ComputeCtAB(A, *spatialbasis, *spatialbasis, invReducedA);
            invReducedA.inverse();

            CAROM::Vector B_carom(B.GetData(), B.Size(), true, false);
            CAROM::Vector X_carom(X.GetData(), X.Size(), true, false);
            CAROM::Vector *reducedRHS = spatialbasis->transposeMult(&B_carom);
            CAROM::Vector reducedSol(numColumnRB, false);
            assembleTimer.Stop();

            // 23. solve ROM
            solveTimer.Start();
            invReducedA.mult(*reducedRHS, reducedSol);
            solveTimer.Stop();

            // 24. reconstruct FOM state
            spatialbasis->mult(reducedSol, X_carom);
            delete spatialbasis;
            delete reducedRHS;
        }

        // 25. Recover the parallel grid function corresponding to X. This is the
        //     local finite element solution on each processor.
        a.RecoverFEMSolution(X, *b, x);

        // 26. Save the refined mesh and the solution in parallel. This output can
        //     be viewed later using GLVis: "glvis -np <np> -m mesh -g sol".
        if (fom || offline)
        {
            ostringstream mesh_name, sol_name;
            mesh_name << "Mesh" << curr_basis_identifier << myid;
            sol_name << "Sol" << curr_basis_identifier << myid;

            ofstream mesh_ofs(mesh_name.str().c_str());
            mesh_ofs.precision(precision);
            pmesh.Print(mesh_ofs);

            ofstream sol_ofs(sol_name.str().c_str());
            if (myid == 0) cout << "Saving solution to: " << sol_name.str() << endl;

            sol_ofs.precision(precision);
            for (int i = 0; i < X.Size(); ++i)
            {
                sol_ofs << X[i] << endl;
            }
        }

        // 27. Save data in the VisIt format.
        DataCollection *dc = NULL;
        if (visit)
        {
            if(offline || fom) dc = new VisItDataCollection("Example1", &pmesh);
            else if(online) dc = new VisItDataCollection("Example1_rom", &pmesh);
            dc->SetPrecision(precision);
            dc->RegisterField("solution", &x);
            dc->Save();
            delete dc;
        }

        // 28. Send the solution by socket to a GLVis server.
        if (visualization)
        {
            char vishost[] = "localhost";
            int  visport   = 19916;
            socketstream sol_sock(vishost, visport);
            sol_sock << "parallel " << num_procs << " " << myid << "\n";
            sol_sock.precision(precision);
            sol_sock << "solution\n" << pmesh << x << flush;
        }

        double curr_error = 0;

        // 29. Calculate the relative error as commanded by the greedy algorithm.
        if (calc_rel_error)
        {
            Vector true_solution(X.Size());
            ifstream solution_file;
            std::string solution_filename = "Sol" + curr_basis_identifier + to_string(myid);
            if (myid == 0) cout << "Comparing current run to solution at: "
                                    << solution_filename << endl;
            solution_file.open(solution_filename);
            true_solution.Load(solution_file, X.Size());
            solution_file.close();
            Vector residual(X.Size());
            subtract(X, true_solution, residual);

            const double abs_error = sqrt(InnerProduct(MPI_COMM_WORLD, residual, residual));
            const double sol_norm = sqrt(InnerProduct(MPI_COMM_WORLD, true_solution,
                                         true_solution));

            curr_error = abs_error / sol_norm;
            if (myid == 0) cout << "The relative error is: " << curr_error << endl;

            greedy_sampler->setPointRelativeError(curr_error);
        }
        // 29. Or calculate the error indicator as commanded by the greedy algorithm.
        else if (calc_err_indicator)
        {
            Vector AX(X.Size());
            A.Mult(X, AX);
            Vector residual(X.Size());
            subtract(B, AX, residual);

            curr_error = sqrt(InnerProduct(MPI_COMM_WORLD, residual, residual));
            if (myid == 0) cout << "The error indicator is: "
                                    << curr_error << endl;
            greedy_sampler->setPointErrorIndicator(curr_error, 1);
        }
        // 29. Or if not using the greedy algorithm, output the relative error
        //     the error indicator of the point. (just for extra information)
        else if (!build_database && online)
        {
            Vector true_solution(X.Size());
            ifstream solution_file;
            std::string solution_filename = "Sol" + curr_basis_identifier + to_string(myid);
            if (myid == 0) cout << "Comparing current run to solution at: "
                                    << solution_filename << endl;
            solution_file.open(solution_filename);
            true_solution.Load(solution_file, X.Size());
            solution_file.close();
            Vector residual(X.Size());
            subtract(X, true_solution, residual);

            const double abs_error = sqrt(InnerProduct(MPI_COMM_WORLD, residual, residual));
            const double sol_norm = sqrt(InnerProduct(MPI_COMM_WORLD, true_solution,
                                         true_solution));

            curr_error = abs_error / sol_norm;
            if (myid == 0) cout << "The relative error is: " << curr_error << endl;

            Vector AX(X.Size());
            A.Mult(X, AX);
            subtract(B, AX, residual);
            curr_error = sqrt(InnerProduct(MPI_COMM_WORLD, residual, residual));
            if (myid == 0) cout << "The error indicator is: "
                                    << curr_error << endl;
        }

        // 30. If calculating the relative error, or after we sampled our first point,
        //     create a global ROM basis.
        if (calc_rel_error || (offline && basisIdentifiers.size() == 1))
        {
            mergeTimer.Start();
            std::unique_ptr<CAROM::BasisGenerator> basis_generator;
            options = new CAROM::Options(fespace.GetTrueVSize(), max_num_snapshots, 1,
                                         update_right_SV);
            generator = new CAROM::BasisGenerator(*options, isIncremental, loadBasisName);
            for (int i = 0; i < basisIdentifiers.size(); ++i)
            {
                std::string snapshot_filename = basisIdentifiers[i] + "_snapshot";
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
        }

        // 31. print timing info
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

        // 32. Free the used memory.
        if (delete_fec)
        {
            delete fec;
        }
    } while(build_database);
    MPI_Finalize();

    return 0;
}

// 33. define spatially varying righthand side function
double rhs(const Vector &x)
{
    if (dim == 3)
    {
        return sin(kappa * (x(0) + x(1) + x(2)));
    }
    else
    {
        return sin(kappa * (x(0) + x(1)));
    }

}
