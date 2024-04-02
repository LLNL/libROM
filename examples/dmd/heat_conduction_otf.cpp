/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

//                       libROM MFEM Example: Heat_Conduction (adapted from ex16p.cpp)
//
// Compile with: make heat_conduction
//
// =================================================================================
//
// Sample runs and results for DMD:
//
// Command 1:
//   mpirun -np 8 heat_conduction -s 1 -a 0.0 -k 1.0 -visit
//
// Output 1:
//   Relative error of DMD temperature (u) at t_final: 0.5 is 0.00049906934
//
// Command 2:
//   mpirun -np 8 heat_conduction -s 3 -a 0.5 -k 0.5 -o 4 -tf 0.7 -vs 1 -visit
//
// Output 2:
//   Relative error of DMD temperature (u) at t_final: 0.7 is 0.00082289823
//
// =================================================================================
//
// Sample runs and results for NonuniformDMD:
//
// Command 1:
//   mpirun -np 8 heat_conduction -s 1 -a 0.0 -k 1.0 -nonunif -visit
//
// Output 1:
//   Relative error of DMD temperature (u) at t_final: 0.5 is 0.00071958947
//
// Command 2:
//   mpirun -np 8 heat_conduction -s 3 -a 0.5 -k 0.5 -o 4 -tf 0.7 -vs 1 -nonunif -visit
//
// Output 2:
//   Relative error of DMD temperature (u) at t_final: 0.7 is 0.00013450754
//
// =================================================================================
//
// Description:  This example solves a time dependent nonlinear heat equation
//               problem of the form du/dt = C(u), with a non-linear diffusion
//               operator C(u) = \nabla \cdot (\kappa + \alpha u) \nabla u.
//
//               The example demonstrates the use of nonlinear operators (the
//               class ConductionOperator defining C(u)), as well as their
//               implicit time integration. Note that implementing the method
//               ConductionOperator::ImplicitSolve is the only requirement for
//               high-order implicit (SDIRK) time integration. Optional saving
//               with ADIOS2 (adios2.readthedocs.io) is also illustrated.

#include "mfem.hpp"
#include "algo/DMD.h"
#include "algo/NonuniformDMD.h"
#include "algo/IncrementalDMD.h"
#include "linalg/Vector.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

/** After spatial discretization, the conduction model can be written as:
 *
 *     du/dt = M^{-1}(-Ku)
 *
 *  where u is the vector representing the temperature, M is the mass matrix,
 *  and K is the diffusion operator with diffusivity depending on u:
 *  (\kappa + \alpha u).
 *
 *  Class ConductionOperator represents the right-hand side of the above ODE.
 */
class ConductionOperator : public TimeDependentOperator
{
protected:
    ParFiniteElementSpace &fespace;
    Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.

    ParBilinearForm *M;
    ParBilinearForm *K;

    HypreParMatrix Mmat;
    HypreParMatrix Kmat;
    HypreParMatrix *T; // T = M + dt K
    double current_dt;

    CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
    HypreSmoother M_prec; // Preconditioner for the mass matrix M

    CGSolver T_solver;    // Implicit solver for T = M + dt K
    HypreSmoother T_prec; // Preconditioner for the implicit solver

    double alpha, kappa;

    mutable Vector z; // auxiliary vector

public:
    ConductionOperator(ParFiniteElementSpace &f, double alpha, double kappa,
                       const Vector &u, double tol_cg);

    virtual void Mult(const Vector &u, Vector &du_dt) const;
    /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
        This is the only requirement for high-order SDIRK implicit integration.*/
    virtual void ImplicitSolve(const double dt, const Vector &u, Vector &k);

    /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
    void SetParameters(const Vector &u);
    
    Vector* Residual(const double dt, const Vector u0, const Vector u1);

    virtual ~ConductionOperator();
};

double InitialTemperature(const Vector &x);

int main(int argc, char *argv[])
{
    // 1. Initialize MPI.
    int num_procs, myid;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // 2. Parse command-line options.
    const char *mesh_file = "../data/star.mesh";
    int ser_ref_levels = 2;
    int par_ref_levels = 1;
    int order = 2;
    int ode_solver_type = 1;
    double t_final = 5.0;
    double dt = 1.0e-2;
    double alpha = 1.0e-2;
    double kappa = 0.5;
    double ef = 0.9999;
    int rdim = 1;
    bool visualization = true;
    bool visit = false;
    int vis_steps = 5;
    bool use_nonuniform = false;
    bool adios2 = false;

    double tol_cg = 1e-12; // tolerance for CG solver
    double tol_dmd = 1e-2; // tolerance for DMD accuracy
    double tol_dep = 1e-6; // tolerance for linear dependence
    
    bool use_incremental = true; // use incremental SVD

    int precision = 8;
    cout.precision(precision);

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                   "Number of times to refine the mesh uniformly in serial.");
    args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                   "Number of times to refine the mesh uniformly in parallel.");
    args.AddOption(&order, "-o", "--order",
                   "Order (degree) of the finite elements.");
    args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                   "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3,\n\t"
                   "\t   11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.");
    args.AddOption(&t_final, "-tf", "--t-final",
                   "Final time; start time is 0.");
    args.AddOption(&dt, "-dt", "--time-step",
                   "Time step.");
    args.AddOption(&alpha, "-a", "--alpha",
                   "Alpha coefficient.");
    args.AddOption(&kappa, "-k", "--kappa",
                   "Kappa coefficient offset.");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                   "--no-visualization",
                   "Enable or disable GLVis visualization.");
    args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                   "--no-visit-datafiles",
                   "Save data files for VisIt (visit.llnl.gov) visualization.");
    args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                   "Visualize every n-th timestep.");
    args.AddOption(&adios2, "-adios2", "--adios2-streams", "-no-adios2",
                   "--no-adios2-streams",
                   "Save data using adios2 streams.");
    args.AddOption(&ef, "-ef", "--energy_fraction",
                   "Energy fraction for DMD.");
    args.AddOption(&rdim, "-rdim", "--rdim",
                   "Reduced dimension for DMD.");
    args.AddOption(&use_nonuniform, "-nonunif", "--nonunif", "-no-nonunif",
                   "--no-nonunif",
                   "Use NonuniformDMD.");
    
    args.AddOption(&tol_cg, "-tolcg", "--tolcg", "CG tolerance.");
    args.AddOption(&tol_dmd, "-toldmd", "--toldmd", "DMD accuracy.");
    args.AddOption(&tol_dep, "-toldep", "--toldep", "Linear dependence of snapshots.");
    
    args.AddOption(&use_incremental, "-inc", "--inc", 
		    		     "-noinc", "--noinc",
				     "Incremental SVD.");
    
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

    // 3. Read the serial mesh from the given mesh file on all processors. We can
    //    handle triangular, quadrilateral, tetrahedral and hexahedral meshes
    //    with the same code.
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    // 4. Define the ODE solver used for time integration. Several implicit
    //    singly diagonal implicit Runge-Kutta (SDIRK) methods, as well as
    //    explicit Runge-Kutta methods are available.
    ODESolver *ode_solver;
    switch (ode_solver_type)
    {
    // Implicit L-stable methods
    case 1:
        ode_solver = new BackwardEulerSolver;
        break;
    case 2:
        ode_solver = new SDIRK23Solver(2);
        break;
    case 3:
        ode_solver = new SDIRK33Solver;
        break;
    // Explicit methods
    case 11:
        ode_solver = new ForwardEulerSolver;
        break;
    case 12:
        ode_solver = new RK2Solver(0.5);
        break; // midpoint method
    case 13:
        ode_solver = new RK3SSPSolver;
        break;
    case 14:
        ode_solver = new RK4Solver;
        break;
    case 15:
        ode_solver = new GeneralizedAlphaSolver(0.5);
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
        cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
        delete mesh;
        return 3;
    }

    // 5. Refine the mesh in serial to increase the resolution. In this example
    //    we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
    //    a command-line parameter.
    for (int lev = 0; lev < ser_ref_levels; lev++)
    {
        mesh->UniformRefinement();
    }

    // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
    //    this mesh further in parallel to increase the resolution. Once the
    //    parallel mesh is defined, the serial mesh can be deleted.
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
    for (int lev = 0; lev < par_ref_levels; lev++)
    {
        pmesh->UniformRefinement();
    }

    // 7. Define the vector finite element space representing the current and the
    //    initial temperature, u_ref.
    H1_FECollection fe_coll(order, dim);
    ParFiniteElementSpace fespace(pmesh, &fe_coll);

    int fe_size = fespace.GlobalTrueVSize();
    if (myid == 0)
    {
        cout << "Number of temperature unknowns: " << fe_size << endl;
    }

    ParGridFunction u_gf(&fespace);

    // 8. Set the initial conditions for u. All boundaries are considered
    //    natural.
    FunctionCoefficient u_0(InitialTemperature);
    u_gf.ProjectCoefficient(u_0);
    Vector u;
    u_gf.GetTrueDofs(u);

    // 9. Initialize the conduction operator and the VisIt visualization.
    ConductionOperator oper(fespace, alpha, kappa, u, tol_cg);

    VisItDataCollection visit_dc("Heat_Conduction", pmesh);
    visit_dc.RegisterField("temperature", &u_gf);
    if (visit)
    {
        visit_dc.SetCycle(0);
        visit_dc.SetTime(0.0);
        visit_dc.Save();
    }

    // Optionally output a BP (binary pack) file using ADIOS2. This can be
    // visualized with the ParaView VTX reader.
#ifdef MFEM_USE_ADIOS2
    ADIOS2DataCollection* adios2_dc = NULL;
    if (adios2)
    {
        std::string postfix(mesh_file);
        postfix.erase(0, std::string("../data/").size() );
        postfix += "_o" + std::to_string(order);
        postfix += "_solver" + std::to_string(ode_solver_type);
        const std::string collection_name = "heat_conduction-p-" + postfix + ".bp";

        adios2_dc = new ADIOS2DataCollection(MPI_COMM_WORLD, collection_name, pmesh);
        adios2_dc->SetParameter("SubStreams", std::to_string(num_procs/2) );
        adios2_dc->RegisterField("temperature", &u_gf);
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
            sout.precision(precision);
            sout << "solution\n" << *pmesh << u_gf;
            sout << flush;
        }
    }

    StopWatch fom_timer, dmd_training_timer, dmd_prediction_timer;

    fom_timer.Start();

    // 10. Perform time-integration (looping over the time iterations, ti, with a
    //     time-step dt).
    ode_solver->Init(oper);
    double t = 0.0;
    vector<double> ts;

    fom_timer.Stop();

    dmd_training_timer.Start();

    // 11. Create DMD object and take initial sample.
    u_gf.SetFromTrueDofs(u);

    CAROM::DMD* dmd_u;
    if (use_incremental) {
	   CAROM::Options svd_options(u.Size(), 100, -1, true);
	   svd_options.setIncrementalSVD(1e-12, dt, 1e-6, 10.0, false, true, false);
	   svd_options.setMaxBasisDimension(100);
	   svd_options.setStateIO(false, false);
	   svd_options.setDebugMode(false);
	   std::string svd_base_file_name = "";
	   dmd_u = new CAROM::IncrementalDMD(u.Size(), dt,
					  svd_options,
					  svd_base_file_name,
					  false);
    }
    else {
    	dmd_u = new CAROM::DMD(u.Size(), dt);
    }

    // Push the initial snapshot.
    dmd_u->takeSample(u.GetData(), t);
    ts.push_back(t);
  
    // 13. Time intergration with on-the-fly DMD construction.
    if (myid == 0 && rdim != -1 && ef != -1)
    {
        std::cout << "Both rdim and ef are set. ef will be ignored." << std::endl;
    }

    bool dmd_trained = false;
    bool last_step = false;
    double res, err_proj_norm;
    Vector u_dmd(u.Size()); // DMD prediction
    Vector* res_vec = NULL;
    
    for (int ti = 1; !last_step; ti++)
    {
	StopWatch total_timer;
	total_timer.Start();

    if (myid == 0) std::cout << "\nIteration " << ti << std::endl;
	
	if (t + dt >= t_final - dt/2)
	{
	    last_step = true;
	}

	// ROM solve
	if (dmd_trained)
	{
	    if (myid == 0){
		   std::cout << "DMD exists: make prediction" << std::endl;
	    }
	    
	    CAROM::Vector* u_pre = new CAROM::Vector(u.GetData(), u.Size(),
			    			     true, true); // previous solution
    	
        if (use_incremental) {
           // Incremental DMD
	       u_dmd = dmd_u->predict_dt(u_pre)->getData();
	    }
	    else {
           // Vanilla DMD
		   dmd_u->projectInitialCondition(u_pre); // offset to u_pre
		   CAROM::Vector* u_dmd_pred = dmd_u->predict(dt);
		   u_dmd = u_dmd_pred->getData(); // copy data to MFEM::Vector
		   delete u_dmd_pred;
	    }

	    res_vec = oper.Residual(dt, u, u_dmd);
	    res = res_vec->Norml2()*res_vec->Norml2(); // Residual (L2): not distributed
	    delete res_vec;
	    MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	    if (myid == 0){
		   std::cout << "Residual of DMD solution is: " << res << std::endl;
	    }

	    if (res < tol_dmd) {
           // DMD prediction is accurate
		   if (myid == 0){
		      std::cout << "DMD prediction accurate: use it" << std::endl;
		   }
		   u = u_dmd; // use DMD solution as new snapshot
		   t += dt;
	    }
	    else {
           // FOM solve
		   if (myid == 0){
		      std::cout << "DMD prediction not accurate: call FOM" << std::endl;
		   }

		   fom_timer.Start();
		   ode_solver->Step(u, t, dt);
		   fom_timer.Stop();
	    }
	    
	    if (!use_incremental){ 
	       // If using vanilla DMD,
           // check the linear dependence of the new snapshot
	       CAROM::Vector* u_new = new CAROM::Vector(u.GetData(), u.Size(),
				    			                     true, false);
	       dmd_u->projectInitialCondition(u_pre);
		
		   CAROM::Vector* u_dmd_pred_0 = dmd_u->predict(0.0);
	       Vector u_new_proj(u_dmd_pred_0->getData(), u.Size()); // projection
	       Vector err_proj(u_pre->getData(), u.Size());
	       err_proj -= u_new_proj; // error = u - proj(u)
	       err_proj_norm = err_proj.Norml2();// / u.Norml2(); // relative norm
	       std::cout << "Projection error: " << err_proj_norm << std::endl;
		   delete u_dmd_pred_0;
	    	
	       // Increment no. of modes if solution is inaccurate and
	       // the new snapshot is linearly independent
	       if (res > tol_dmd && err_proj_norm > tol_dep){
		      rdim += 1;
	       }
		   delete u_new;
	    }

	    delete u_pre;

	}
	else // No constructed DMD: FOM solve
	{
	    if (myid == 0){
		   std::cout << "No DMD found: solve FOM" << std::endl;
	    }

	    fom_timer.Start();
	    ode_solver->Step(u, t, dt);
	    fom_timer.Stop();
	}

	// Append new snapshot
	dmd_u->takeSample(u.GetData(), t);
	ts.push_back(t);
	
	// DMD training
    dmd_training_timer.Start();
    
	if (dmd_u->getNumSamples() > rdim)
	{
        if (myid == 0) {
           std::cout << "Creating DMD with rdim: " << rdim << std::endl;
        }
        
        StopWatch train_timer;
	    train_timer.Start();
	    dmd_u->train(rdim);
	    train_timer.Stop();
	    
        if (myid == 0) std::cout << "Time train:" << train_timer.RealTime() << std::endl;
	    
        dmd_trained = true;

	}
	

	// Visualize solutions
        u_gf.SetFromTrueDofs(u);
        
	if (last_step || (ti % vis_steps) == 0)
        {
            if (myid == 0)
            {
                cout << "step " << ti << ", t = " << t << endl;
            }

            u_gf.SetFromTrueDofs(u);
            if (visualization)
            {
                sout << "parallel " << num_procs << " " << myid << "\n";
                sout << "solution\n" << *pmesh << u_gf << flush;
            }

            if (visit)
            {
                visit_dc.SetCycle(ti);
                visit_dc.SetTime(t);
                visit_dc.Save();
            }

#ifdef MFEM_USE_ADIOS2
            if (adios2)
            {
                adios2_dc->SetCycle(ti);
                adios2_dc->SetTime(t);
                adios2_dc->Save();
            }
#endif
        }
        oper.SetParameters(u);

	total_timer.Stop();
	std::cout << "Time total:" << total_timer.RealTime() << std::endl;
    }

#ifdef MFEM_USE_ADIOS2
    if (adios2)
    {
        delete adios2_dc;
    }
#endif

    dmd_training_timer.Stop();
    if (myid == 0) {
    	std::cout << "Total training time:" << dmd_training_timer.RealTime() << std::endl;
    }

    // 12. Save the final solution in parallel. This output can be viewed later
    //     using GLVis: "glvis -np <np> -m heat_conduction-mesh -g heat_conduction-final".
    {
        ostringstream sol_name;
        sol_name << "heat_conduction-final." << setfill('0') << setw(6) << myid;
        ofstream osol(sol_name.str().c_str());
        osol.precision(precision);
        u_gf.Save(osol);
    }

    // 13. Free the used memory.
    delete ode_solver;
    delete pmesh;

    MPI_Finalize();

    return 0;
}

ConductionOperator::ConductionOperator(ParFiniteElementSpace &f, double al,
                                       double kap, const Vector &u, double tol_cg)
    : TimeDependentOperator(f.GetTrueVSize(), 0.0), fespace(f), M(NULL), K(NULL),
      T(NULL), current_dt(0.0),
      M_solver(f.GetComm()), T_solver(f.GetComm()), z(height)
{
    const double rel_tol = tol_cg;

    M = new ParBilinearForm(&fespace);
    M->AddDomainIntegrator(new MassIntegrator());
    M->Assemble(0); // keep sparsity pattern of M and K the same
    M->FormSystemMatrix(ess_tdof_list, Mmat);

    M_solver.iterative_mode = false;
    M_solver.SetRelTol(rel_tol);
    M_solver.SetAbsTol(0.0);
    M_solver.SetMaxIter(1000);
    M_solver.SetPrintLevel(0);
    M_prec.SetType(HypreSmoother::Jacobi);
    M_solver.SetPreconditioner(M_prec);
    M_solver.SetOperator(Mmat);

    alpha = al;
    kappa = kap;

    T_solver.iterative_mode = false;
    T_solver.SetRelTol(rel_tol);
    T_solver.SetAbsTol(0.0);
    T_solver.SetMaxIter(1000);
    T_solver.SetPrintLevel(0);
    T_solver.SetPreconditioner(T_prec);

    SetParameters(u);
}

void ConductionOperator::Mult(const Vector &u, Vector &du_dt) const
{
    // Compute:
    //    du_dt = M^{-1}*-K(u)
    // for du_dt
    Kmat.Mult(u, z);
    z.Neg(); // z = -z
    M_solver.Mult(z, du_dt);
}

void ConductionOperator::ImplicitSolve(const double dt,
                                       const Vector &u, Vector &du_dt)
{
    // Solve the equation:
    //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
    // for du_dt
    if (!T)
    {
        T = Add(1.0, Mmat, dt, Kmat);
        current_dt = dt;
        T_solver.SetOperator(*T);
    }
    MFEM_VERIFY(dt == current_dt, ""); // SDIRK methods use the same dt
    Kmat.Mult(u, z);
    z.Neg();
    T_solver.Mult(z, du_dt);
}

void ConductionOperator::SetParameters(const Vector &u)
{
    ParGridFunction u_alpha_gf(&fespace);
    u_alpha_gf.SetFromTrueDofs(u);
    for (int i = 0; i < u_alpha_gf.Size(); i++)
    {
        u_alpha_gf(i) = kappa + alpha*u_alpha_gf(i);
    }

    delete K;
    K = new ParBilinearForm(&fespace);

    GridFunctionCoefficient u_coeff(&u_alpha_gf);

    K->AddDomainIntegrator(new DiffusionIntegrator(u_coeff));
    K->Assemble(0); // keep sparsity pattern of M and K the same
    K->FormSystemMatrix(ess_tdof_list, Kmat);
    delete T;
    T = NULL; // re-compute T on the next ImplicitSolve
}

Vector* ConductionOperator::Residual(const double dt,
				    const Vector u0, const Vector u1)
{
    //
    // Compute res(u0, u1) = [(M+dt*K(u0))*u1 - M*u0] / dt / norm(K*u0)
    //
    T = Add(1.0, Mmat, dt, Kmat); // T=M+dt*K
    Vector* res = new Vector(u1.Size());
    T->Mult(u1, *res); // (M+dt*K)*u^n
    Kmat.Mult(u0, z); // z=K*u^{n-1}
    
    double norm = z.Norml2()*z.Norml2(); // norm(dt*K*u^{n-1}): local
    MPI_Allreduce(MPI_IN_PLACE,
		  &norm,
		  1,
		  MPI_DOUBLE,
		  MPI_SUM,
		  MPI_COMM_WORLD);
    norm = sqrt(norm) * dt;
    Mmat.Mult(u0, z); // z=M*u^{n-1}
    delete T;
    T = NULL;
    
    *res -= z;
    *res /= norm;
    
    return res; // local
}

ConductionOperator::~ConductionOperator()
{
    delete T;
    delete M;
    delete K;
}

double InitialTemperature(const Vector &x)
{
    if (x.Norml2() < 0.5)
    {
        return 2.0;
    }
    else
    {
        return 1.0;
    }
}
