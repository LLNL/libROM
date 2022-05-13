//                                MFEM Example 16
//
// Compile with: make ex16
//
// Sample runs:  ex16
//               ex16 -m ../data/inline-tri.mesh
//               ex16 -m ../data/disc-nurbs.mesh -tf 2
//               ex16 -s 1 -a 0.0 -k 1.0
//               ex16 -s 2 -a 1.0 -k 0.0
//               ex16 -s 3 -a 0.5 -k 0.5 -o 4
//               ex16 -s 14 -dt 1.0e-4 -tf 4.0e-2 -vs 40
//               ex16 -m ../data/fichera-q2.mesh
//               ex16 -m ../data/fichera-mixed.mesh
//               ex16 -m ../data/escher.mesh
//               ex16 -m ../data/beam-tet.mesh -tf 10 -dt 0.1
//               ex16 -m ../data/amr-quad.mesh -o 4 -r 0
//               ex16 -m ../data/amr-hex.mesh -o 2 -r 0
//
// Description:  This example solves a time dependent nonlinear heat equation
//               problem of the form du/dt = C(u), with a non-linear diffusion
//               operator C(u) = \nabla \cdot (\kappa + \alpha u) \nabla u.
//
//               The example demonstrates the use of nonlinear operators (the
//               class ConductionOperator defining C(u)), as well as their
//               implicit time integration. Note that implementing the method
//               ConductionOperator::ImplicitSolve is the only requirement for
//               high-order implicit (SDIRK) time integration.
//
//               We recommend viewing examples 2, 9 and 10 before viewing this
//               example.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;


// parameters for initial condition
double h = 0.5; // height of initial condition
double w = 0.2; // width of initial condition
double x1 = 0.5, x2 = 0.5; // center coordinates of initial condition

// parameters related to residual calculation
bool res_flag = false; // flag to indicate whether residual is calculated from given solutions
double res_ns = 0.1; // percentage of the time steps for residual evaluation
int Tmax_iter = 100; // maximum CG iterations in Tsolver
const char *u_file = "./ex16-u_pred.gf"; // file path of given solutions
const char *res_file = "./ex16-residual.gf"; // file path to save calculated residual


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
   FiniteElementSpace &fespace;
   Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.

   BilinearForm *M;
   BilinearForm *K;

   SparseMatrix Mmat, Kmat;
   SparseMatrix *T; // T = M + dt K
   double current_dt;

   CGSolver M_solver; // Krylov solver for inverting the mass matrix M
   DSmoother M_prec;  // Preconditioner for the mass matrix M

   CGSolver T_solver; // Implicit solver for T = M + dt K
   DSmoother T_prec;  // Preconditioner for the implicit solver

   double alpha, kappa;

   mutable Vector z; // auxiliary vector

public:
   ConductionOperator(FiniteElementSpace &f, double alpha, double kappa,
                      const Vector &u);

   virtual void Mult(const Vector &u, Vector &du_dt) const;
   /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const double dt, const Vector &u, Vector &k);

   /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
   void SetParameters(const Vector &u);

   virtual ~ConductionOperator();
};

double InitialTemperature(const Vector &x);

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "../data/star.mesh";
   int ref_levels = 2;
   int order = 2;
   int ode_solver_type = 3;
   double t_final = 0.5;
   double dt = 1.0e-2;
   double alpha = 1.0e-2;
   double kappa = 0.5;
   bool visualization = false;
   bool visit = false;
   int vis_steps = 5;

   int precision = 18;
   cout.precision(precision);
    
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&u_file, "-uf", "--sol_file",
                  "Input file of predicted solution");
   args.AddOption(&res_file, "-rf", "--res_file",
                  "Output file of residual");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly.");
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
   args.AddOption(&h, "-h0", "--height",
                  "Height of initial condition");
   args.AddOption(&w, "-w0", "--width",
                  "Width of initial condition");
   args.AddOption(&x1, "-x1", "--x1",
                  "x coordinate of the center of initial condition");
   args.AddOption(&x2, "-x2", "--x2",
                  "y coordinate of the center of initial condition");
   args.AddOption(&res_ns, "-res_ns", "--residual_num_steps",
                  "Percentage of the time steps for residual evaluation");
   args.AddOption(&Tmax_iter, "-Tmax_iter", "--Tsolver_max_iter",
                  "Maximum CG iterations in Tsolver");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                  "--no-visit-datafiles",
                  "Save data files for VisIt (visit.llnl.gov) visualization.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
   args.AddOption(&res_flag, "-res", "--residual", "-no-res", "--no-residual",
                  "Calculate residual from given solutions");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
//    args.PrintOptions(cout);

   // 2. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral and hexahedral meshes with the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   // 3. Define the ODE solver used for time integration. Several implicit
   //    singly diagonal implicit Runge-Kutta (SDIRK) methods, as well as
   //    explicit Runge-Kutta methods are available.
   ODESolver *ode_solver;
   switch (ode_solver_type)
   {
      // Implicit L-stable methods
      case 1:  ode_solver = new BackwardEulerSolver; break;
      case 2:  ode_solver = new SDIRK23Solver(2); break;
      case 3:  ode_solver = new SDIRK33Solver; break;
      // Explicit methods
      case 11: ode_solver = new ForwardEulerSolver; break;
      case 12: ode_solver = new RK2Solver(0.5); break; // midpoint method
      case 13: ode_solver = new RK3SSPSolver; break;
      case 14: ode_solver = new RK4Solver; break;
      case 15: ode_solver = new GeneralizedAlphaSolver(0.5); break;
      // Implicit A-stable methods (not L-stable)
      case 22: ode_solver = new ImplicitMidpointSolver; break;
      case 23: ode_solver = new SDIRK23Solver; break;
      case 24: ode_solver = new SDIRK34Solver; break;
      default:
         cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
         delete mesh;
         return 3;
   }

   // 4. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement, where 'ref_levels' is a
   //    command-line parameter.
   for (int lev = 0; lev < ref_levels; lev++)
   {
      mesh->UniformRefinement();
   }

   // 5. Define the vector finite element space representing the current and the
   //    initial temperature, u_ref.
   H1_FECollection fe_coll(order, dim);
   FiniteElementSpace fespace(mesh, &fe_coll);

   int fe_size = fespace.GetTrueVSize();
//    cout << "Number of temperature unknowns: " << fe_size << endl;

   GridFunction u_gf(&fespace);

   // 6. Set the initial conditions for u. All boundaries are considered
   //    natural.
   FunctionCoefficient u_0(InitialTemperature);
   u_gf.ProjectCoefficient(u_0);
   Vector u;
   u_gf.GetTrueDofs(u);

   // 7. Initialize the conduction operator and the visualization.
   ConductionOperator oper(fespace, alpha, kappa, u);

   u_gf.SetFromTrueDofs(u);
    
   if (!res_flag)
   {
      ofstream omesh("ex16.mesh");
      omesh.precision(precision);
      mesh->Print(omesh);
      ofstream osol("ex16-init.gf");
      osol.precision(precision);
      u_gf.Save(osol);
   }

   VisItDataCollection visit_dc("Example16", mesh);
   visit_dc.RegisterField("temperature", &u_gf);
   if (visit)
   {
      visit_dc.SetCycle(0);
      visit_dc.SetTime(0.0);
      visit_dc.Save();
   }

   
   socketstream sout;
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      sout.open(vishost, visport);
      if (!sout)
      {
         cout << "Unable to connect to GLVis server at "
              << vishost << ':' << visport << endl;
         visualization = false;
         cout << "GLVis visualization disabled.\n";
      }
      else
      {
         sout.precision(precision);
         sout << "solution\n" << *mesh << u_gf;
         sout << "pause\n";
         sout << flush;
         cout << "GLVis visualization paused."
              << " Press space (in the GLVis window) to resume it.\n";
      }
   }

   ode_solver->Init(oper);
    
    
   if (!res_flag)
   {
    // output the initial condition of u and dudt
    GridFunction tmp_u(&fespace), tmp_dudt(&fespace);
    tmp_u.SetFromTrueDofs(u);
    ofstream osol("ex16-u_"+to_string(0)+".gf");
    osol.precision(precision);
    tmp_u.Save(osol);
    osol.close();

    oper.Mult(u, tmp_dudt);
    osol.open("ex16-du_"+to_string(0)+".gf");
    osol.precision(precision);
    tmp_dudt.Save(osol);
    osol.close();
       
       
   // 8. Perform time-integration (looping over the time iterations, ti, with a
   //    time-step dt).
   double t = 0.0;
   bool last_step = false;

   for (int ti = 1; !last_step; ti++)
   {
      if (t + dt >= t_final - dt/2)
      {
         last_step = true;
      }

      ode_solver->Step(u, t, dt);
  
      if (last_step || (ti % vis_steps) == 0)
      {
         cout << "step " << ti << ", t = " << t << endl;

         u_gf.SetFromTrueDofs(u);
         if (visualization)
         {
            sout << "solution\n" << *mesh << u_gf << flush;
         }

         if (visit)
         {
            visit_dc.SetCycle(ti);
            visit_dc.SetTime(t);
            visit_dc.Save();
         }
      }
      oper.SetParameters(u);
       
      // output u and dudt at every time step
      GridFunction tmp_u(&fespace), tmp_dudt(&fespace);
      tmp_u.SetFromTrueDofs(u);
      ofstream osol("ex16-u_"+to_string(ti)+".gf");
      osol.precision(precision);
      tmp_u.Save(osol);
      osol.close();

      oper.Mult(u, tmp_dudt);
      osol.open("ex16-du_"+to_string(ti)+".gf");
      osol.precision(precision);
      tmp_dudt.Save(osol);
      osol.close();

   }

   // 9. Save the final solution. This output can be viewed later using GLVis:
   //    "glvis -m ex16.mesh -g ex16-final.gf".
   {
      ofstream osol("ex16-final_u.gf");
      osol.precision(precision);
      u_gf.Save(osol);
      osol.close();
       
      GridFunction dudt(&fespace);
      oper.Mult(u, dudt);
      osol.open("ex16-final_dudt.gf");
      osol.precision(precision);
      dudt.Save(osol);
      osol.close();
   }
   }
   else // calculate residual from given solution
   {
   ifstream infile(u_file);
   infile.precision(precision); 
   int d1, d2;
   infile >> d1;
   infile >> d2;
//    cout << d1 << ' ' << d2 << endl;

   GridFunction res(&fespace), u1(&fespace), u2(&fespace);
   double t = 0, res_norm = 0.0, res_norm_tmp = 0.0;
   u1.Load(infile, d2); // u at t=0
   for (int i = 0; i < int(d1*res_ns)-1; i++){
       oper.SetParameters(u1);
       ode_solver->Step(u1, t, dt); // true next u
       u2.Load(infile, d2); // predicted next step
       subtract(u1, u2, res);
       res_norm_tmp = res.Norml2();
       res_norm += res_norm_tmp;
       cout << i+1 << " - " << "t: " << t << "  residual: " << res_norm_tmp << endl;
       u1 = u2;
   }
   cout << "final residual: " << res_norm << endl;
       
   // output the residual norm
   ofstream outfile(res_file);
   outfile.precision(precision);
   outfile << res_norm << endl;
//    cout << "output done" << endl;
   }
 
   // 10. Free the used memory.
   delete ode_solver;
   delete mesh;

   return 0;
}

ConductionOperator::ConductionOperator(FiniteElementSpace &f, double al,
                                       double kap, const Vector &u)
   : TimeDependentOperator(f.GetTrueVSize(), 0.0), fespace(f), M(NULL), K(NULL),
     T(NULL), current_dt(0.0), z(height)
{
   const double rel_tol = 1e-8;

   M = new BilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator());
   M->Assemble();
   M->FormSystemMatrix(ess_tdof_list, Mmat);

   M_solver.iterative_mode = false;
   M_solver.SetRelTol(rel_tol);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(30);
   M_solver.SetPrintLevel(0);
   M_solver.SetPreconditioner(M_prec);
   M_solver.SetOperator(Mmat);

   alpha = al;
   kappa = kap;

   T_solver.iterative_mode = false;
   T_solver.SetRelTol(rel_tol);
   T_solver.SetAbsTol(0.0);
   T_solver.SetMaxIter(Tmax_iter);
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
   GridFunction u_alpha_gf(&fespace);
   u_alpha_gf.SetFromTrueDofs(u);
   for (int i = 0; i < u_alpha_gf.Size(); i++)
   {
      u_alpha_gf(i) = kappa + alpha*u_alpha_gf(i);
   }

   delete K;
   K = new BilinearForm(&fespace);

   GridFunctionCoefficient u_coeff(&u_alpha_gf);

   K->AddDomainIntegrator(new DiffusionIntegrator(u_coeff));
   K->Assemble();
   K->FormSystemMatrix(ess_tdof_list, Kmat);
   delete T;
   T = NULL; // re-compute T on the next ImplicitSolve
}

ConductionOperator::~ConductionOperator()
{
   delete T;
   delete M;
   delete K;
}

double InitialTemperature(const Vector &x)
{
//     double rd = sqrt(pow(x(0)-x1,2) + pow(x(1)-x2,2));
//    if (rd < w)
//    {
//       return h;
//    }
//    else
//    {
//       return 0.0;
//    }
    
//     return h * sin(w*(x(0)+x(1))) + h;
    return h * sin(w*x.Norml2()) + h;
//     return h * sin(w*x(0)) + h;
//     return h * (sin(w*x(0)) + cos(w*x(1))) + h;
}
