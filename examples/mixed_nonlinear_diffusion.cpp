//               libROM MFEM Example: parametric ROM for nonlinear diffusion problem
//
// Compile with: ./scripts/compile.sh -m 
//
// Description:  This example solves a time dependent nonlinear diffusion equation
//               du/dt + div(v) = f, grad u = -a(u) v. After discretization by mixed FEM,
//               M(u) v + B^T u = 0
//               B v - C u_t = -f
//               where 
//               M(u) = \int_\Omega a(u) w_h \cdot v_h d\Omega   w_h, v_h \in R_h
//               B = -\int_\Omega \div w_h q_h d\Omega   w_h \in R_h, q_h \in W_h
//               C = \int_\Omega q_h p_h d\Omega   p_h \in W_h, q_h \in W_h
//               Here, R_h is a Raviart-Thomas finite element subspace of H(div),
//               and W_h is a finite element subspace of L2.
//               The first equation allows the substitution v = -M(u)^{-1} B^T u, so
//               C u_t + B M(u)^{-1} B^T u = f
//               For the purpose of using an ODE solver, this can be expressed as
//               u_t = C^{-1} (f - B M(u)^{-1} B^T u) = F(u)

// Sample runs:
//               Analytic test
//               mpirun -n 1 ./mixed_nonlinear_diffusion -m ../../mfem/data/inline-quad.mesh -offline
//               mpirun -n 1 ./mixed_nonlinear_diffusion -m ../../mfem/data/inline-quad.mesh -merge -ns 1
//               mpirun -n 1 ./mixed_nonlinear_diffusion -m ../../mfem/data/inline-quad.mesh -online -rrdim 8 -rwdim 8
//
//               Initial step test
//               mpirun -n 1 ./mixed_nonlinear_diffusion -m ../../mfem/data/inline-quad.mesh -offline -p 1
//               mpirun -n 1 ./mixed_nonlinear_diffusion -m ../../mfem/data/inline-quad.mesh -merge -ns 1 -p 1
//               mpirun -n 1 ./mixed_nonlinear_diffusion -m ../../mfem/data/inline-quad.mesh -online -rrdim 8 -rwdim 8 -p 1
//
//               Initial step parametric test
//               mpirun -n 1 ./mixed_nonlinear_diffusion -m ../../mfem/data/inline-quad.mesh -offline -p 1 -id 0 -sh 0.25
//               mpirun -n 1 ./mixed_nonlinear_diffusion -m ../../mfem/data/inline-quad.mesh -offline -p 1 -id 1 -sh 0.15
//               mpirun -n 1 ./mixed_nonlinear_diffusion -m ../../mfem/data/inline-quad.mesh -offline -p 1 -id 2 -sh 0.35
//               mpirun -n 1 ./mixed_nonlinear_diffusion -m ../../mfem/data/inline-quad.mesh -merge -ns 3 -p 1
//               mpirun -n 1 ./mixed_nonlinear_diffusion -m ../../mfem/data/inline-quad.mesh -offline -p 1 -id 3 -sh 0.3
//               mpirun -n 1 ./mixed_nonlinear_diffusion -m ../../mfem/data/inline-quad.mesh -online -rrdim 8 -rwdim 8 -p 1 -sh 0.3 -id 3

#include "mfem.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>

#include "BasisGenerator.h"
#include "BasisReader.h"
#include "DEIM.h"
#include "SampleMesh.hpp"


typedef enum {ANALYTIC, INIT_STEP} PROBLEM;

//#define USE_GNAT

using namespace mfem;
using namespace std;

static bool nonlinear_problem;
static int problem;
static double diffusion_c, step_half;

class NonlinearDiffusionOperator;
class RomOperator;

class NonlinearDiffusionGradientOperator : public Operator
{
private:
  friend class NonlinearDiffusionOperator;

  void SetParameters(Operator *C_solver_, Operator *M_solver_, 
		     Operator *M_prime_, Operator *B_, const double dt_);

  Operator *C_solver;
  Operator *M_solver;
  Operator *M_prime;
  Operator *B;

  double dt;

  mutable Vector zR; // auxiliary vector
  mutable Vector yR; // auxiliary vector
  mutable Vector xR; // auxiliary vector
  mutable Vector zW; // auxiliary vector
  
public:
  NonlinearDiffusionGradientOperator(const int sizeR, const int sizeW);

  void Mult(const Vector &x, Vector &y) const;
};

NonlinearDiffusionGradientOperator::NonlinearDiffusionGradientOperator(const int sizeR,
								       const int sizeW)
  : Operator(sizeW), zW(sizeW), zR(sizeR), yR(sizeR), xR(sizeR),
    C_solver(NULL), M_solver(NULL), M_prime(NULL), B(NULL), dt(0.0)
{
}

void NonlinearDiffusionGradientOperator::SetParameters(Operator *C_solver_, Operator *M_solver_, 
						       Operator *M_prime_, Operator *B_, const double dt_)
{
  C_solver = C_solver_;
  M_solver = M_solver_;
  M_prime = M_prime_;
  B = B_;

  dt = dt_;
}

void NonlinearDiffusionGradientOperator::Mult(const Vector &x, Vector &y) const
{
  // Gradient is I + dt C^{-1} B M(u)^{-1} B^T - dt C^{-1} B M(u)^{-1} M(a'(u)) M(u)^{-1} B^T
  //   = I + dt C^{-1} B (I - M(u)^{-1} M(a'(u))) M(u)^{-1} B^T

  // Apply M(u)^{-1} B^T
  B->MultTranspose(x, zR);
  M_solver->Mult(zR, yR);

  // Apply (I - M(u)^{-1} M(a'(u))) to yR

  M_prime->Mult(yR, zR);
  M_solver->Mult(zR, xR);

  yR.Add(-1.0, xR);

  // Apply C^{-1} B
  B->Mult(yR, zW);
  C_solver->Mult(zW, y);

  zW = y;
  y = x;
  y.Add(dt, zW);
}

class NonlinearDiffusionOperator : public TimeDependentOperator
{
  
private:
  bool SchurComplement;

protected:
  friend class RomOperator;
  
  Array<int> ess_Rdof_list; // this list remains empty for pure essential b.c.

  mutable ParBilinearForm *M;
  mutable HypreParMatrix *Cmat, *Mmat;
  mutable ParBilinearForm *Mprime;
  ParBilinearForm *C;

  HypreParMatrix *Bmat;
  HypreParMatrix *BTmat;
  mutable HypreParMatrix Mprimemat;

  mutable BlockOperator *fullOp;
  mutable BlockOperator *fullGradient;
  mutable BlockDiagonalPreconditioner *fullPrec;

  Array<int> block_trueOffsets;

  double current_dt;

  mutable CGSolver *M_solver;    // Krylov solver for inverting the R mass matrix M
  mutable HypreSmoother M_prec;  // Preconditioner for the R mass matrix M

  mutable CGSolver C_solver;    // Krylov solver for inverting the W mass matrix C
  HypreSmoother C_prec; // Preconditioner for the W mass matrix C

  GMRESSolver *J_gmres;

  NewtonSolver newton_solver;

  NonlinearDiffusionGradientOperator *gradient;

  double linear_solver_rel_tol;

  Vector u0;
  Vector dudt_prev;

  mutable Vector zR; // auxiliary vector
  mutable Vector yR; // auxiliary vector
  mutable Vector zW; // auxiliary vector

public:
  NonlinearDiffusionOperator(ParFiniteElementSpace &fR, ParFiniteElementSpace &fW,
			     const double rel_tol, const double abs_tol,
			     const int iter, const Vector &u, const bool SchurComplement_);

  virtual void Mult(const Vector &u, Vector &du_dt) const;

  void Mult_Mmat(const Vector &u, Vector &Mu) const
  {
    Mmat->Mult(u, Mu);
  }
  
  void Mult_FullSystem(const Vector &u, Vector &du_dt) const;
  void Mult_SchurComplement(const Vector &u, Vector &du_dt) const;
  void SetBTV(const CAROM::Matrix *V, CAROM::Matrix *BTV) const;
  
  void GetSource(Vector& s) const;
  
  /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
      This is the only requirement for high-order SDIRK implicit integration.*/
  virtual void ImplicitSolve(const double dt, const Vector &u, Vector &du_dt);

  virtual Operator &GetGradient(const Vector &u) const;

  /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
  void SetParameters(const Vector &u) const;

  void CopyDuDt(Vector &dudt) const
  {
    dudt = dudt_prev;
  }

  void CopyDuDt_W(Vector &dudt) const
  {
    Vector dudt_W(dudt_prev.GetData() + zR.Size(), zW.Size());
    
    dudt = dudt_W;
  }

  virtual ~NonlinearDiffusionOperator();

  ParFiniteElementSpace &fespace_R;  // RT finite element space
  ParFiniteElementSpace &fespace_W;  // L2 discontinuous scalar finite element space

  bool newtonFailure;
};

class RomOperator : public TimeDependentOperator
{
private:
  int rrdim, rwdim, nldim;
  int nsamp_R, nsamp_S;
  double current_dt;
  NewtonSolver newton_solver;
  GMRESSolver *J_gmres;
  CAROM::Matrix *BRsp, *BWsp;
  CAROM::Vector *usp_librom, *usp_R_librom, *usp_W_librom;
  Vector *usp;
  Vector *usp_R;
  Vector *usp_W;
  mutable Vector zR;
  mutable CAROM::Vector zY;
  mutable CAROM::Vector zN;
  vector<int> s2sp;
  const CAROM::Matrix *Vsinv;

  // Data for source function
  vector<int> s2sp_S;
  const CAROM::Matrix *Ssinv;
  mutable CAROM::Vector zS;
  mutable CAROM::Vector zT;  
  const CAROM::Matrix *S;
  
  mutable DenseMatrix J;

  bool hyperreduce, hyperreduce_source;
  bool sourceFOM;

  CAROM::Vector *ufom_librom, *ufom_R_librom, *ufom_W_librom;
  Vector *ufom;
  Vector *ufom_R;
  Vector *ufom_W;
  mutable Vector zfomR;
  mutable Vector zfomW;
  CAROM::Vector *zfomR_librom;
  mutable CAROM::Vector VtzR;

  void PrintFDJacobian(const Vector &u) const;
  
protected:
  CAROM::Matrix* BR;
  CAROM::Matrix* CR;
  const CAROM::Matrix* U_R;
  Vector y0;
  Vector dydt_prev;
  NonlinearDiffusionOperator *fom;
  NonlinearDiffusionOperator *fomSp;
  
public:
  RomOperator(NonlinearDiffusionOperator *fom_, NonlinearDiffusionOperator *fomSp_,
	      const int rrdim_, const int rwdim_, const int nldim_,
	      const CAROM::Matrix* V_R_, const CAROM::Matrix* U_R_, const CAROM::Matrix* V_W_,
	      const CAROM::Matrix *Bsinv, const int N1,
	      const double newton_rel_tol, const double newton_abs_tol, const int newton_iter, const vector<int>& s2sp,
	      const CAROM::Matrix* S_, const vector<int>& s2sp_S_, const CAROM::Matrix *Ssinv_,
	      const vector<int>& st2sp, const vector<int>& sprows, const vector<int>& all_sprows, const int myid,
	      const bool hyperreduce_source);

  virtual void Mult(const Vector &y, Vector &dy_dt) const;
  void Mult_Hyperreduced(const Vector &y, Vector &dy_dt) const;
  void Mult_FullOrder(const Vector &y, Vector &dy_dt) const;

  /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
      This is the only requirement for high-order SDIRK implicit integration.*/
  virtual void ImplicitSolve(const double dt, const Vector &y, Vector &dy_dt);

  virtual Operator &GetGradient(const Vector &y) const;

  CAROM::Matrix V_W, V_R, VTU_R;

  CAROM::Matrix VTCS_W;
  
  virtual ~RomOperator();
};

// TODO: remove this by making online computation serial?
void BroadcastUndistributedRomVector(CAROM::Vector* v)
{
  const int N = v->dim();

  MFEM_VERIFY(N > 0, "");
  
  double *d = new double[N];

  MFEM_VERIFY(d != 0, "");
  
  for (int i=0; i<N; ++i)
    d[i] = (*v)(i);

  MPI_Bcast(d, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for (int i=0; i<N; ++i)
    (*v)(i) = d[i];

  delete d;
}

double InitialTemperature(const Vector &x);
double SourceFunction(const Vector &x, const double t);
double ExactSolution(const Vector &x, const double t);
double NonlinearCoefficient(const double u);
double NonlinearCoefficientDerivative(const double u);

// TODO: move this to the library?
CAROM::Matrix* GetFirstColumns(const int N, const CAROM::Matrix* A)
{
  CAROM::Matrix* S = new CAROM::Matrix(A->numRows(), std::min(N, A->numColumns()), A->distributed());
  for (int i=0; i<S->numRows(); ++i)
    {
      for (int j=0; j<S->numColumns(); ++j)
	(*S)(i,j) = (*A)(i,j);
    }

  // delete A;  // TODO: find a good solution for this.
  return S;
}

void MergeBasis(const int dimFOM, const int nparam, const int max_num_snapshots, std::string name)
{
  MFEM_VERIFY(nparam > 0, "Must specify a positive number of parameter sets");

  bool update_right_SV = false;
  bool isIncremental = false;

  CAROM::Options options(dimFOM, nparam * max_num_snapshots, 1, update_right_SV);
  CAROM::BasisGenerator generator(options, isIncremental, "basis" + name);
      
  for (int paramID=0; paramID<nparam; ++paramID)
    {
      std::string snapshot_filename = "basis" + std::to_string(paramID) + "_" + name + "_snapshot";
      generator.loadSamples(snapshot_filename,"snapshot");
    }

  generator.endSamples(); // save the merged basis file
}

int main(int argc, char *argv[])
{
  // 1. Initialize MPI.
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // 2. Parse command-line options.
  nonlinear_problem = true;
  problem = ANALYTIC;
  diffusion_c = 2.0;
  step_half = 0.25;
  const char *mesh_file = "../data/star.mesh";
  int ser_ref_levels = 2;
  int par_ref_levels = 1;
  int order = 0;
  double t_final = 1.0e-1;
  double maxdt = 1.0e-3;
  double dt = maxdt;
  bool visualization = false;
  bool visit = false;
  int vis_steps = 5;
  double newton_rel_tol = 1e-4;
  double newton_abs_tol = 1e-10;
  int newton_iter = 10;

  // ROM parameters
  bool offline = false;
  bool merge = false;
  bool online = false;

  int nsets = 0;

  int id_param = 0;
  
  int rdim = -1;  // number of basis vectors to use
  int rrdim = -1;  // number of basis vectors to use
  int rwdim = -1;  // number of basis vectors to use

  int precision = 16;
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
  args.AddOption(&id_param, "-id", "--id", "Parametric index");
  args.AddOption(&problem, "-p", "--problem", "Problem setup to use");
  args.AddOption(&step_half, "-sh", "--stephalf", "Initial step function half-width");
  args.AddOption(&diffusion_c, "-dc", "--diffusion-constant", "Diffusion coefficient constant term");
  args.AddOption(&rrdim, "-rrdim", "--rrdim",
		 "Basis dimension for H(div) vector finite element space.");
  args.AddOption(&rwdim, "-rwdim", "--rwdim",
		 "Basis dimension for L2 scalar finite element space.");
  args.AddOption(&t_final, "-tf", "--t-final",
		 "Final time; start time is 0.");
  args.AddOption(&dt, "-dt", "--time-step",
		 "Time step.");
  args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
		 "--no-visualization",
		 "Enable or disable GLVis visualization.");
  args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
		 "--no-visit-datafiles",
		 "Save data files for VisIt (visit.llnl.gov) visualization.");
  args.AddOption(&vis_steps, "-vs", "--visualization-steps",
		 "Visualize every n-th timestep.");
  args.AddOption(&nsets, "-ns", "--nset", "Number of parametric snapshot sets");
  args.AddOption(&offline, "-offline", "--offline", "-no-offline", "--no-offline",
		 "Enable or disable the offline phase.");
  args.AddOption(&online, "-online", "--online", "-no-online", "--no-online",
		 "Enable or disable the online phase.");
  args.AddOption(&merge, "-merge", "--merge", "-no-merge", "--no-merge",
		 "Enable or disable the merge phase.");

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

  const bool check = (offline && !merge && !online) || (!offline && merge && !online) || (!offline && !merge && online);
  MFEM_VERIFY(check, "only one of offline, merge, or online must be true!");

  const bool hyperreduce_source = (problem != INIT_STEP);

  StopWatch solveTimer, totalTimer;
  totalTimer.Start();

  // 3. Read the serial mesh from the given mesh file on all processors. We can
  //    handle triangular, quadrilateral, tetrahedral and hexahedral meshes
  //    with the same code.
  Mesh *mesh = new Mesh(mesh_file, 1, 1);
  const int dim = mesh->Dimension();

  // 4. Define the ODE solver used for time integration. For this example,
  //    only backward Euler is currently supported.
  BackwardEulerSolver ode_solver;

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

  // 7. Define the mixed finite element spaces.

  RT_FECollection hdiv_coll(order, dim);
  L2_FECollection l2_coll(order, dim);

  ParFiniteElementSpace R_space(pmesh, &hdiv_coll);
  ParFiniteElementSpace W_space(pmesh, &l2_coll);
  
  HYPRE_Int dimW = W_space.GlobalTrueVSize();
  HYPRE_Int dimR = R_space.GlobalTrueVSize();

  if (myid == 0)
    {
      cout << "Global number of L2 unknowns: " << dimW << endl;
      cout << "Global number of RT unknowns: " << dimR << endl;
    }

  bool update_right_SV = false;
  bool isIncremental = false;
  const std::string basisFileName = "basis" + std::to_string(id_param);
  int max_num_snapshots = int(t_final/dt) + 2;

  // The merge phase
  if (merge) 
    {
      totalTimer.Clear();
      totalTimer.Start();

      MergeBasis(R_space.GetTrueVSize(), nsets, max_num_snapshots, "R");
      MergeBasis(R_space.GetTrueVSize(), nsets, max_num_snapshots, "FR");
      MergeBasis(W_space.GetTrueVSize(), nsets, max_num_snapshots, "W");

      if (hyperreduce_source)
	MergeBasis(W_space.GetTrueVSize(), nsets, max_num_snapshots, "S");

      totalTimer.Stop();
      if (myid == 0) 
      {
        printf("Elapsed time for merging and building ROM basis: %e second\n", totalTimer.RealTime());
      }
      MPI_Finalize();
      return 0;
    }
  
  ParGridFunction u_gf(&W_space);

  // 8. Set the initial conditions for u.

  const bool SchurComplement = false;

  FunctionCoefficient u_0(InitialTemperature);
  u_gf.ProjectCoefficient(u_0);
  Vector u, uprev, dudt, source;
  Vector *u_W = &u;
  ParGridFunction* sp_u_gf = 0;
  Vector sp_u;
  Vector *sp_u_W = &sp_u;

  Vector *wMFEM = 0;
  
  CAROM::Vector *u_librom = 0;
  CAROM::Vector *u_W_librom = 0;
  CAROM::Vector* w_W = 0;
  CAROM::Vector* w = 0;
  
  const int N1 = R_space.GetTrueVSize();
  const int N2 = W_space.GetTrueVSize();
  const int fdim = N1 + N2;

  cout << myid << ": Local number of L2 unknowns: " << N2 << endl;
  cout << myid << ": Local number of RT unknowns: " << N1 << endl;
  
  if (SchurComplement)
    u_gf.GetTrueDofs(u);
  else
    {
      u_librom = new CAROM::Vector(fdim, true);
      u.SetDataAndSize(&((*u_librom)(0)), fdim);
      u_W_librom = new CAROM::Vector(&((*u_librom)(N1)), N2, true, false);

      u = 0.0;
      u_W = new Vector(u.GetData() + N1, N2);
      u_gf.GetTrueDofs(*u_W);

      source.SetSize(N2);
    }

  // 9. Initialize the diffusion operator and the VisIt visualization.
  NonlinearDiffusionOperator oper(R_space, W_space, newton_rel_tol, newton_abs_tol, newton_iter, u, SchurComplement);  // FOM operator
  NonlinearDiffusionOperator *soper = 0;  // Sample mesh operator
  
  if (SchurComplement)
    u_gf.SetFromTrueDofs(*u_W);

  if (offline)
  {
    ostringstream mesh_name, sol_name;
    mesh_name << "nldiff-mesh." << setfill('0') << setw(6) << myid;
    sol_name << "nldiff-init" << id_param << "." << setfill('0') << setw(6) << myid;
    ofstream omesh(mesh_name.str().c_str());
    omesh.precision(precision);
    pmesh->Print(omesh);
    ofstream osol(sol_name.str().c_str());
    osol.precision(precision);
    u_gf.Save(osol);
  }

  VisItDataCollection * visit_dc = NULL;
  if (visit)
    {
      if (offline)
	visit_dc = new VisItDataCollection("nldiff-fom", pmesh);
      else
	visit_dc = new VisItDataCollection("nldiff-rom", pmesh);

      visit_dc->RegisterField("temperature", &u_gf);
      visit_dc->SetCycle(0);
      visit_dc->SetTime(0.0);
      visit_dc->Save();
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
	  sout.precision(precision);
	  sout << "solution\n" << *pmesh << u_gf;
	  sout << "pause\n";
	  sout << flush;
	  if (myid == 0)
	    {
	      cout << "GLVis visualization paused."
		   << " Press space (in the GLVis window) to resume it.\n";
	    }
	}
    }

  CAROM::BasisGenerator *basis_generator_R = 0;  // For the solution component in vector H(div)
  CAROM::BasisGenerator *basis_generator_W = 0;  // For the solution component in scalar L2
  CAROM::BasisGenerator *basis_generator_FR = 0; // For the nonlinear term M(u)v with u in L2, v in H(div)
  CAROM::BasisGenerator *basis_generator_S = 0;  // For the source in scalar L2

  if (offline) {
    CAROM::Options options_R(R_space.GetTrueVSize(), max_num_snapshots, 1, update_right_SV);
    CAROM::Options options_W(W_space.GetTrueVSize(), max_num_snapshots, 1, update_right_SV);

    if (hyperreduce_source)
      basis_generator_S = new CAROM::BasisGenerator(options_W, isIncremental, basisFileName + "_S");

    basis_generator_R = new CAROM::BasisGenerator(options_R, isIncremental, basisFileName + "_R");
    basis_generator_W = new CAROM::BasisGenerator(options_W, isIncremental, basisFileName + "_W");

    basis_generator_FR = new CAROM::BasisGenerator(options_R, isIncremental, basisFileName + "_FR");
  }
  
  vector<int> s2sp;   // mapping from sample dofs in original mesh (s) to stencil dofs in sample mesh (s+)
  vector<int> st2sp;  // mapping from stencil dofs in original mesh (st) to stencil dofs in sample mesh (s+)
  ParMesh* sample_pmesh = 0;
  RomOperator *romop = 0;

  const CAROM::Matrix* B_librom = 0;
  const CAROM::Matrix* BR_librom = 0;
  const CAROM::Matrix* FR_librom = 0;
  const CAROM::Matrix* BW_librom = 0;
  const CAROM::Matrix* S_librom = 0;

  int nsamp_R = -1;
  int nsamp_S = -1;

  if (online)
    {
      CAROM::BasisReader readerR("basisR");
      BR_librom = readerR.getSpatialBasis(0.0);
      if (rrdim == -1)
	rrdim = BR_librom->numColumns();
      else
	BR_librom = GetFirstColumns(rrdim, BR_librom);  // TODO: reduce rrdim if too large

      MFEM_VERIFY(BR_librom->numRows() == N1, "");
      
      if (myid == 0)
	printf("reduced R dim = %d\n",rrdim);

      CAROM::BasisReader readerW("basisW");
      BW_librom = readerW.getSpatialBasis(0.0);
      if (rwdim == -1)
	rwdim = BW_librom->numColumns();
      else
	BW_librom = GetFirstColumns(rwdim, BW_librom);

      MFEM_VERIFY(BW_librom->numRows() == N2, "");

      if (myid == 0)
	printf("reduced W dim = %d\n",rwdim);

      // TODO: To get the basis U_R, considering the relation M(u) v + B^T u = 0 in the FOM, we can just use
      // U_R = B^T V_W. Note that V_W and V_R may have different numbers of columns, which is fine.
      // TODO: maybe we need POD applied to B^T multiplied by W-snapshots, or just a basis generator for
      // snapshots of M(u) v. This could be different from B^T multiplied by POD results for the W solutions.

      /*
      FR_librom = new CAROM::Matrix(N1, rwdim, true);
      oper.SetBTV(BW_librom, FR_librom);
      */
      
      CAROM::BasisReader readerFR("basisFR");
      FR_librom = readerFR.getSpatialBasis(0.0);

      // Compute sample points using DEIM, for hyperreduction

      // TODO: reduce this?
      const int nldim = FR_librom->numColumns(); // rwdim;

      cout << "FR dim " << FR_librom->numColumns() << endl;
      
      MFEM_VERIFY(FR_librom->numRows() == N1 && FR_librom->numColumns() >= nldim, "");

      if (FR_librom->numColumns() > nldim)
	FR_librom = GetFirstColumns(nldim, FR_librom);
      
      vector<int> num_sample_dofs_per_proc(num_procs);

      nsamp_R = nldim;

#ifdef USE_GNAT
      vector<int> sample_dofs(nsamp_R);  // Indices of the sampled rows
      CAROM::Matrix *Bsinv = new CAROM::Matrix(nsamp_R, nldim, false);
      CAROM::GNAT(FR_librom,
                  nldim,
                  &sample_dofs[0],
                  &num_sample_dofs_per_proc[0],
                  *Bsinv,
                  myid,
                  num_procs,
		  nsamp_R);
#else
      // Now execute the DEIM algorithm to get the sampling information.
      CAROM::Matrix *Bsinv = new CAROM::Matrix(nldim, nldim, false);
      vector<int> sample_dofs(nldim);  // Indices of the sampled rows
      CAROM::DEIM(FR_librom,
                  nldim,
                  &sample_dofs[0],
                  &num_sample_dofs_per_proc[0],
                  *Bsinv,
                  myid,
                  num_procs);
#endif

      vector<int> sample_dofs_withS;  // Indices of the sampled rows
      int nsdim = 0;
      int *allNR = 0;
      CAROM::Matrix *Ssinv = 0;
      vector<int> num_sample_dofs_per_proc_withS;
      CAROM::BasisReader *readerS = 0;
      if (hyperreduce_source)
	{
	  readerS = new CAROM::BasisReader("basisS");
	  S_librom = readerS->getSpatialBasis(0.0);

	  // Compute sample points using DEIM

	  nsdim = S_librom->numColumns();

	  cout << "S dim " << nsdim << endl;
      
	  vector<int> num_sample_dofs_per_proc_S(num_procs);

	  // Now execute the DEIM algorithm to get the sampling information.
	  nsamp_S = nsdim;

#ifdef USE_GNAT
	  Ssinv = new CAROM::Matrix(nsamp_S, nsdim, false);
	  vector<int> sample_dofs_S(nsamp_S);  // Indices of the sampled rows

	  CAROM::GNAT(S_librom,
		      nsdim,
		      &sample_dofs_S[0],
		      &num_sample_dofs_per_proc_S[0],
		      *Ssinv,
		      myid,
		      num_procs,
		      nsamp_S);
#else
	  Ssinv = new CAROM::Matrix(nsdim, nsdim, false);
	  vector<int> sample_dofs_S(nsdim);  // Indices of the sampled rows
	  CAROM::DEIM(S_librom,
		      nsdim,
		      &sample_dofs_S[0],
		      &num_sample_dofs_per_proc_S[0],
		      *Ssinv,
		      myid,
		      num_procs);
#endif

	  // Merge sample_dofs and sample_dofs_S
	  allNR = new int [num_procs];
	  MPI_Allgather(&N1, 1, MPI_INT, allNR, 1, MPI_INT, MPI_COMM_WORLD);

	  sample_dofs_withS.resize(sample_dofs.size() + sample_dofs_S.size());  // Indices of the sampled rows
	  num_sample_dofs_per_proc_withS.resize(num_procs);
	  int offset = 0;
	  int offset_S = 0;
	  for (int p=0; p<num_procs; ++p)
	    {
	      for (int i=0; i<num_sample_dofs_per_proc[p]; ++i)
		sample_dofs_withS[offset + offset_S + i] = sample_dofs[offset + i];

	      offset += num_sample_dofs_per_proc[p];
	  
	      for (int i=0; i<num_sample_dofs_per_proc_S[p]; ++i)
		sample_dofs_withS[offset + offset_S + i] = allNR[p] + sample_dofs_S[offset_S + i];

	      offset_S += num_sample_dofs_per_proc_S[p];
	      num_sample_dofs_per_proc_withS[p] = num_sample_dofs_per_proc[p] + num_sample_dofs_per_proc_S[p];
	    }
	}

      // Construct sample mesh

      // TODO: put this in CreateSampleMesh!
      // Define a superfluous finite element space, merely to get global vertex indices for the sample mesh construction.
      H1_FECollection h1_coll(1, dim);  // Must be first order, to get a bijection between vertices and DOF's.
      ParFiniteElementSpace H1_space(pmesh, &h1_coll);  // This constructor effectively sets vertex (DOF) global indices.

      ParFiniteElementSpace *sp_R_space, *sp_W_space;

      MPI_Comm rom_com;
      int color = myid != 0;
      const int status = MPI_Comm_split(MPI_COMM_WORLD, color, myid, &rom_com);
      MFEM_VERIFY(status == MPI_SUCCESS,
                  "Construction of hyperreduction comm failed");

      vector<int> sprows;
      vector<int> all_sprows;
      vector<int> s2sp_S;

      if (hyperreduce_source)
	{
	  vector<int> s2sp_withS;
	  CAROM::CreateSampleMesh(*pmesh, H1_space, R_space, W_space, hdiv_coll, l2_coll, rom_com, sample_dofs_withS, num_sample_dofs_per_proc_withS,
				  sample_pmesh, sprows, all_sprows, s2sp_withS, st2sp, sp_R_space, sp_W_space);

	  if (myid == 0)
	    {
	      s2sp_S.resize(nsamp_S);

#ifndef USE_GNAT
	      MFEM_VERIFY(s2sp_withS.size() == nldim + nsdim, "");
#endif
	      MFEM_VERIFY(s2sp_withS.size() == nsamp_R + nsamp_S, "");

	      ofstream sfile("smesh");
	      sample_pmesh->Print(sfile);
	      sfile.close();

	      const int NRsp = sp_R_space->GetTrueVSize();

	      s2sp.resize(nsamp_R);

	      int count = 0;
	      int count_S = 0;

	      int offset = 0;
	      for (int p=0; p<num_procs; ++p)
		{
		  for (int i=0; i<num_sample_dofs_per_proc_withS[p]; ++i)
		    {
		      if (sample_dofs_withS[offset + i] >= allNR[p])
			{
			  s2sp_S[count_S] = s2sp_withS[offset + i] - NRsp;
			  count_S++;
			}
		      else
			{
			  s2sp[count] = s2sp_withS[offset + i];
			  count++;
			}
		    }

		  offset += num_sample_dofs_per_proc_withS[p];
		}
	  
	      //MFEM_VERIFY(count == nldim && count_S == nsdim, "");
	      MFEM_VERIFY(count == nsamp_R && count_S == nsamp_S, "");
	    }
	}
      else
	{
	  CAROM::CreateSampleMesh(*pmesh, H1_space, R_space, W_space, hdiv_coll, l2_coll, rom_com, sample_dofs, num_sample_dofs_per_proc,
				  sample_pmesh, sprows, all_sprows, s2sp, st2sp, sp_R_space, sp_W_space);
	}

      w = new CAROM::Vector(rrdim + rwdim, false);
      w_W = new CAROM::Vector(rwdim, false);
      
      // Initialize w = B_W^T u.
      BW_librom->transposeMult(*u_W_librom, *w_W);

      for (int i=0; i<rrdim; ++i)
	(*w)(i) = 0.0;

      for (int i=0; i<rwdim; ++i)
	(*w)(rrdim + i) = (*w_W)(i);

      // Note that some of this could be done only on the ROM solver process, but it is tricky, since RomOperator assembles Bsp in parallel.
      wMFEM = new Vector(&((*w)(0)), rrdim + rwdim);
	  
      if (myid == 0)
	{
	  // Initialize sp_u with initial conditions.
	  if (SchurComplement)
	    {
	      MFEM_VERIFY(false, "Schur complement formulation cannot be used with ROM.");
	    }
	  else
	    {
	      sp_u_gf = new ParGridFunction(sp_W_space);
	      sp_u_gf->ProjectCoefficient(u_0);

	      sp_u.SetSize(sp_R_space->GetTrueVSize() + sp_W_space->GetTrueVSize());
	      sp_u = 0.0;
	      sp_u_W = new Vector(sp_u.GetData() + sp_R_space->GetTrueVSize(), sp_W_space->GetTrueVSize());
	      sp_u_gf->GetTrueDofs(*sp_u_W);
	    }
 
	  soper = new NonlinearDiffusionOperator(*sp_R_space, *sp_W_space, newton_rel_tol, newton_abs_tol, newton_iter, sp_u, SchurComplement);
	}

      romop = new RomOperator(&oper, soper, rrdim, rwdim, nldim,
			      BR_librom, FR_librom, BW_librom,
			      Bsinv, N1, newton_rel_tol, newton_abs_tol, newton_iter, s2sp,
			      S_librom, s2sp_S, Ssinv,
			      st2sp, sprows, all_sprows, myid, hyperreduce_source);
      
      ode_solver.Init(*romop);

      delete readerS;
    }
  else  // fom
    ode_solver.Init(oper);

  // 10. Perform time-integration (looping over the time iterations, ti, with a time-step dt).
  double t = 0.0;
  int lastStepChange = 1;

  ConstantCoefficient coeff0(0.0);

  oper.newtonFailure = false;

  solveTimer.Start();

  bool last_step = false;
  for (int ti = 1; !last_step; ti++)
    {
      if (offline && !oper.newtonFailure)
	{
	  const bool sampleW = basis_generator_W->isNextSample(t);
	  
	  if (sampleW && hyperreduce_source) // TODO: Instead, basis_generator_S->isNextSample(t) could be used if dS/dt were computed.
	    {
	      oper.GetSource(source);
	      basis_generator_S->takeSample(source.GetData(), t, dt);
	      // TODO: dfdt? In this example, one can implement the exact formula.
	      //   In general, one can use finite differences in time (dudt is computed that way). 
	      //basis_generator_S->computeNextSampleTime(u.GetData(), dfdt.GetData(), t);
	    }

	  if (basis_generator_R->isNextSample(t))
	    {
	      oper.CopyDuDt(dudt);

	      basis_generator_R->takeSample(u.GetData(), t, dt);
	      basis_generator_R->computeNextSampleTime(u.GetData(), dudt.GetData(), t);

	      Vector u_R(u.GetData(), N1);
	      Vector Mu(N1);
	      oper.SetParameters(u);
	      oper.Mult_Mmat(u_R, Mu);
	      basis_generator_FR->takeSample(Mu.GetData(), t, dt);
	    }

	  if (sampleW)
	    {
	      oper.CopyDuDt_W(dudt);

	      basis_generator_W->takeSample(u_W->GetData(), t, dt);
	      basis_generator_W->computeNextSampleTime(u_W->GetData(), dudt.GetData(), t);
	    }
	}

      if (online)
	{
	  if (myid == 0)
	    {
	      ode_solver.Step(*wMFEM, t, dt);
	    }

	  MPI_Bcast(&t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
      else  // fom
	{
	  oper.newtonFailure = false;
	  uprev = u;  // Save solution, to reset in case of a Newton failure.
	  const double tprev = t;
	  ode_solver.Step(u, t, dt);

	  if (oper.newtonFailure)
	    {
	      // Reset and retry.
	      u = uprev;
	      t = tprev;
	      dt *= 0.5;
	      cout << "step " << ti << ", t = " << t << " had a Newton failure, cutting dt to " << dt << endl;
	      lastStepChange = ti;
	      ti--;
	      continue;
	    }
	  else if (ti - lastStepChange > 10)
	    {
	      dt *= 2.0;
	      dt = min(dt, maxdt);
	      lastStepChange = ti;
	    }
	}

      if (myid == 0)
	{
	  cout << "step " << ti << ", t = " << t << " == " << oper.GetTime() << ", dt " << dt << endl;

	  if (online)
	    cout << "rom time " << romop->GetTime() << endl;
	}
      
      if (t >= t_final - dt/2)
	last_step = true;

      if (last_step || (ti % vis_steps) == 0)
	{
	  FunctionCoefficient exsol(ExactSolution);
	  exsol.SetTime(oper.GetTime());

	  if (online)  // Lift ROM solution to FOM.
	    {
	      exsol.SetTime(t);
	      
	      BroadcastUndistributedRomVector(w);
	      
	      for (int i=0; i<rwdim; ++i)
		(*w_W)(i) = (*w)(rrdim + i);
	      
	      romop->V_W.mult(*w_W, *u_W_librom);

	      u_gf.SetFromTrueDofs(*u_W);

	      if (last_step)
		{
		  // Calculate the relative l2 error between the final ROM solution and FOM solution, using id_param for FOM solution.
		  Vector fom_solution(N2);
		  ifstream solution_file;
		  ostringstream solution_filename, rom_filename;
		  solution_filename << "nldiff-fom-values-final" << id_param << "." << setfill('0') << setw(6) << myid;
		  rom_filename << "nldiff-rom-final" << id_param << "." << setfill('0') << setw(6) << myid;

		  if (myid == 0) std::cout << "Comparing current run to solution at: " << solution_filename.str() << " with offline parameter index " << id_param << std::endl;
		  solution_file.open(solution_filename.str());
		  fom_solution.Load(solution_file, N2);
		  solution_file.close();
		  const double fomNorm = sqrt(InnerProduct(MPI_COMM_WORLD, fom_solution, fom_solution));
		  //const double romNorm = sqrt(InnerProduct(MPI_COMM_WORLD, *u_W, *u_W));
		  fom_solution -= *u_W;
		  const double diffNorm = sqrt(InnerProduct(MPI_COMM_WORLD, fom_solution, fom_solution));
		  if (myid == 0) std::cout << "Relative l2 error of ROM solution " << diffNorm / fomNorm << std::endl;

		  ofstream osol(rom_filename.str().c_str());
		  osol.precision(precision);
		  u_gf.Save(osol);
		}
	    }
	  else
	    u_gf.SetFromTrueDofs(*u_W);

	  if (problem == ANALYTIC)
	    {
	      const double l2err = u_gf.ComputeL2Error(exsol);
	      const double l2nrm = u_gf.ComputeL2Error(coeff0);

	      if (myid == 0)
		cout << "L2 norm of exact error: " << l2err << ", FEM solution norm " << l2nrm << ", relative norm " << l2err / l2nrm << endl;
	    }

	  if (visualization)
	    {
	      sout << "parallel " << num_procs << " " << myid << "\n";
	      sout << "solution\n" << *pmesh << u_gf << flush;
	    }

	  if (visit)
	    {
	      visit_dc->SetCycle(ti);
	      visit_dc->SetTime(t);
	      visit_dc->Save();
	    }
	}
    }  // timestep loop

  solveTimer.Stop();
  if (myid == 0) cout << "Elapsed time for time integration loop " << solveTimer.RealTime() << endl;

  if (visit)
    delete visit_dc;

  if (offline)
    {
      // Sample final solution, to prevent extrapolation in ROM between the last sample and the end of the simulation.

      oper.CopyDuDt(dudt);

      // R space
      basis_generator_R->takeSample(u.GetData(), t, dt);

      Vector u_R(u.GetData(), N1);
      Vector Mu(N1);
      oper.SetParameters(u);
      oper.Mult_Mmat(u_R, Mu);
      basis_generator_FR->takeSample(Mu.GetData(), t, dt);
      
      // Terminate the sampling and write out information.
      basis_generator_R->writeSnapshot();
      basis_generator_FR->writeSnapshot();

      int rom_dim = basis_generator_R->getSpatialBasis()->numColumns();
      cout << "R-space ROM dimension = " << rom_dim << endl;

      const CAROM::Vector* sing_vals_R = basis_generator_R->getSingularValues();

      if (myid == 0)
	{
            
	  cout << "Singular Values:" << endl;
	  for (int sv = 0; sv < sing_vals_R->dim(); ++sv) {
            double this_sv = (*sing_vals_R)(sv);
            cout << this_sv << endl;
	  }
	}

      // W space

      // TODO: why call computeNextSampleTime if you just do takeSample on every step anyway?
      basis_generator_W->takeSample(u_W->GetData(), t, dt);
      basis_generator_W->writeSnapshot();

      oper.GetSource(source);

      if (hyperreduce_source)
	{
	  basis_generator_S->takeSample(source.GetData(), t, dt);
	  basis_generator_S->writeSnapshot();
	}
      
      rom_dim = basis_generator_W->getSpatialBasis()->numColumns();
      cout << "W-space ROM dimension = " << rom_dim << endl;

      const CAROM::Vector* sing_vals_W = basis_generator_W->getSingularValues();

      if (myid == 0)
	{
	  cout << "Singular Values:" << endl;
	  for (int sv = 0; sv < sing_vals_W->dim(); ++sv) {
            double this_sv = (*sing_vals_W)(sv);
            cout << this_sv << endl;
	  }
	}

      if (hyperreduce_source)
	{
	  rom_dim = basis_generator_S->getSpatialBasis()->numColumns();
	  const CAROM::Vector* sing_vals_S = basis_generator_S->getSingularValues();

	  if (myid == 0)
	    {
	      cout << "S-space ROM dimension = " << rom_dim << endl;

	      cout << "Singular Values:" << endl;
	      for (int sv = 0; sv < sing_vals_S->dim(); ++sv) {
		double this_sv = (*sing_vals_S)(sv);
		cout << this_sv << endl;
	      }
	    }
	}

      delete basis_generator_R;
      delete basis_generator_FR;
      delete basis_generator_W;
      delete basis_generator_S;      
    }

  // 11. Save the final solution in parallel. This output can be viewed later
  //     using GLVis: "glvis -np <np> -m nldiff-mesh -g nldiff-final".
  if (offline)
  {
    ostringstream sol_name, fomsol_name;
    sol_name << "nldiff-final" << id_param << "." << setfill('0') << setw(6) << myid;
    ofstream osol(sol_name.str().c_str());
    osol.precision(precision);
    u_gf.Save(osol);

    fomsol_name << "nldiff-fom-values-final" << id_param << "." << setfill('0') << setw(6) << myid;
    ofstream fomsol(fomsol_name.str().c_str());
    fomsol.precision(precision);
    for (int i = 0; i < N2; ++i)
      {
	fomsol << (*u_W)[i] << std::endl;
      }
  }

  // 12. Free the used memory.
  delete pmesh;
  delete romop;
  
  if (!SchurComplement)
    delete u_W;

  totalTimer.Stop();
  if (myid == 0) cout << "Elapsed time for entire simulation " << totalTimer.RealTime() << endl;

  MPI_Finalize();
  return 0;
}

NonlinearDiffusionOperator::NonlinearDiffusionOperator(ParFiniteElementSpace &fR, ParFiniteElementSpace &fW,
						       const double rel_tol, const double abs_tol, 
						       const int iter, const Vector &u, const bool SchurComplement_)
  : TimeDependentOperator(SchurComplement_ ? fW.GetTrueVSize() : fR.GetTrueVSize() + fW.GetTrueVSize(), 0.0), 
    fespace_R(fR), fespace_W(fW), M(NULL), C(NULL), Bmat(NULL), BTmat(NULL), Mprime(NULL), current_dt(0.0), 
    newton_solver(fW.GetComm()), M_solver(NULL), C_solver(fW.GetComm()), zW(fW.GetTrueVSize()), yR(fR.GetTrueVSize()),
    zR(fR.GetTrueVSize()), u0(height), dudt_prev(height),
    SchurComplement(SchurComplement_), fullOp(NULL), fullGradient(NULL), fullPrec(NULL)
{
  gradient = new NonlinearDiffusionGradientOperator(fR.GetTrueVSize(), fW.GetTrueVSize());

  linear_solver_rel_tol = 1.0e-14;

  C = new ParBilinearForm(&fespace_W);
  C->AddDomainIntegrator(new MassIntegrator());
  C->Assemble();
  C->Finalize();
  Cmat = C->ParallelAssemble();

  C_solver.iterative_mode = false;
  C_solver.SetRelTol(linear_solver_rel_tol);
  C_solver.SetAbsTol(0.0);
  C_solver.SetMaxIter(1000);
  C_solver.SetPrintLevel(0);
  C_prec.SetType(HypreSmoother::Jacobi);
  C_solver.SetPreconditioner(C_prec);
  C_solver.SetOperator(*Cmat);

  ParMixedBilinearForm *bVarf(new ParMixedBilinearForm(&fespace_R, &fespace_W));
  bVarf->AddDomainIntegrator(new VectorFEDivergenceIntegrator);
  bVarf->Assemble();
  bVarf->Finalize();
  Bmat = bVarf->ParallelAssemble();
  (*Bmat) *= -1.0;

  BTmat = Bmat->Transpose();

  delete bVarf;

  // See ex5p.cpp for an explanation of block_trueOffsets
  block_trueOffsets.SetSize(3); // number of variables + 1
  block_trueOffsets[0] = 0;
  block_trueOffsets[1] = fespace_R.GetTrueVSize();
  block_trueOffsets[2] = fespace_W.GetTrueVSize();
  block_trueOffsets.PartialSum();

  // SetParameters(u);

  M_prec.SetType(HypreSmoother::Jacobi);

  // Set the newton solve parameters

  //Solver *J_prec = new DSmoother(1);
  J_gmres = new GMRESSolver(MPI_COMM_WORLD);
  J_gmres->SetRelTol(linear_solver_rel_tol);
  J_gmres->SetAbsTol(0.0);
  J_gmres->SetMaxIter(1000);
  J_gmres->SetPrintLevel(2);
  // TODO: precondition, for an efficient FOM solver.
  // J_gmres->SetPreconditioner(*J_prec);

  newton_solver.iterative_mode = true;
  newton_solver.SetSolver(*J_gmres);
  newton_solver.SetOperator(*this);
  newton_solver.SetPrintLevel(1);
  newton_solver.SetRelTol(rel_tol);
  newton_solver.SetAbsTol(abs_tol);
  newton_solver.SetMaxIter(iter);

  dudt_prev = 0.0;
}

void NonlinearDiffusionOperator::SetBTV(const CAROM::Matrix *V, CAROM::Matrix *BTV) const 
{
  const int ncol = BTV->numColumns();
  const int nw = zW.Size();
  const int nr = zR.Size();
  
  MFEM_VERIFY(V->numRows() == nw && BTV->numRows() == nr && V->numColumns() >= ncol, "");
  
  for (int k=0; k<ncol; ++k)
    {
      for (int i=0; i<nw; ++i)
	zW[i] = (*V)(i,k);

      BTmat->Mult(zW, zR);

      for (int i=0; i<nr; ++i)
	(*BTV)(i,k) = zR[i];
    }
}
 
void NonlinearDiffusionOperator::Mult(const Vector &du_dt, Vector &res) const
{
  if (SchurComplement)
    Mult_SchurComplement(du_dt, res);
  else 
    Mult_FullSystem(du_dt, res);
}

void NonlinearDiffusionOperator::Mult_SchurComplement(const Vector &du_dt, Vector &res) const
{
  // Compute:
  //    du_dt - C^{-1} (f - B M(u)^{-1} B^T u), with u = u0 + dt*du_dt

  // Set grid function for f
  ParGridFunction f_gf(&fespace_W);

  FunctionCoefficient f(SourceFunction);
  f.SetTime(GetTime());
  f_gf.ProjectCoefficient(f);
  Vector fproj;
  f_gf.GetTrueDofs(fproj);

  Vector u(u0);
  u.Add(current_dt, du_dt);

  SetParameters(u);  // Create M(a(u)), M(aprime(u))

  // Compute C^{-1} (f - B M(u)^{-1} B^T u)
  Bmat->MultTranspose(u, zR);
  M_solver->Mult(zR, yR);
  Bmat->Mult(yR, zW);

  res = du_dt;
  res.Add(-1.0, fproj);

  C_solver.Mult(zW, fproj);

  res.Add(1.0, fproj);
}

void NonlinearDiffusionOperator::GetSource(Vector& s) const
{
  // Set grid function for f
  ParGridFunction f_gf(&fespace_W);

  FunctionCoefficient f(SourceFunction);
  f.SetTime(GetTime());
  f_gf.ProjectCoefficient(f);
  f_gf.GetTrueDofs(s);
}

void NonlinearDiffusionOperator::Mult_FullSystem(const Vector &du_dt, Vector &res) const
{
  // Compute:
  //    [   Mv + B^T u   ], with u = u0 + dt*du_dt
  //    [C du_dt - Bv - f]       v = v0 + dt*dv_dt

  GetSource(zW);
  
  Vector u(u0);
  u.Add(current_dt, du_dt);

  SetParameters(u);  // Create fullOp

  fullOp->Mult(u, res);  // Sets the first block row of res. The second block row is computed below.

  Vector u_R(u.GetData() + block_trueOffsets[0], block_trueOffsets[1]-block_trueOffsets[0]);
  Vector ut_W(du_dt.GetData() + block_trueOffsets[1], block_trueOffsets[2]-block_trueOffsets[1]);
  Vector res_W(res.GetData() + block_trueOffsets[1], block_trueOffsets[2]-block_trueOffsets[1]);

  res_W = ut_W;
  res_W.Add(-1.0, zW);  // -= f
  Cmat->Mult(res_W, zW);  // = C du_dt - Cf

  res_W = zW;

  Bmat->Mult(u_R, zW);  // Bv
  res_W.Add(-1.0, zW);  // -= Bv
}

void NonlinearDiffusionOperator::ImplicitSolve(const double dt,
					       const Vector &u, Vector &du_dt)
{
  // Solve the equation:
  //    du_dt = C^{-1} (f - B M(u + dt du_dt)^{-1} B^T (u + dt du_dt)), in the Schur complement case
  // for du_dt

  current_dt = dt;
  // MFEM_VERIFY(dt == current_dt, ""); // SDIRK methods use the same dt

  u0 = u;
  
  // Set the initial guess for du_dt, to be used by newton_solver.
  //du_dt = 0.0;
  du_dt = dudt_prev;

  Vector zero; // empty vector is interpreted as zero r.h.s. by NewtonSolver
  newton_solver.Mult(zero, du_dt);

  // MFEM_VERIFY(newton_solver.GetConverged(), "Newton solver did not converge.");
  if (newton_solver.GetConverged())
    dudt_prev = du_dt;
  else
    {
      du_dt = 0.0;  // Zero update in SDIRK Step() function.
      newtonFailure = true;
    }
}

Operator &NonlinearDiffusionOperator::GetGradient(const Vector &u) const
{
  // Note that if a matrix A depends on a parameter t, then dA^{-1}/dt = -A^{-1} dA/dt A^{-1}.
  // (d/du) M(u)^{-1} = -M(u)^{-1} M(a'(u)) M(u)^{-1}

  // Gradient is C^{-1} B M(a(u))^{-1} B^T - C^{-1} B M(u)^{-1} M(a'(u)) M(u)^{-1} B^T, Schur complement case.

  if (SchurComplement)
    return *gradient;
  else
    return *fullGradient;
}

void NonlinearDiffusionOperator::SetParameters(const Vector &u) const 
{
  // Set grid function for a(u)
  ParGridFunction u_gf(&fespace_W);
  ParGridFunction a_gf(&fespace_W);
  ParGridFunction aprime_gf(&fespace_W);
  ParGridFunction a_plus_aprime_gf(&fespace_W);

  if (SchurComplement)
    u_gf.SetFromTrueDofs(u);
  else
    {
      Vector u_W(u.GetData() + block_trueOffsets[1], block_trueOffsets[2]-block_trueOffsets[1]);
      u_gf.SetFromTrueDofs(u_W);
    }

  for (int i = 0; i < a_gf.Size(); i++)
    {
      a_gf(i) = NonlinearCoefficient(u_gf(i));
      aprime_gf(i) = NonlinearCoefficientDerivative(u_gf(i));
      a_plus_aprime_gf(i) = a_gf(i) + aprime_gf(i);
    }

  GridFunctionCoefficient a_coeff(&a_gf);
  GridFunctionCoefficient aprime_coeff(&aprime_gf);
  GridFunctionCoefficient a_plus_aprime_coeff(&a_plus_aprime_gf);

  delete M;
  M = new ParBilinearForm(&fespace_R);

  M->AddDomainIntegrator(new VectorFEMassIntegrator(a_coeff));
  M->Assemble();
  M->Finalize();
  Mmat = M->ParallelAssemble();

  delete Mprime;
  Mprime = new ParBilinearForm(&fespace_R);

  if (SchurComplement)
    Mprime->AddDomainIntegrator(new VectorFEMassIntegrator(aprime_coeff));
  else
    Mprime->AddDomainIntegrator(new VectorFEMassIntegrator(a_plus_aprime_coeff));

  Mprime->Assemble();
  Mprime->FormSystemMatrix(ess_Rdof_list, Mprimemat);

  Mprimemat *= current_dt;

  delete M_solver;
  M_solver = new CGSolver(fespace_R.GetComm());

  M_solver->iterative_mode = false;
  M_solver->SetRelTol(linear_solver_rel_tol);
  M_solver->SetAbsTol(0.0);
  M_solver->SetMaxIter(1000);
  M_solver->SetPrintLevel(-1);
  M_solver->SetPreconditioner(M_prec);

  if (SchurComplement)
    {
      M_solver->SetOperator(*Mmat);
      gradient->SetParameters(&C_solver, M_solver, &Mprimemat, Bmat, current_dt);
    }
  else
    {
      M_solver->SetOperator(Mprimemat);

      delete fullOp;
      fullOp = new BlockOperator(block_trueOffsets);
      fullOp->SetBlock(0, 0, Mmat);
      fullOp->SetBlock(0, 1, BTmat);
      fullOp->SetBlock(1, 0, Bmat, -1.0);
      fullOp->SetBlock(1, 1, Cmat);

      delete fullGradient;
      fullGradient = new BlockOperator(block_trueOffsets);
      fullGradient->SetBlock(0, 0, &Mprimemat);
      fullGradient->SetBlock(0, 1, BTmat, current_dt);
      fullGradient->SetBlock(1, 0, Bmat, -current_dt);
      fullGradient->SetBlock(1, 1, Cmat);

      delete fullPrec;
      fullPrec = new BlockDiagonalPreconditioner(block_trueOffsets);
      fullPrec->SetDiagonalBlock(0, M_solver);
      fullPrec->SetDiagonalBlock(1, &C_solver);
      J_gmres->SetPreconditioner(*fullPrec);
    }
}

NonlinearDiffusionOperator::~NonlinearDiffusionOperator()
{
  delete M;
  delete Mprime;
  delete C;
  delete Bmat;
  delete J_gmres;
}

void Compute_CtAB(const HypreParMatrix* A,
		  const CAROM::Matrix& B,  // Distributed matrix.
		  const CAROM::Matrix& C,  // Distributed matrix.
		  CAROM::Matrix* CtAB)     // Non-distributed (local) matrix, computed identically and redundantly on every process.
{
  MFEM_VERIFY(B.distributed() && C.distributed() && !CtAB->distributed(), "");

  int num_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  const int num_rows = B.numRows();
  const int num_cols = B.numColumns();
  const int num_rows_A = A->GetNumRows();

  MFEM_VERIFY(C.numRows() == num_rows_A, "");

  Vector Bvec(num_rows);
  Vector ABvec(num_rows_A);
  
  CAROM::Matrix AB(num_rows_A, num_cols, true);
  
  for (int i = 0; i < num_cols; ++i) {
    for (int j = 0; j < num_rows; ++j) {
      Bvec[j] = B(j, i);
    }
    A->Mult(Bvec, ABvec);
    for (int j = 0; j < num_rows_A; ++j) {
      AB(j, i) = ABvec[j];
    }
  }
  
  C.transposeMult(AB, CtAB);
}

RomOperator::RomOperator(NonlinearDiffusionOperator *fom_, NonlinearDiffusionOperator *fomSp_, const int rrdim_, const int rwdim_,
			 const int nldim_,
			 const CAROM::Matrix* V_R_, const CAROM::Matrix* U_R_, const CAROM::Matrix* V_W_,
			 const CAROM::Matrix *Bsinv, const int N1,
			 const double newton_rel_tol, const double newton_abs_tol, const int newton_iter, const vector<int>& s2sp_,
			 const CAROM::Matrix* S_, const vector<int>& s2sp_S_, const CAROM::Matrix *Ssinv_,
			 const vector<int>& st2sp, const vector<int>& sprows, const vector<int>& all_sprows, const int myid,
			 const bool hyperreduce_source_)
  : TimeDependentOperator(rrdim_ + rwdim_, 0.0),
    newton_solver(),
    fom(fom_), fomSp(fomSp_), BR(NULL), rrdim(rrdim_), rwdim(rwdim_), nldim(nldim_),
    nsamp_R(s2sp_.size()), nsamp_S(s2sp_S_.size()),
    V_R(*V_R_), U_R(U_R_), V_W(*V_W_), VTU_R(rrdim_, nldim_, false),
    // TODO: get rid of the zero dimension hacks here
    y0(height), dydt_prev(height), zY(nldim, false), zN(std::max(nsamp_R, 1), false), s2sp(s2sp_), Vsinv(Bsinv), J(height),
    s2sp_S(s2sp_S_), zS(std::max(s2sp_S_.size(), std::size_t(1)), false), zT(std::max(s2sp_S_.size(), std::size_t(1)), false), Ssinv(Ssinv_),
    VTCS_W(rwdim, std::max(s2sp_S_.size(), std::size_t(1)), false), S(S_),
    VtzR(rrdim_, false), hyperreduce_source(hyperreduce_source_)
{
  dydt_prev = 0.0;

  if (myid == 0)
    {
      zR.SetSize(fomSp_->zR.Size());
      BRsp = new CAROM::Matrix(fomSp->zR.Size(), rrdim, false);
      BWsp = new CAROM::Matrix(fomSp->zW.Size(), rwdim, false);
    }
  
  V_R.transposeMult(*U_R, VTU_R);

  GatherDistributedMatrixRows(V_R, V_W, rrdim, rwdim, fom->fespace_R.GetVSize(), fom->fespace_R, fom->fespace_W, st2sp, sprows, all_sprows, *BRsp, *BWsp);
  
  // Compute BR = V_W^t B V_R and CR = V_W^t C V_W, and store them throughout the simulation.

  BR = new CAROM::Matrix(rwdim, rrdim, false);
  CR = new CAROM::Matrix(rwdim, rwdim, false);
  Compute_CtAB(fom->Bmat, V_R, V_W, BR);
  Compute_CtAB(fom->Cmat, V_W, V_W, CR);

  // The ROM residual is
  // [ V_{R,s}^{-1} M(a(Pst V_W u)) Pst V_R v + V_R^t B^T V_W u ]
  // [ V_W^t C V_W du_dt - V_W^t B V_R v - V_W^t f ]
  // or, with [v, u] = [V_R yR, V_W yW],
  // [ V_{R,s}^{-1} M(a(Pst V_W yW)) Pst V_R yR + BR^T yW ]
  // [ CR dyW_dt - BR yR - V_W^t f ]
  // The Jacobian with respect to [dyR_dt, dyW_dt], with [yR, yW] = [yR0, yW0] + dt * [dyR_dt, dyW_dt], is
  // [ dt V_{R,s}^{-1} M(a'(Pst V_W yW)) Pst V_R  dt BR^T ]
  // [                 -dt BR                        CR   ]

  if (myid == 0)
    {
      const double linear_solver_rel_tol = 1.0e-14;
  
      J_gmres = new GMRESSolver;
      J_gmres->SetRelTol(linear_solver_rel_tol);
      J_gmres->SetAbsTol(0.0);
      J_gmres->SetMaxIter(1000);
      J_gmres->SetPrintLevel(1);
  
      newton_solver.iterative_mode = true;
      newton_solver.SetSolver(*J_gmres);
      newton_solver.SetOperator(*this);
      newton_solver.SetPrintLevel(1);
      newton_solver.SetRelTol(newton_rel_tol);
      newton_solver.SetAbsTol(newton_abs_tol);
      newton_solver.SetMaxIter(newton_iter);

      const int spdim = fomSp->Height();

      usp_librom = new CAROM::Vector(spdim, false);
      usp = new Vector(&((*usp_librom)(0)), spdim);

      // Define sub-vectors of usp. 
      usp_R = new Vector(usp->GetData(), fomSp->zR.Size());
      usp_W = new Vector(usp->GetData() + fomSp->zR.Size(), fomSp->zW.Size());

      usp_R_librom = new CAROM::Vector(usp_R->GetData(), usp_R->Size(), false, false);
      usp_W_librom = new CAROM::Vector(usp_W->GetData(), usp_W->Size(), false, false);
  
      MFEM_VERIFY(nsamp_R == s2sp.size(), "");
    }
  
  hyperreduce = true;
  sourceFOM = false;

  if (!hyperreduce || sourceFOM)
    {
      const int fdim = fom->Height();

      ufom_librom = new CAROM::Vector(fdim, false);
      ufom = new Vector(&((*ufom_librom)(0)), fdim);

      // Define sub-vectors of ufom. 
      ufom_R = new Vector(ufom->GetData(), fom->zR.Size());
      ufom_W = new Vector(ufom->GetData() + fom->zR.Size(), fom->zW.Size());

      ufom_R_librom = new CAROM::Vector(ufom_R->GetData(), ufom_R->Size(), false, false);
      ufom_W_librom = new CAROM::Vector(ufom_W->GetData(), ufom_W->Size(), false, false);

      zfomR.SetSize(fom->zR.Size());
      zfomR_librom = new CAROM::Vector(zfomR.GetData(), zfomR.Size(), false, false);

      zfomW.SetSize(fom->zW.Size());
    }

  if (hyperreduce_source)
    Compute_CtAB(fom->Cmat, *S, V_W, &VTCS_W);
}

RomOperator::~RomOperator()
{
  delete BR;
  delete CR;
}

void RomOperator::Mult_Hyperreduced(const Vector &dy_dt, Vector &res) const
{
  MFEM_VERIFY(dy_dt.Size() == rrdim + rwdim && res.Size() == rrdim + rwdim, "");
  
  Vector y(y0);
  y.Add(current_dt, dy_dt);

  // Evaluate the ROM residual:
  // [ V_R^T U_R U_{R,s}^{-1} M(a(Pst V_W yW)) Pst V_R yR + BR^T yW ]
  // [ CR dyW_dt - BR yR - V_W^t C f ]

  CAROM::Vector y_librom(y.GetData(), y.Size(), false, false);
  CAROM::Vector yR_librom(y.GetData(), rrdim, false, false);
  CAROM::Vector yW_librom(y.GetData() + rrdim, rwdim, false, false);

  CAROM::Vector resR_librom(res.GetData(), rrdim, false, false);
  CAROM::Vector resW_librom(res.GetData() + rrdim, rwdim, false, false);

  CAROM::Vector dyW_dt_librom(dy_dt.GetData() + rrdim, rwdim, false, false);
  
  // 1. Lift u_s+ = B_s+ y
  BRsp->mult(yR_librom, *usp_R_librom);
  BWsp->mult(yW_librom, *usp_W_librom);
  
  fomSp->SetParameters(*usp);

  fomSp->Mmat->Mult(*usp_R, zR);  // M(a(Pst V_W yW)) Pst V_R yR

  // Select entries out of zR.
  for (int i=0; i<nsamp_R; ++i)
    zN(i) = zR[s2sp[i]];

  // Note that it would be better to just store VTU_R * Vsinv, but these are small matrices.
#ifdef USE_GNAT
  Vsinv->transposeMult(zN, zY);
#else
  Vsinv->mult(zN, zY);
#endif

  BR->transposeMult(yW_librom, resR_librom);
  VTU_R.multPlus(resR_librom, zY, 1.0);

  // Apply V_W^t C to fsp

  if (sourceFOM)
    {
      fom->GetSource(zfomW);
      zfomW.Neg();

      fom->Cmat->Mult(zfomW, *ufom_W);
  
      V_W.transposeMult(*ufom_W_librom, resW_librom);

      CR->multPlus(resW_librom, dyW_dt_librom, 1.0);
      BR->multPlus(resW_librom, yR_librom, -1.0);
    }
  else
    {
      CR->mult(dyW_dt_librom, resW_librom);
      BR->multPlus(resW_librom, yR_librom, -1.0);

      fomSp->GetSource(fomSp->zW);

      if (hyperreduce_source)
	{
	  // Select entries
	  for (int i=0; i<nsamp_S; ++i)
	    zT(i) = fomSp->zW(s2sp_S[i]);

#ifdef USE_GNAT
	  Ssinv->transposeMult(zT, zS);
#else
	  Ssinv->mult(zT, zS);
#endif

	  // Multiply by the f-basis, followed by C, followed by V_W^T. This is stored in VTCS_W = V_W^T CS.
	  VTCS_W.multPlus(resW_librom, zS, -1.0);
	}
      else
	{
	  fomSp->Cmat->Mult(fomSp->zW, *usp_W);

	  const int nRsp = fomSp->zR.Size();
	  const int nWsp = fomSp->zW.Size();
	  for (int i=0; i<rwdim; ++i)
	    for (int j=0; j<nWsp; ++j)
	      res[rrdim + i] -= (*BWsp)(j, i) * (*usp_W)[j];
	}
    }
}

void RomOperator::Mult_FullOrder(const Vector &dy_dt, Vector &res) const
{
  MFEM_VERIFY(dy_dt.Size() == rrdim + rwdim && res.Size() == rrdim + rwdim, "");
  
  Vector y(y0);
  y.Add(current_dt, dy_dt);

  // Evaluate the unreduced ROM residual:
  // [ V_R^T M(a(V_W yW)) V_R yR + BR^T yW ]
  // [ CR dyW_dt - BR yR - V_W^t Cf ]
  
  CAROM::Vector y_librom(y.GetData(), y.Size(), false, false);
  CAROM::Vector yR_librom(y.GetData(), rrdim, false, false);
  CAROM::Vector yW_librom(y.GetData() + rrdim, rwdim, false, false);

  CAROM::Vector resR_librom(res.GetData(), rrdim, false, false);
  CAROM::Vector resW_librom(res.GetData() + rrdim, rwdim, false, false);

  CAROM::Vector dyW_dt_librom(dy_dt.GetData() + rrdim, rwdim, false, false);
  
  // 1. Lift u_fom = [V_R^T V_W^T]^T y
  V_R.mult(yR_librom, *ufom_R_librom);
  V_W.mult(yW_librom, *ufom_W_librom);
  
  fom->SetParameters(*ufom);

  fom->Mmat->Mult(*ufom_R, zfomR);  // M(a(V_W yW)) V_R yR
  V_R.transposeMult(*zfomR_librom, VtzR);  // V_R^T M(a(V_W yW)) V_R yR

  BR->transposeMult(yW_librom, resR_librom);
  resR_librom.plusEqAx(1.0, VtzR);

  // Apply V_W^t C to f
  fom->GetSource(zfomW);
  zfomW.Neg();

  fom->Cmat->Mult(zfomW, *ufom_W);
  
  V_W.transposeMult(*ufom_W_librom, resW_librom);
  
  CR->multPlus(resW_librom, dyW_dt_librom, 1.0);
  BR->multPlus(resW_librom, yR_librom, -1.0);
}

void RomOperator::Mult(const Vector &dy_dt, Vector &res) const
{
  if (hyperreduce)
    Mult_Hyperreduced(dy_dt, res);
  else
    Mult_FullOrder(dy_dt, res);
}

void RomOperator::ImplicitSolve(const double dt, const Vector &y, Vector &dy_dt)
{
  y0 = y;

  current_dt = dt;
  fomSp->SetTime(GetTime());
  fomSp->current_dt = dt;

  if (!hyperreduce || sourceFOM)
    {
      fom->SetTime(GetTime());
      fom->current_dt = dt;
    }
  
  // Set the initial guess for du_dt, to be used by newton_solver.
  //du_dt = 0.0;
  dy_dt = dydt_prev;

  Vector zero; // empty vector is interpreted as zero r.h.s. by NewtonSolver
  newton_solver.Mult(zero, dy_dt);

  // MFEM_VERIFY(newton_solver.GetConverged(), "Newton solver did not converge.");
  if (newton_solver.GetConverged())
    dydt_prev = dy_dt;
  else
    {
      dy_dt = 0.0;  // Zero update in SDIRK Step() function.
      //newtonFailure = true;
      MFEM_VERIFY(false, "ROM Newton convergence failure!");
    }
}

// Debugging tool for checking Jacobian
void RomOperator::PrintFDJacobian(const Vector &u) const
{
  const int N = u.Size();
  Vector up(N);
  Vector r(N);
  Vector r0(N);

  const double d = 1.0e-8;
  
  DenseMatrix JFD(N);

  Mult(u, r0);

  for (int j=0; j<N; ++j)
    {
      up = u;
      up[j] += d;

      Mult(up, r);

      r -= r0;
      r /= d;

      for (int i=0; i<N; ++i)
	JFD(i,j) = r[i];
    }

  JFD.Print(cout);
}

Operator &RomOperator::GetGradient(const Vector &u) const
{
  // The Jacobian with respect to [dyR_dt, dyW_dt], with [yR, yW] = [yR0, yW0] + dt * [dyR_dt, dyW_dt], is
  // [ dt V_{R,s}^{-1} M(a'(Pst V_W yW)) Pst V_R  dt BR^T ]
  // [                 -dt BR                        CR   ]

  // Compute JR = V_{R,s}^{-1} M(a'(Pst V_W yW)) Pst V_R, assuming M(a'(Pst V_W yW)) is already stored in fomSp->Mprimemat,
  // which is computed in fomSp->SetParameters, which was called by RomOperator::Mult, which was called by newton_solver
  // before this call to GetGradient. Note that V_R restricted to the sample matrix is already stored in Bsp.

  CAROM::Vector r(nldim, false);
  CAROM::Vector c(rrdim, false);
  CAROM::Vector z(nsamp_R, false);
  
  for (int i=0; i<rrdim; ++i)
    {
      if (hyperreduce)
	{
	  // Compute the i-th column of M(a'(Pst V_W yW)) Pst V_R.
	  for (int j=0; j<usp_R->Size(); ++j)
	    (*usp_R)[j] = (*BRsp)(j,i);
      
	  fomSp->Mprimemat.Mult(*usp_R, zR);

	  for (int j=0; j<nsamp_R; ++j)
	    z(j) = zR[s2sp[j]];

	  // Note that it would be better to just store VTU_R * Vsinv, but these are small matrices. 
      
#ifdef USE_GNAT
	  Vsinv->transposeMult(z, r);
#else
	  Vsinv->mult(z, r);
#endif

	  VTU_R.mult(r, c);
	}
      else
	{
	  // Compute the i-th column of V_R^T M(a'(V_W yW)) V_R.
	  for (int j=0; j<ufom_R->Size(); ++j)
	    (*ufom_R)[j] = V_R(j,i);
      
	  fom->Mprimemat.Mult(*ufom_R, zfomR);
	  V_R.transposeMult(*zfomR_librom, c);  // V_R^T M(a'(V_W yW)) V_R(:,i)
	}
      
      for (int j=0; j<rrdim; ++j)
	J(j, i) = c(j);  // This already includes a factor of current_dt, from Mprimemat.
	  
      for (int j=0; j<rwdim; ++j)
	{
	  J(rrdim + j, i) = -current_dt * (*BR)(j,i);
	  J(i, rrdim + j) = current_dt * (*BR)(j,i);
	}
    }

  for (int i=0; i<rwdim; ++i)
    {
      for (int j=0; j<rwdim; ++j)
	{
	  J(rrdim + j, rrdim + i) = (*CR)(j,i);
	}
    }

  // TODO: define Jacobian block-wise rather than entry-wise?
  /*
  gradient->SetBlock(0, 0, JR, current_dt);
  gradient->SetBlock(0, 1, BRT, current_dt);
  gradient->SetBlock(1, 0, BR, -current_dt);
  gradient->SetBlock(1, 1, CR);
  */

  // PrintFDJacobian(u);

  return J;
}

double InitialTemperature(const Vector &x)
{
  if (problem == INIT_STEP)
    {
      if (0.5 - step_half < x[0] && x[0] < 0.5 + step_half && 0.5 - step_half < x[1] && x[1] < 0.5 + step_half)
	return 1.0;
      else
	return 0.0;
    }
  else
    return 0.0;
}

double ExactSolution(const Vector &x, const double t)
{
  const double pi = acos(-1.0);
  const double pi2 = 2.0 * acos(-1.0);
  return sin(4.0 * pi * t) * sin(pi2 * x[0]) * sin(pi2 * x[1]);
}

double SourceFunction_linear(const Vector &x, const double t)
{
  // du/dt + div(v) = f, grad u = -a(u) v
  // u(x,y) = sin(2 pi x) sin(2 pi y) sin(4 pi t)
  // a(u) = 1
  // Set cx = cos(2 pi x), sy = sin(2 pi y), etc. Then v = -2 pi [cx sy st, sx cy st]
  // div(v) = 8 pi^2 sx sy st
  
  // du/dt + div(v) = 4 pi sin(2 pi x) sin(2 pi y) cos(4 pi t) + 8 pi^2 sx sy st

  //return 0.0;

  const double pi = acos(-1.0);
  const double pi2 = 2.0 * pi;

  const double sx = sin(pi2 * x[0]);
  const double sy = sin(pi2 * x[1]);

  const double st = sin(4.0 * pi * t);
  const double ct = cos(4.0 * pi * t);

  return (4.0 * pi * ct * sx * sy) + (st * 8.0 * pi * pi * sx * sy);
}

// a(u) = c + u, where c = diffusion_c
double SourceFunction_cpu(const Vector &x, const double t)
{
  // du/dt + div(v) = f, grad u = -a(u) v
  // u(x,y) = sin(2 pi x) sin(2 pi y) sin(4 pi t)
  // a(u) = c + sin(2 pi x) sin(2 pi y) sin(4 pi t)
  // Set cx = cos(2 pi x), sy = sin(2 pi y), etc. Then v = -2 pi st/(2+sx sy st) * [cx sy, sx cy] = -2 pi st/a * [cx sy, sx cy]

  // div(v) = 4 pi^2 cx^2 sy^2 st^2 / a^2 + 4 pi^2 st / a * sx sy         (dv_x/dx)
  //            + 4 pi^2 sx^2 cy^2 st^2 / a^2 + 4 pi^2 st / a * sx sy     (dv_y/dy)
  //        = 4 pi^2 st / a * [(cx sy)^2 st / a + 2 sx sy + (sx cy)^2 st / a]
  
  // du/dt + div(v) = 4 pi sx sy ct + 4 pi^2 st / a * [(cx sy)^2 st / a + 2 sx sy + (sx cy)^2 st / a]

  //return 0.0;

  const double pi = acos(-1.0);
  const double pi2 = 2.0 * acos(-1.0);

  const double cx = cos(pi2 * x[0]);
  const double cy = cos(pi2 * x[1]);

  const double sx = sin(pi2 * x[0]);
  const double sy = sin(pi2 * x[1]);

  const double st = sin(4.0 * pi * t);
  const double ct = cos(4.0 * pi * t);

  const double a = diffusion_c + (sx * sy * st);

  return (4.0 * pi * ct * sx * sy) + (st * 4.0 * pi * pi * ((cx*sy*cx*sy*st/a) + (2.0*sx*sy) + (sx*cy*sx*cy*st/a)) / a);
}

double SourceFunction(const Vector &x, const double t)
{
  if (nonlinear_problem)
    {
      if (problem == INIT_STEP)
	return 0.0;
      else
	return SourceFunction_cpu(x, t);
    }
  else
    return SourceFunction_linear(x, t);
}

double NonlinearCoefficient(const double u)
{
  if (nonlinear_problem)
    return diffusion_c + u;
  else
    return 1.0;
}

double NonlinearCoefficientDerivative(const double u)
{
  if (nonlinear_problem)
    return 1.0;
  else
    return 0.0;
}
