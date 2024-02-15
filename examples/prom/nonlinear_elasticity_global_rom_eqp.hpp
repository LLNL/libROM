/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Functions used by nonlinear_elastic_global_rom.cpp with EQP.
#include "mfem/Utilities.hpp"
#include <cmath>
using namespace mfem;
using namespace std;

class HyperelasticOperator : public TimeDependentOperator
{

protected:
    ParBilinearForm *M, *S;

    CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
    HypreSmoother M_prec; // Preconditioner for the mass matrix M

public:
    HyperelasticOperator(ParFiniteElementSpace &f, Array<int> &ess_tdof_list_,
                         double visc, double mu, double K);

    /// Compute the right-hand side of the ODE system.
    virtual void Mult(const Vector &vx, Vector &dvx_dt) const;

    double ElasticEnergy(const ParGridFunction &x) const;
    double KineticEnergy(const ParGridFunction &v) const;
    void GetElasticEnergyDensity(const ParGridFunction &x,
                                 ParGridFunction &w) const;

    mutable Vector H_sp;
    mutable Vector dvxdt_sp;

    ParFiniteElementSpace &fespace;
    double viscosity;
    Array<int> ess_tdof_list;
    ParNonlinearForm *H;
    HyperelasticModel *model;
    mutable Vector z;     // auxiliary vector
    mutable Vector z2;    // auxiliary vector
    HypreParMatrix *Mmat; // Mass matrix from ParallelAssemble()
    HypreParMatrix Smat;

    virtual ~HyperelasticOperator();
};


struct ElemMatrices
{
    DenseMatrix DSh;
    DenseMatrix DS;
    DenseMatrix Jrt;
    DenseMatrix Jpt;
    DenseMatrix P;
    DenseMatrix P_f;
    DenseMatrix PMatI;
    DenseMatrix PMatO;
    Vector elvect;
    // Constructor for matrices struct
    ElemMatrices(int dof, int dim) : DSh(dof, dim),
        DS(dof, dim),
        Jrt(dim),
        Jpt(dim),
        P(dim),
        P_f(dim),
        PMatI(),
        PMatO(),
        elvect(dof * dim)
    {
        // Set dimensions for PMatI and PMatO
        PMatO.UseExternalData(elvect.GetData(), dof, dim);
    }
};

// Compute coefficients of the reduced integrator with respect to inputs Q and x
// in HyperelasticNLFIntegrator_ComputeReducedEQP.
void GetEQPCoefficients_HyperelasticNLFIntegrator(ParFiniteElementSpace *fesR,
        std::vector<double> const &rw, std::vector<int> const &qp,
        const IntegrationRule *ir, NeoHookeanModel *model,
        CAROM::Matrix const &V_v, const int rank, Vector &coef, Vector &DS_coef,
        ElemMatrices *em);

// Perform hyperreduction with EQP
void HyperelasticNLFIntegrator_ComputeReducedEQP(ParFiniteElementSpace *fesR,
        std::vector<double> const &rw, std::vector<int> const &qp,
        const IntegrationRule *ir, NeoHookeanModel *model, const Vector *x0,
        CAROM::Matrix const &V_v, CAROM::Vector const &x, const int rank, Vector &res,
        ElemMatrices *em, const CAROM::Matrix eqp_lifting,
        const std::vector<int> eqp_liftDOFs,CAROM::Vector eqp_lifted);

// Optimized EQP hyperreduction routine with preallocated arrays
void HyperelasticNLFIntegrator_ComputeReducedEQP_Fast(ParFiniteElementSpace
        *fesR,std::vector<double> const &rw, std::vector<int> const &qp,
        const IntegrationRule *ir, NeoHookeanModel *model, const Vector *x0,
        const int rvdim, CAROM::Vector const &x, Vector const &coef,
        Vector const &DS_coef, const int rank, Vector &res, ElemMatrices *em,
        const CAROM::Matrix eqp_lifting, const std::vector<int> eqp_liftDOFs,
        CAROM::Vector eqp_lifted);

// Compute a row in the G matrix which corresponds to a given FE element
void ComputeElementRowOfG(const IntegrationRule *ir, Array<int> const &vdofs,
                          Vector const &ve_j, NeoHookeanModel *model, Vector const &elfun,
                          FiniteElement const &fe, ElementTransformation &Trans, Vector &r, const int dof,
                          const int dim,
                          ElemMatrices &em);

void SolveNNLS(const int rank, const double nnls_tol, const int maxNNLSnnz,
               CAROM::Vector const &w, CAROM::Matrix &Gt,
               CAROM::Vector &sol);

// Compute EQP solution from constraints on snapshots.
void SetupEQP_snapshots(const IntegrationRule *ir0, const int rank,
                        ParFiniteElementSpace *fespace_X,
                        const int nsets, const CAROM::Matrix *BV,
                        const CAROM::Matrix *BX_snapshots, const Vector x0,
                        const bool precondition, const double nnls_tol,
                        const int maxNNLSnnz, NeoHookeanModel *model,
                        CAROM::Vector &sol, CAROM::Vector *window_ids, const int snap_step);

void WriteMeshEQP(ParMesh *pmesh, const int myid, const int nqe,
                  CAROM::Vector const &eqpSol);

// Function to compute the indices at which to switch time windows
void get_window_ids(int n_step, int n_window, CAROM::Vector *ids);

// Helper function to check if a file exists
bool fileExists(const std::string &filename);

// Load a EQP solution
void get_EQPsol(const int current_window, CAROM::Vector *load_eqpsol);
