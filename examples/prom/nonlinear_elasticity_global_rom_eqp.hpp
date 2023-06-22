/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Functions used by mixed_nonlinear_diffusion.cpp with EQP.
#include "mfem/Utilities.hpp"
// #include "fem/nonlininteg.hpp"

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

// Compute coefficients of the reduced integrator with respect to inputs Q and x
// in HyperelasticNLFIntegrator_ComputeReducedEQP.
void GetEQPCoefficients_HyperelasticNLFIntegrator(ParFiniteElementSpace *fesR,
                                                  std::vector<double> const &rw, std::vector<int> const &qp,
                                                  const IntegrationRule *ir,
                                                  CAROM::Matrix const &V, Vector &res)
{
    MFEM_ABORT("TODO");
}

void HyperelasticNLFIntegrator_ComputeReducedEQP(ParFiniteElementSpace *fesR,
                                                 std::vector<double> const &rw, std::vector<int> const &qp,
                                                 const IntegrationRule *ir, NeoHookeanModel *model, const Vector *x0, 
                                                 CAROM::Matrix const &V_x,CAROM::Matrix const &V_v, CAROM::Vector const &x, const int rank, Vector &res)
{
    const int rxdim = V_x.numColumns();
    const int rvdim = V_v.numColumns();
    const int fomdim = V_v.numRows();
    MFEM_VERIFY(rw.size() == qp.size(), "");
    MFEM_VERIFY(x.dim() == rxdim, "");
    MFEM_VERIFY(V_x.numRows() == fesR->GetTrueVSize(), "");

    MFEM_VERIFY(rank == 0,
                "TODO: generalize to parallel. This uses full dofs in V, which has true dofs");

    const int nqe = ir->GetWeights().Size();

    ElementTransformation *eltrans;
    DofTransformation *doftrans;
    const FiniteElement *fe = NULL;
    Array<int> vdofs;

    res.SetSize(rxdim);
    res = 0.0;

    int eprev = -1;
    int dof = 0;
    int dim = 0;

    // Get prolongation matrix
    const Operator *P = fesR->GetProlongationMatrix();

    // Vectors to be prolongated
    CAROM::Vector *Vx_librom_temp = new CAROM::Vector(fomdim, true);
    Vector *Vx_temp = new Vector(&((*Vx_librom_temp)(0)), fomdim);
    Vector Vx(fomdim);
    Vector vj(fomdim);

    // Prolongated vectors
    Vector p_Vx(P->Height());
    Vector p_vj(P->Height());

    // Element vectors
    Vector Vx_e;
    Vector vj_e;

    // Lift x, add x0 and prolongate result
    V_x.mult(x, Vx_librom_temp);
    add(*Vx_temp, *x0, Vx);
    P->Mult(Vx, p_Vx);

    // For every basis vector
    for (int j = 0; j < rvdim; ++j)

    {
        // Get basis vector and prolongate
        for (int k = 0; k < V_v.numRows(); ++k)
            vj[k] = V_v(k, j);
        P->Mult(vj, p_vj);
        res[j] = 0.0;

        eprev = -1;

        // For every quadrature weight
        for (int i = 0; i < rw.size(); ++i) // NOTE: i < 9
        {
            const int e = qp[i] / nqe; // Element index
            // Local (element) index of the quadrature point
            const int qpi = qp[i] - (e * nqe);  
            const IntegrationPoint &ip = ir->IntPoint(qpi);

            if (e != eprev) // Update element transformation
            {
                doftrans = fesR->GetElementVDofs(e, vdofs);
                fe = fesR->GetFE(e);
                eltrans = fesR->GetElementTransformation(e);

                dof = fe->GetDof(); // Get number of dofs in element
                dim = fe->GetDim();

                if (doftrans)
                {
                    MFEM_ABORT("TODO");
                }

                // Get element vectors
                p_Vx.GetSubVector(vdofs, Vx_e);
                p_vj.GetSubVector(vdofs, vj_e);
                eprev = e;

            }

            // Integration at ip

            // Initialize nonlinear operator matrices (there is probably a better way)
            DenseMatrix DSh(dof, dim);
            DenseMatrix DS(dof, dim);
            DenseMatrix Jrt(dim);
            DenseMatrix Jpt(dim);
            DenseMatrix P_f(dim);
            DenseMatrix PMatI; // Extract element dofs
            PMatI.UseExternalData(Vx_e.GetData(), dof, dim);
            DenseMatrix PMatO;
            Vector elvect(dof * dim);
            PMatO.UseExternalData(elvect.GetData(), dof, dim);

            elvect = 0.0;

            // Set integration point in the element transformation
            eltrans->SetIntPoint(&ip);
            model->SetTransformation(*eltrans);

            // Get the transformation weight
            double t = eltrans->Weight();

            // Compute action of nonlinear operator
            CalcInverse(eltrans->Jacobian(), Jrt);
            fe->CalcDShape(ip, DSh);
            Mult(DSh, Jrt, DS);
            MultAtB(PMatI, DS, Jpt);
            model->EvalP(Jpt, P_f);
            P_f *= (t * rw[i]); // NB: Not by ip.weight //Or is it rw[qpi]?
            AddMultABt(DS, P_f, PMatO);

            // Calculate r[i] = ve_j^T * elvect
            for (int k = 0; k < elvect.Size(); k++)
            {
                res[j] += vj_e[k] * elvect[k];
            }
        }
    }
}

/* TODO if time...
void HyperelasticNLFIntegrator_ComputeReducedEQP_Fast(ParFiniteElementSpace *fesR,
        std::vector<int> const& qp, const IntegrationRule *ir,
        Coefficient *Q, CAROM::Vector const& x,
        Vector const& coef, Vector & res)
{

}
*/

void ComputeElementRowOfG(const IntegrationRule *ir, Array<int> const &vdofs,
                          Vector const &ve_j, NeoHookeanModel *model, Vector const &elfun,
                          FiniteElement const &fe, ElementTransformation &Trans, Vector &r)
{
    MFEM_VERIFY(r.Size() == ir->GetNPoints(), "");
    int dof = fe.GetDof(); // Get number of dofs in element
    int dim = fe.GetDim();

    // Initialize nonlinear operator matrices (there is probably a better way)
    DenseMatrix DSh(dof, dim);
    DenseMatrix DS(dof, dim);
    DenseMatrix Jrt(dim);
    DenseMatrix Jpt(dim);
    DenseMatrix P(dim);
    DenseMatrix PMatI; // Extract element dofs
    PMatI.UseExternalData(elfun.GetData(), dof, dim);
    DenseMatrix PMatO;
    Vector elvect(dof * dim);
    PMatO.UseExternalData(elvect.GetData(), dof, dim);

    model->SetTransformation(Trans);

    // For each integration point
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        elvect = 0.0;
        // Get integration point
        const IntegrationPoint &ip = ir->IntPoint(i);

        // Set integration point in the element transformation
        Trans.SetIntPoint(&ip);

        // Get the transformation weight
        double t = Trans.Weight();

        // Compute action of nonlinear operator
        CalcInverse(Trans.Jacobian(), Jrt);
        fe.CalcDShape(ip, DSh);
        Mult(DSh, Jrt, DS);
        MultAtB(PMatI, DS, Jpt);
        model->EvalP(Jpt, P);
        P *= t; // NB: Not by ip.weight
        AddMultABt(DS, P, PMatO);

        r[i] = 0.0;

        // Calculate r[i] = ve_j^T * elvect
        for (int k = 0; k < elvect.Size(); k++)
        {
            r[i] += ve_j[k] * elvect[k];
        }
    }
}

void SolveNNLS(const int rank, const double nnls_tol, const int maxNNLSnnz,
               CAROM::Vector const &w, CAROM::Matrix &Gt,
               CAROM::Vector &sol)
{
    CAROM::NNLSSolver nnls(nnls_tol, 0, maxNNLSnnz, 2);

    CAROM::Vector rhs_ub(Gt.numColumns(), false);
    // G.mult(w, rhs_ub);  // rhs = Gw
    // rhs = Gw. Note that by using Gt and multTranspose, we do parallel communication.
    Gt.transposeMult(w, rhs_ub);

    CAROM::Vector rhs_lb(rhs_ub);
    CAROM::Vector rhs_Gw(rhs_ub);

    const double delta = 1.0e-11;
    for (int i = 0; i < rhs_ub.dim(); ++i)
    {
        rhs_lb(i) -= delta;
        rhs_ub(i) += delta;
    }

    // nnls.normalize_constraints(Gt, rhs_lb, rhs_ub);
    nnls.solve_parallel_with_scalapack(Gt, rhs_lb, rhs_ub, sol);

    int nnz = 0;
    for (int i = 0; i < sol.dim(); ++i)
    {
        if (sol(i) != 0.0)
        {
            nnz++;
        }
    }

    cout << rank << ": Number of nonzeros in NNLS solution: " << nnz
         << ", out of " << sol.dim() << endl;

    MPI_Allreduce(MPI_IN_PLACE, &nnz, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0)
        cout << "Global number of nonzeros in NNLS solution: " << nnz << endl;

    // Check residual of NNLS solution
    CAROM::Vector res(Gt.numColumns(), false);
    Gt.transposeMult(sol, res);

    const double normGsol = res.norm();
    const double normRHS = rhs_Gw.norm();

    res -= rhs_Gw;
    const double relNorm = res.norm() / std::max(normGsol, normRHS);
    cout << rank << ": relative residual norm for NNLS solution of Gs = Gw: " << relNorm << endl;
}

// Compute EQP solution from constraints on snapshots.
void SetupEQP_snapshots(const IntegrationRule *ir0, const int rank,
                        ParFiniteElementSpace *fespace_X,
                        const int nsets, const CAROM::Matrix *BV,
                        const CAROM::Matrix *BX_snapshots,
                        const bool precondition, const double nnls_tol,
                        const int maxNNLSnnz, NeoHookeanModel *model,
                        CAROM::Vector &sol)
{
    const int nqe = ir0->GetNPoints();
    const int ne = fespace_X->GetNE();
    const int NB = BV->numColumns();
    const int NQ = ne * nqe;
    const int nsnap = BX_snapshots->numColumns();

    MFEM_VERIFY(nsnap == BX_snapshots->numColumns() ||
                    nsnap + nsets == BX_snapshots->numColumns(), // Q: nsets?
                "");
    MFEM_VERIFY(BV->numRows() == BX_snapshots->numRows(), "");
    MFEM_VERIFY(BV->numRows() == fespace_X->GetTrueVSize(), "");

    const bool skipFirstW = (nsnap + nsets == BX_snapshots->numColumns());

    // Compute G of size (NB * nsnap) x NQ, but only store its transpose Gt.
    CAROM::Matrix Gt(NQ, NB * nsnap, true);

    // For 0 <= j < NB, 0 <= i < nsnap, 0 <= e < ne, 0 <= m < nqe,
    // G(j + (i*NB), (e*nqe) + m)
    // is the coefficient of v_j^T M(p_i) V v_i at point m of element e,
    // with respect to the integration rule weight at that point,
    // where the "exact" quadrature solution is ir0->GetWeights().

    Vector x_i(BX_snapshots->numRows());
    Vector v_j(BV->numRows());

    Vector r(nqe);

    int skip = 0;
    const int nsnapPerSet = nsnap / nsets;
    if (skipFirstW)
    {
        MFEM_VERIFY(nsets * nsnapPerSet == nsnap, "");
        skip = 1;
    }

    // Get prolongation matrix
    const Operator *P = fespace_X->GetProlongationMatrix();
    if (!P)
    {
        MFEM_ABORT("P is null, generalize to serial case")
    }

    // For every snapshot
    for (int i = 0; i < nsnap; ++i)
    {
        // Set the sampled dofs from the snapshot matrix
        for (int j = 0; j < BX_snapshots->numRows(); ++j)
            x_i[j] = (*BX_snapshots)(j, i);

        // Get prolongated dofs
        Vector px_i;
        Vector elfun;

        px_i.SetSize(P->Height());
        P->Mult(x_i, px_i);

        if (skipFirstW && i > 0 && i % nsnapPerSet == 0)
            skip++;

        // For each basis vector
        for (int j = 0; j < NB; ++j)
        {

            // Get basis vector
            for (int k = 0; k < BV->numRows(); ++k)
                v_j[k] = (*BV)(k, j);

            // Get prolongated dofs
            Vector pv_j;
            Vector ve_j;

            pv_j.SetSize(P->Height());
            P->Mult(v_j, pv_j);

            // TODO: is it better to make the element loop the outer loop?
            // For each element
            for (int e = 0; e < ne; ++e)
            {
                // Get element and its dofs and transformation.
                Array<int> vdofs;
                DofTransformation *doftrans = fespace_X->GetElementVDofs(e, vdofs);
                const FiniteElement &fe = *fespace_X->GetFE(e);
                ElementTransformation *eltrans = fespace_X->GetElementTransformation(e);
                px_i.GetSubVector(vdofs, elfun);
                pv_j.GetSubVector(vdofs, ve_j);
                if (doftrans)
                {
                    MFEM_ABORT("Doftrans is true, make corresponding edits")
                }

                // Compute the row of G corresponding to element e, store in r
                ComputeElementRowOfG(ir0, vdofs, ve_j, model, elfun, fe, *eltrans, r);

                for (int m = 0; m < nqe; ++m)
                    Gt((e * nqe) + m, j + (i * NB)) = r[m];
            }
        }

        if (precondition)
        {
            MFEM_ABORT("TODO");
        }
    } // Loop (i) over snapshots

    Array<double> const &w_el = ir0->GetWeights();
    MFEM_VERIFY(w_el.Size() == nqe, "");

    CAROM::Vector w(ne * nqe, true);

    for (int i = 0; i < ne; ++i)
    {
        for (int j = 0; j < nqe; ++j)
            w((i * nqe) + j) = w_el[j];
    }
    bool stop = true;

    SolveNNLS(rank, nnls_tol, maxNNLSnnz, w, Gt, sol);
}

void WriteMeshEQP(ParMesh *pmesh, const int myid, const int nqe,
                  CAROM::Vector const &eqpSol)
{
    // Find the elements with quadrature points in eqpSol.
    std::set<int> elements;
    Vector elemCount(pmesh->GetNE());
    elemCount = 0.0;

    for (int i = 0; i < eqpSol.dim(); ++i)
    {
        if (eqpSol(i) > 1.0e-12)
        {
            const int e = i / nqe; // Element index
            elements.insert(e);
            elemCount[e] += 1.0;
        }
    }

    // Empty sets, since EQP on samples inside elements.
    std::set<int> faces, edges, vertices;
    CAROM::SampleVisualization(pmesh, elements, elements, faces, edges,
                               vertices, "EQPvis", &elemCount);
}
