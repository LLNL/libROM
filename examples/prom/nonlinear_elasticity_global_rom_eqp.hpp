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
    MFEM_ABORT("TODO")
}

void HyperelasticNLFIntegrator_ComputeReducedEQP(ParFiniteElementSpace *fesR,
                                                 std::vector<double> const &rw, std::vector<int> const &qp,
                                                 const IntegrationRule *ir,
                                                 CAROM::Matrix const &V, CAROM::Vector const &x, const int rank, Vector &res)
{
    MFEM_ABORT("TODO");
}
/*
{
    const int rdim = V.numColumns();
    MFEM_VERIFY(rw.size() == qp.size(), "");
    MFEM_VERIFY(x.dim() == rdim, "");
    MFEM_VERIFY(V.numRows() == fesR->GetTrueVSize(), "");

    MFEM_VERIFY(rank == 0,
                "TODO: generalize to parallel. This uses full dofs in V, which has true dofs");

    const int nqe = ir->GetWeights().Size();

    ElementTransformation *eltrans;
    DofTransformation *doftrans;
    const FiniteElement *fe = NULL;
    Array<int> vdofs;

    DenseMatrix trial_vshape;

    Vector Vx;
    Vector Vj;
    res.SetSize(rdim);
    res = 0.0;

    int eprev = -1;
    int dof = 0;
    int spaceDim = 0;

    // Note that the alternative version VectorFEMassIntegrator_ComputeReducedEQP_Fast
    // of this function is optimized by storing some intermediate computations.

    for (int i = 0; i < rw.size(); ++i)
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

            if (doftrans)
            {
                MFEM_ABORT("TODO");
            }

            dof = fe->GetDof();
            spaceDim = eltrans->GetSpaceDim();
            trial_vshape.SetSize(dof, spaceDim);
            Vx.SetSize(spaceDim);
            Vj.SetSize(spaceDim);

            // Initialize nonlinear operator matrices (there is probably a better way)
            int dim = fe->GetDim();
            DenseMatrix DSh(dof, dim);
            DenseMatrix DS(dof, dim);
            DenseMatrix Jrt(dim);
            DenseMatrix Jpt(dim);
            DenseMatrix P(dim);
            DenseMatrix PMatI;
            PMatI.UseExternalData(elfun.GetData(), dof, dim); // Q: How to replace?
            // PMatO.SetSize(dof, dim); // Q: should this be here?

            eprev = e;
        }

        // Integrate at the current point
        eltrans->SetIntPoint(&ip);
        fe->CalcVShape(*eltrans, trial_vshape);

        double w = eltrans->Weight() * rw[i]; // using rw[i] instead of ip.weight

        // Compute action of nonlinear operator
        CalcInverse(eltrans->Jacobian(), Jrt);
        fe.CalcDShape(ip, DSh);
        Mult(DSh, Jrt, DS);
        MultAtB(PMatI, DS, Jpt);
        model->EvalP(Jpt, P);
        // AddMultABt(DS, P, PMatO); //Q: Should this be here?

        // Lift Vx at ip.
        Vx = 0.0;
        for (int k = 0; k < dof; ++k)
        {
            const int dofk = (vdofs[k] >= 0) ? vdofs[k] : -1 - vdofs[k];

            double Vx_k = 0.0;
            for (int j = 0; j < rdim; ++j)
            {
                Vx_k += V(dofk, j) * x(j);
            }

            if (vdofs[k] < 0) // Q:why this?
                Vx_k = -Vx_k;

            for (int j = 0; j < spaceDim; ++j)
                Vx[j] += Vx_k * trial_vshape(k, j);
        }

        Vj = 0.0;
        for (int k = 0; k < spaceDim; ++k)
        {
            double Vjk = 0.0;
            for (int l = 0; l < dof; ++l)
            {
                // Q: Should these ones be the same?
                const int dofl = (vdofs[l] >= 0) ? vdofs[l] : -1 - vdofs[l];
                const double s = (vdofs[l] >= 0) ? 1.0 : -1.0;
                Vjk += s * V(dofl, j) * trial_vshape(l, k);
            }

            Vj[k] = Vjk;
        }

        // Calculate r[i] = Vj^T * P * Vx
        // Perform the vector-matrix-vector multiplication: a^T * B * c
        Vector temp(dim);
        P.Mult(Vx, temp);
        w *= Vj * temp;

        // Perform final scaling
        for (int j = 0; j < rdim; ++j)
        {
            res[j] += w
        }
    }
}
*/
/*
void HyperelasticNLFIntegrator_ComputeReducedEQP(ParFiniteElementSpace *fesR,
                                                 std::vector<double> const &rw, std::vector<int> const &qp,
                                                 const IntegrationRule *ir, Coefficient *Q,
                                                 CAROM::Matrix const &V, CAROM::Vector const &x, const int rank, Vector &res)
{

    const int rdim = V.numColumns();
    MFEM_VERIFY(V.numRows() == fesR->GetTrueVSize(), "");
    MFEM_VERIFY(rw.size() == qp.size(), "");

    const int ne = fesR->GetNE();
    const int nqe = ir->GetWeights().Size();

    ElementTransformation *eltrans;
    DofTransformation *doftrans;
    const FiniteElement *fe = NULL;
    Array<int> vdofs;

    DenseMatrix trial_vshape;

    Vector Vx;
    res.SetSize(rdim * rdim * rw.size());
    res = 0.0;

    // For the parallel case, we must get all DOFs of V on sampled elements.
    CAROM::Matrix Vs;

    // Since V only has rows for true DOFs, we use a ParGridFunction to get all DOFs.
    int eprev = -1;
    int dof = 0;
    int elemCount = 0;

    // First, find all sampled elements.
    for (int i = 0; i < rw.size(); ++i)
    {
        const int e = qp[i] / nqe; // Element index

        if (e != eprev) // Update element transformation
        {
            doftrans = fesR->GetElementVDofs(e, vdofs);
            if (dof > 0)
            {
                MFEM_VERIFY(dof == vdofs.Size(), "All elements must have same DOF size");
            }
            dof = vdofs.Size();
            eprev = e;
            elemCount++;
        }
    }

    // Now set Vs.
    // Note that, ideally, these FOM data structures and operations should be
    // avoided, but this only done in setup.
    Vs.setSize(elemCount * dof, rdim);
    ParGridFunction v_gf(fesR);
    Vector vtrue(fesR->GetTrueVSize());

    for (int j = 0; j < rdim; ++j)
    {
        eprev = -1;
        elemCount = 0;

        for (int i = 0; i < vtrue.Size(); ++i)
            vtrue[i] = V(i, j);

        v_gf.SetFromTrueDofs(vtrue);

        for (int i = 0; i < rw.size(); ++i)
        {
            const int e = qp[i] / nqe; // Element index

            if (e != eprev) // Update element transformation
            {
                doftrans = fesR->GetElementVDofs(e, vdofs);

                for (int k = 0; k < dof; ++k)
                { // Q: Why this code?
                    const int dofk = (vdofs[k] >= 0) ? vdofs[k] : -1 - vdofs[k];
                    const double vk = (vdofs[k] >= 0) ? v_gf[dofk] : -v_gf[dofk];
                    Vs((elemCount * dof) + k, j) = vk;
                }

                eprev = e;
                elemCount++;
            }
        }
    }

    eprev = -1;
    elemCount = 0;
    int spaceDim = 0;

    for (int i = 0; i < rw.size(); ++i)
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

            if (doftrans)
            {
                MFEM_ABORT("TODO");
            }

            dof = fe->GetDof();
            spaceDim = eltrans->GetSpaceDim();
            trial_vshape.SetSize(dof, spaceDim);
            Vx.SetSize(spaceDim);

            // Initialize nonlinear operator matrices (there is probably a better way)
            int dim = el.GetDim();
            DenseMatrix DSh(dof, dim);
            DenseMatrix DS(dof, dim);
            DenseMatrix Jrt(dim);
            DenseMatrix Jpt(dim);
            DenseMatrix P(dim);
            DenseMatrix PMatI;
            PMatI.UseExternalData(elfun.GetData(), dof, dim); // Q: How to replace?
            // PMatO.SetSize(dof, dim); // Q: should this be here?

            eprev = e;
            elemCount++;
        }

        // Integrate at the current point

        eltrans->SetIntPoint(&ip);
        fe->CalcVShape(*eltrans, trial_vshape);

        double w = eltrans->Weight() * rw[i]; // using rw[i] instead of ip.weight

        // Compute action of nonlinear operator
        CalcInverse(eltrans.Jacobian(), Jrt);
        fe.CalcDShape(ip, DSh);
        Mult(DSh, Jrt, DS);
        MultAtB(PMatI, DS, Jpt);
        model->EvalP(Jpt, P);

        // AddMultABt(DS, P, PMatO); //Q: Should this be here?

        for (int jx = 0; jx < rdim; ++jx)
        {
            // Lift Vx = V_{jx} at ip, where x = e_{jx}.
            Vx = 0.0;
            for (int k = 0; k < dof; ++k)
            {
                const double Vx_k = Vs(((elemCount - 1) * dof) + k, jx);

                for (int j = 0; j < spaceDim; ++j)
                    Vx[j] += Vx_k * trial_vshape(k, j);
            }

            for (int j = 0; j < rdim; ++j)
            {
                double rj = 0.0;
                for (int k = 0; k < spaceDim; ++k)
                {
                    double Vjk = 0.0;
                    for (int l = 0; l < dof; ++l)
                    {
                        Vjk += Vs(((elemCount - 1) * dof) + l, j) * trial_vshape(l, k);
                    }

                    rj += Vx[k] * Vjk;
                }

                res[j + (jx * rdim) + (i * rdim * rdim)] = w * rj;
            }
        }
    }
}

*/

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

    nnls.normalize_constraints(Gt, rhs_lb, rhs_ub);
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
                        ParFiniteElementSpace *fespace_H,
                        const int nsets, const CAROM::Matrix *BH,
                        const CAROM::Matrix *BH_snapshots,
                        const bool precondition, const double nnls_tol,
                        const int maxNNLSnnz, HyperelasticOperator oper, NeoHookeanModel *model,
                        CAROM::Vector &sol)
{
    const int nqe = ir0->GetNPoints();
    const int ne = fespace_H->GetNE();
    const int NB = BH->numColumns();
    const int NQ = ne * nqe;
    const int nsnap = BH_snapshots->numColumns();

    MFEM_VERIFY(nsnap == BH_snapshots->numColumns() ||
                    nsnap + nsets == BH_snapshots->numColumns(), // Q: nsets?
                "");
    MFEM_VERIFY(BH->numRows() == BH_snapshots->numRows(), "");
    MFEM_VERIFY(BH->numRows() == fespace_H->GetTrueVSize(), "");

    const bool skipFirstW = (nsnap + nsets == BH_snapshots->numColumns());

    // Compute G of size (NB * nsnap) x NQ, but only store its transpose Gt.
    CAROM::Matrix Gt(NQ, NB * nsnap, true);

    // For 0 <= j < NB, 0 <= i < nsnap, 0 <= e < ne, 0 <= m < nqe,
    // G(j + (i*NB), (e*nqe) + m)
    // is the coefficient of v_j^T M(p_i) V v_i at point m of element e,
    // with respect to the integration rule weight at that point,
    // where the "exact" quadrature solution is ir0->GetWeights().

    Vector h_i(BH_snapshots->numRows());
    Vector v_j(BH->numRows());

    Vector r(nqe);

    int skip = 0;
    const int nsnapPerSet = nsnap / nsets;
    if (skipFirstW)
    {
        MFEM_VERIFY(nsets * nsnapPerSet == nsnap, "");
        skip = 1;
    }

    // Get prolongation matrix
    const Operator *P = fespace_H->GetProlongationMatrix();
    if (!P)
    {
        MFEM_ABORT("P is null, generalize to serial case")
    }

    // For every snapshot
    for (int i = 0; i < nsnap; ++i)
    {
        // Set the sampled dofs from the snapshot matrix
        for (int j = 0; j < BH_snapshots->numRows(); ++j)
            h_i[j] = (*BH_snapshots)(j, i);

        // Get prolongated dofs
        Vector ph_i;
        Vector elfun;

        ph_i.SetSize(P->Height());
        P->Mult(h_i, ph_i);

        if (skipFirstW && i > 0 && i % nsnapPerSet == 0)
            skip++;

        // For each basis vector
        for (int j = 0; j < NB; ++j)
        {

            // Get basis vector
            for (int k = 0; k < BH->numRows(); ++k)
                v_j[k] = (*BH)(k, j);

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
                DofTransformation *doftrans = fespace_H->GetElementVDofs(e, vdofs);
                const FiniteElement &fe = *fespace_H->GetFE(e);
                ElementTransformation *eltrans = fespace_H->GetElementTransformation(e);
                ph_i.GetSubVector(vdofs, elfun);
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
