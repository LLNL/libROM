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

using namespace mfem;
using namespace std;

// Element matrix assembly, copying element loop from BilinearForm::Assemble(int skip_zeros)
// and element matrix computation from HyperelasticNLFIntegrator::AssembleElementMatrix.
void AssembleElementMatrix_HyperelasticNLFIntegrator(Coefficient *Q,
                                                     const FiniteElement &el,
                                                     ElementTransformation &Trans,
                                                     DenseMatrix &elmat)
{
    int dof = el.GetDof(), dim = el.GetDim();

    DSh.SetSize(dof, dim);
    DS.SetSize(dof, dim);
    Jrt.SetSize(dim);
    Jpt.SetSize(dim);
    P.SetSize(dim);
    PMatI.UseExternalData(elfun.GetData(), dof, dim);
    elvect.SetSize(dof * dim);
    PMatO.UseExternalData(elvect.GetData(), dof, dim);

    const IntegrationRule *ir = IntRule;
    if (!ir)
    {
        ir = &(IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + 3)); // <---
    }

    elvect = 0.0;
    model->SetTransformation(Ttr);
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        Ttr.SetIntPoint(&ip);
        CalcInverse(Ttr.Jacobian(), Jrt);

        el.CalcDShape(ip, DSh);
        Mult(DSh, Jrt, DS);
        MultAtB(PMatI, DS, Jpt);

        model->EvalP(Jpt, P);

        P *= ip.weight * Ttr.Weight();
        AddMultABt(DS, P, PMatO);
    }
}

// Compute coefficients of the reduced integrator with respect to inputs Q and x
// in HyperelasticNLFIntegrator_ComputeReducedEQP.
void GetEQPCoefficients_HyperelasticNLFIntegrator(ParFiniteElementSpace *fesR,
                                                  std::vector<double> const &rw, std::vector<int> const &qp,
                                                  const IntegrationRule *ir,
                                                  CAROM::Matrix const &V, Vector &res)
{
     const int rdim = V.numColumns();
    MFEM_VERIFY(V.numRows() == fesR->GetTrueVSize(), "");
    MFEM_VERIFY(rw.size() == qp.size(), "");

    const int ne = fesR->GetNE();
    const int nqe = ir->GetWeights().Size();

    ElementTransformation *eltrans;
    DofTransformation * doftrans;
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
    for (int i=0; i<rw.size(); ++i)
    {
        const int e = qp[i] / nqe;  // Element index

        if (e != eprev)  // Update element transformation
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

}

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
    DofTransformation * doftrans;
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
    for (int i=0; i<rw.size(); ++i)
    {
        const int e = qp[i] / nqe;  // Element index

        if (e != eprev)  // Update element transformation
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

    for (int j=0; j<rdim; ++j)
    {
        eprev = -1;
        elemCount = 0;

        for (int i=0; i<vtrue.Size(); ++i)
            vtrue[i] = V(i,j);

        v_gf.SetFromTrueDofs(vtrue);

        for (int i=0; i<rw.size(); ++i)
        {
            const int e = qp[i] / nqe;  // Element index

            if (e != eprev)  // Update element transformation
            {
                doftrans = fesR->GetElementVDofs(e, vdofs);

                for (int k=0; k<dof; ++k)
                {
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

    for (int i=0; i<rw.size(); ++i)
    {
        const int e = qp[i] / nqe;  // Element index
        // Local (element) index of the quadrature point
        const int qpi = qp[i] - (e*nqe);
        const IntegrationPoint &ip = ir->IntPoint(qpi);

        if (e != eprev)  // Update element transformation
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

            eprev = e;
            elemCount++;
        }

        // Integrate at the current point

        eltrans->SetIntPoint(&ip);
        // Ttr.SetIntPoint(&ip); From hyperelastic integrator
        fe->CalcVShape(*eltrans, trial_vshape);

        double w = eltrans->Weight() * rw[i]; // using rw[i] instead of ip.weight

        /* CalcInverse(Ttr.Jacobian(), Jrt);

        el.CalcDShape(ip, DSh);
        Mult(DSh, Jrt, DS);
        MultAtB(PMatI, DS, Jpt);

        model->EvalP(Jpt, P);

        P *= ip.weight * Ttr.Weight(); */

        for (int jx=0; jx<rdim; ++jx)
        {
            // Lift Vx = V_{jx} at ip, where x = e_{jx}.
            Vx = 0.0;
            for (int k=0; k<dof; ++k)
            {
                const double Vx_k = Vs(((elemCount-1) * dof) + k, jx);

                for (int j=0; j<spaceDim; ++j)
                    Vx[j] += Vx_k * trial_vshape(k, j);
            }

            for (int j=0; j<rdim; ++j)
            {
                double rj = 0.0;
                for (int k=0; k<spaceDim; ++k)
                {
                    double Vjk = 0.0;
                    for (int l=0; l<dof; ++l)
                    {
                        Vjk += Vs(((elemCount-1) * dof) + l, j) * trial_vshape(l, k);
                    }

                    rj += Vx[k] * Vjk;
                }

                res[j + (jx * rdim) + (i * rdim * rdim)] = w * rj;
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

// Compute a(p) u . v at all quadrature points on the given element. Coefficient Q is a(p).
void ComputeElementRowOfG(const IntegrationRule *ir, Array<int> const &vdofs,
                          Coefficient *Q, Vector const &u, Vector const &v,
                          FiniteElement const &fe, ElementTransformation &Trans, Vector &r)
{
    MFEM_VERIFY(r.Size() == ir->GetNPoints(), "");
    int dof = fe.GetDof();
    int spaceDim = Trans.GetSpaceDim();

    Vector u_i(spaceDim);
    Vector v_i(spaceDim);

    DenseMatrix trial_vshape(dof, spaceDim);

    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);

        Trans.SetIntPoint(&ip);

        fe.CalcVShape(Trans, trial_vshape);

        double w = Trans.Weight();

        u_i = 0.0;
        v_i = 0.0;

        for (int j = 0; j < dof; ++j)
        {
            const int dofj = (vdofs[j] >= 0) ? vdofs[j] : -1 - vdofs[j];
            const double s = (vdofs[j] >= 0) ? 1.0 : -1.0;
            for (int k = 0; k < spaceDim; ++k)
            {
                u_i[k] += s * u[dofj] * trial_vshape(j, k);
                v_i[k] += s * v[dofj] * trial_vshape(j, k);
            }
        }

        if (Q)
        {
            w *= Q->Eval(Trans, ip);
        }

        r[i] = 0.0;
        for (int k = 0; k < spaceDim; ++k)
        {
            r[i] += u_i[k] * v_i[k];
        }

        r[i] *= w;
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
                        ParFiniteElementSpace *fespace_R,
                        const int nsets, const CAROM::Matrix *BR,
                        const CAROM::Matrix *BR_snapshots,
                        const bool precondition, const double nnls_tol,
                        const int maxNNLSnnz,
                        CAROM::Vector &sol)
{
    const int nqe = ir0->GetNPoints();
    const int ne = fespace_R->GetNE();
    const int NB = BR->numColumns();
    const int NQ = ne * nqe;
    const int nsnap = BR_snapshots->numColumns();

    MFEM_VERIFY(nsnap == BR_snapshots->numColumns() ||
                    nsnap + nsets == BR_snapshots->numColumns(),
                "");
    MFEM_VERIFY(BR->numRows() == BR_snapshots->numRows(), "");
    MFEM_VERIFY(BR->numRows() == fespace_R->GetTrueVSize(), "");

    const bool skipFirstW = (nsnap + nsets == BR_snapshots->numColumns());

    // Compute G of size (NB * nsnap) x NQ, but only store its transpose Gt.
    CAROM::Matrix Gt(NQ, NB * nsnap, true);

    // For 0 <= j < NB, 0 <= i < nsnap, 0 <= e < ne, 0 <= m < nqe,
    // G(j + (i*NB), (e*nqe) + m)
    // is the coefficient of v_j^T M(p_i) V v_i at point m of element e,
    // with respect to the integration rule weight at that point,
    // where the "exact" quadrature solution is ir0->GetWeights().

    Vector v_i(BR_snapshots->numRows());
    Vector v_j(BR->numRows());

    Vector r(nqe);

    int skip = 0;
    const int nsnapPerSet = nsnap / nsets;
    if (skipFirstW)
    {
        MFEM_VERIFY(nsets * nsnapPerSet == nsnap, "");
        skip = 1;
    }

    for (int i = 0; i < nsnap; ++i)
    {
        for (int j = 0; j < BR_snapshots->numRows(); ++j)
            v_i[j] = (*BR_snapshots)(j, i);

        if (skipFirstW && i > 0 && i % nsnapPerSet == 0)
            skip++;

        // Set grid function for a(p)
        ParGridFunction p_gf(fespace_R);

        p_gf.SetFromTrueDofs(v_i);

        GridFunctionCoefficient p_coeff(&p_gf);
        TransformedCoefficient a_coeff(&p_coeff, NonlinearCoefficient);

        ParGridFunction vi_gf(fespace_R);
        vi_gf.SetFromTrueDofs(v_i);

        for (int j = 0; j < NB; ++j)
        {
            for (int k = 0; k < BR->numRows(); ++k)
                v_j[k] = (*BR)(k, j);

            ParGridFunction vj_gf(fespace_R);
            vj_gf.SetFromTrueDofs(v_j);

            // TODO: is it better to make the element loop the outer loop?
            for (int e = 0; e < ne; ++e)
            {
                Array<int> vdofs;
                DofTransformation *doftrans = fespace_R->GetElementVDofs(e, vdofs);
                const FiniteElement &fe = *fespace_R->GetFE(e);
                ElementTransformation *eltrans = fespace_R->GetElementTransformation(e);

                ComputeElementRowOfG(ir0, vdofs, &a_coeff, vi_gf, vj_gf, fe, *eltrans, r);

                for (int m = 0; m < nqe; ++m)
                    Gt((e * nqe) + m, j + (i * NB)) = r[m];
            }
        }

        if (precondition)
        {
            // TODO
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
