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
// and element matrix computation from VectorFEMassIntegrator::AssembleElementMatrix.
void AssembleElementMatrix_VectorFEMassIntegrator(Coefficient *Q,
        const FiniteElement &el,
        ElementTransformation &Trans,
        DenseMatrix &elmat)
{
    int dof = el.GetDof();
    int spaceDim = Trans.GetSpaceDim();

    double w;

    DenseMatrix trial_vshape(dof, spaceDim);

    elmat.SetSize(dof);
    elmat = 0.0;

    const IntegrationRule *ir = NULL;
    if (ir == NULL)
    {
        // int order = 2 * el.GetOrder();
        int order = Trans.OrderW() + 2 * el.GetOrder();
        ir = &IntRules.Get(el.GetGeomType(), order);
    }

    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);

        Trans.SetIntPoint(&ip);

        el.CalcVShape(Trans, trial_vshape);

        w = ip.weight * Trans.Weight();

        if (Q)
        {
            w *= Q -> Eval (Trans, ip);
        }
        AddMult_a_AAt (w, trial_vshape, elmat);
    }
}

// Compute coefficients of the reduced integrator with respect to inputs Q and x
// in VectorFEMassIntegrator_ComputeReducedEQP.
void GetEQPCoefficients_VectorFEMassIntegrator(ParFiniteElementSpace *fesR,
        std::vector<double> const& rw, std::vector<int> const& qp,
        const IntegrationRule *ir,
        CAROM::Matrix const& V, Vector & res)
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
        fe->CalcVShape(*eltrans, trial_vshape);

        double w = eltrans->Weight() * rw[i]; // using rw[i] instead of ip.weight

        // Note that the coefficient Q is omitted: w *= Q -> Eval(*eltrans, ip);

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

void VectorFEMassIntegrator_ComputeReducedEQP(ParFiniteElementSpace *fesR,
        std::vector<double> const& rw, std::vector<int> const& qp,
        const IntegrationRule *ir, Coefficient *Q,
        CAROM::Matrix const& V, CAROM::Vector const& x, const int rank, Vector & res)
{
    const int rdim = V.numColumns();
    MFEM_VERIFY(rw.size() == qp.size(), "");
    MFEM_VERIFY(x.dim() == rdim, "");
    MFEM_VERIFY(V.numRows() == fesR->GetTrueVSize(), "");

    MFEM_VERIFY(rank == 0,
                "TODO: generalize to parallel. This uses full dofs in V, which has true dofs");

    const int nqe = ir->GetWeights().Size();

    ElementTransformation *eltrans;
    DofTransformation * doftrans;
    const FiniteElement *fe = NULL;
    Array<int> vdofs;

    DenseMatrix trial_vshape;

    Vector Vx;
    res.SetSize(rdim);
    res = 0.0;

    int eprev = -1;
    int dof = 0;
    int spaceDim = 0;

    // Note that the alternative version VectorFEMassIntegrator_ComputeReducedEQP_Fast
    // of this function is optimized by storing some intermediate computations.

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
        }

        // Integrate at the current point

        eltrans->SetIntPoint(&ip);
        fe->CalcVShape(*eltrans, trial_vshape);

        double w = eltrans->Weight() * rw[i]; // using rw[i] instead of ip.weight

        if (Q)
        {
            w *= Q -> Eval(*eltrans, ip);
        }

        // Lift Vx at ip.
        Vx = 0.0;
        for (int k=0; k<dof; ++k)
        {
            const int dofk = (vdofs[k] >= 0) ? vdofs[k] : -1 - vdofs[k];

            double Vx_k = 0.0;
            for (int j=0; j<rdim; ++j)
            {
                Vx_k += V(dofk, j) * x(j);
            }

            if (vdofs[k] < 0) Vx_k = -Vx_k;

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
                    const int dofl = (vdofs[l] >= 0) ? vdofs[l] : -1 - vdofs[l];
                    const double s = (vdofs[l] >= 0) ? 1.0 : -1.0;
                    Vjk += s * V(dofl, j) * trial_vshape(l, k);
                }

                rj += Vx[k] * Vjk;
            }

            res[j] += w * rj;
        }
    }
}

void VectorFEMassIntegrator_ComputeReducedEQP_Fast(ParFiniteElementSpace *fesR,
        std::vector<int> const& qp, const IntegrationRule *ir,
        Coefficient *Q, CAROM::Vector const& x,
        Vector const& coef, Vector & res)
{
    const int rdim = x.dim();
    const int nqe = ir->GetWeights().Size();
    ElementTransformation *eltrans;

    res.SetSize(rdim);
    res = 0.0;

    int eprev = -1;

    for (int i=0; i<qp.size(); ++i)
    {
        const int e = qp[i] / nqe;  // Element index
        // Local (element) index of the quadrature point
        const int qpi = qp[i] - (e*nqe);
        const IntegrationPoint &ip = ir->IntPoint(qpi);

        if (e != eprev)  // Update element transformation
        {
            eltrans = fesR->GetElementTransformation(e);
            eprev = e;
        }

        eltrans->SetIntPoint(&ip);
        const double q_ip = Q -> Eval(*eltrans, ip);

        for (int j=0; j<rdim; ++j)
        {
            for (int k=0; k<rdim; ++k)
            {
                res[j] += coef[j + (k * rdim) + (i * rdim * rdim)] * q_ip * x(k);
            }
        }
    }
}

// Compute coefficients of the reduced integrator with respect to input Q
// in LinearMassIntegrator_ComputeReducedEQP.
void GetEQPCoefficients_LinearMassIntegrator(ParFiniteElementSpace *fes,
        std::vector<double> const& rw, std::vector<int> const& qp,
        const IntegrationRule *ir, CAROM::Matrix const& V,
        Vector & res)
{
    // Assumption: fes is an L2 space, where true DOFs and full DOFs are the same.
    // This allows for using V(dof,j), where dof is from fes->GetElementVDofs.
    MFEM_VERIFY(fes->GetTrueVSize() == fes->GetTrueVSize(), "");

    const int rdim = V.numColumns();
    MFEM_VERIFY(rw.size() == qp.size(), "");
    MFEM_VERIFY(V.numRows() == fes->GetTrueVSize(), "");

    const int nqe = ir->GetWeights().Size();

    ElementTransformation *eltrans;
    DofTransformation * doftrans;
    const FiniteElement *fe = NULL;
    Array<int> vdofs;

    Vector shape;

    int eprev = -1;
    int dof = 0;

    for (int i=0; i<rw.size(); ++i)
    {
        const int e = qp[i] / nqe;  // Element index
        // Local (element) index of the quadrature point
        const int qpi = qp[i] - (e*nqe);
        const IntegrationPoint &ip = ir->IntPoint(qpi);

        if (e != eprev)  // Update element transformation
        {
            doftrans = fes->GetElementVDofs(e, vdofs);
            fe = fes->GetFE(e);
            eltrans = fes->GetElementTransformation(e);

            if (doftrans)
            {
                MFEM_ABORT("TODO");
            }

            dof = fe->GetDof();

            shape.SetSize(dof);

            eprev = e;

            if (i == 0)
            {
                res.SetSize(rdim * rw.size() * dof);
                res = 0.0;
            }
        }

        // Integrate at the current point

        eltrans->SetIntPoint(&ip);

        fe->CalcPhysShape(*eltrans, shape);

        double w = eltrans->Weight() * rw[i]; // using rw[i] instead of ip.weight

        for (int j=0; j<rdim; ++j)
        {
            double Vj = 0.0;
            for (int l=0; l<dof; ++l)
            {
                const int dofl = (vdofs[l] >= 0) ? vdofs[l] : -1 - vdofs[l];
                const double s = (vdofs[l] >= 0) ? 1.0 : -1.0;
                res[l + (j * dof) + (i * rdim * dof)] = w * s * V(dofl,
                                                        j) * shape(l) * shape(l);
            }
        }
    }
}

// Based on MassIntegrator::AssembleElementMatrix
void LinearMassIntegrator_ComputeReducedEQP(ParFiniteElementSpace *fes,
        std::vector<double> const& rw, std::vector<int> const& qp,
        const IntegrationRule *ir, Coefficient *Q,
        CAROM::Matrix const& V, CAROM::Vector & res)
{
    const int rdim = V.numColumns();
    MFEM_VERIFY(rw.size() == qp.size(), "");
    MFEM_VERIFY(res.dim() == rdim, "");
    MFEM_VERIFY(V.numRows() == fes->GetTrueVSize(), "");

    ParGridFunction f_gf(fes);
    f_gf.ProjectCoefficient(*Q);
    Vector sv(fes->GetTrueVSize());
    f_gf.GetTrueDofs(sv);

    const int nqe = ir->GetWeights().Size();

    ElementTransformation *eltrans;
    DofTransformation * doftrans;
    const FiniteElement *fe = NULL;
    Array<int> vdofs;

    Vector shape;

    res = 0.0;

    int eprev = -1;
    int dof = 0;

    // Note that the alternative version LinearMassIntegrator_ComputeReducedEQP_Fast
    // of this function is optimized by storing some intermediate computations.

    for (int i=0; i<rw.size(); ++i)
    {
        const int e = qp[i] / nqe;  // Element index
        // Local (element) index of the quadrature point
        const int qpi = qp[i] - (e*nqe);
        const IntegrationPoint &ip = ir->IntPoint(qpi);

        if (e != eprev)  // Update element transformation
        {
            doftrans = fes->GetElementVDofs(e, vdofs);
            fe = fes->GetFE(e);
            eltrans = fes->GetElementTransformation(e);

            if (doftrans)
            {
                MFEM_ABORT("TODO");
            }

            dof = fe->GetDof();

            shape.SetSize(dof);

            eprev = e;
        }

        // Integrate at the current point

        eltrans->SetIntPoint(&ip);

        fe->CalcPhysShape(*eltrans, shape);

        double w = eltrans->Weight() * rw[i]; // using rw[i] instead of ip.weight

        //w *= Q -> Eval(*eltrans, ip);

        for (int j=0; j<rdim; ++j)
        {
            double Vj = 0.0;
            for (int l=0; l<dof; ++l)
            {
                const int dofl = (vdofs[l] >= 0) ? vdofs[l] : -1 - vdofs[l];
                const double s = (vdofs[l] >= 0) ? 1.0 : -1.0;
                //Vj += s * V(dofl, j) * shape(l);
                Vj += s * V(dofl, j) * shape(l) * sv[dofl] * shape(l);
            }

            res(j) += w * Vj;
        }
    }
}

void LinearMassIntegrator_ComputeReducedEQP_Fast(ParFiniteElementSpace *fes,
        ParGridFunction & f_gf,
        std::vector<double> const& rw, std::vector<int> const& qp,
        const IntegrationRule *ir, Coefficient *Q,
        CAROM::Matrix const& V, Vector const& coef,
        CAROM::Vector & res)
{
    const int rdim = V.numColumns();
    MFEM_VERIFY(rw.size() == qp.size(), "");
    MFEM_VERIFY(V.numRows() == fes->GetTrueVSize(), "");

    const int nqe = ir->GetWeights().Size();

    DofTransformation * doftrans;
    const FiniteElement *fe = NULL;
    Array<int> vdofs;

    res = 0.0;

    int eprev = -1;
    int dof = 0;

    for (int i=0; i<rw.size(); ++i)
    {
        const int e = qp[i] / nqe;  // Element index
        // Local (element) index of the quadrature point
        const int qpi = qp[i] - (e*nqe);
        const IntegrationPoint &ip = ir->IntPoint(qpi);

        if (e != eprev)  // Update element transformation
        {
            doftrans = fes->GetElementVDofs(e, vdofs);
            fe = fes->GetFE(e);

            if (doftrans)
            {
                MFEM_ABORT("TODO");
            }

            dof = fe->GetDof();

            eprev = e;

            f_gf.ProjectCoefficient(*Q, vdofs);

            if (i == 0)
            {
                MFEM_VERIFY(coef.Size() == rdim * rw.size() * dof, "");
            }
        }

        // Integrate at the current point

        for (int j=0; j<rdim; ++j)
        {
            double Vj = 0.0;
            for (int l=0; l<dof; ++l)
            {
                const int dofl = (vdofs[l] >= 0) ? vdofs[l] : -1 - vdofs[l];
                Vj += f_gf[dofl] * coef[l + (j * dof) + (i * rdim * dof)];
            }

            res(j) += Vj;
        }
    }
}

// Compute a(p) u . v at all quadrature points on the given element. Coefficient Q is a(p).
void ComputeElementRowOfG(const IntegrationRule *ir, Array<int> const& vdofs,
                          Coefficient *Q, Vector const& u, Vector const& v,
                          FiniteElement const& fe, ElementTransformation & Trans, Vector & r)
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

        for (int j=0; j<dof; ++j)
        {
            const int dofj = (vdofs[j] >= 0) ? vdofs[j] : -1 - vdofs[j];
            const double s = (vdofs[j] >= 0) ? 1.0 : -1.0;
            for (int k=0; k<spaceDim; ++k)
            {
                u_i[k] += s * u[dofj] * trial_vshape(j, k);
                v_i[k] += s * v[dofj] * trial_vshape(j, k);
            }
        }

        if (Q)
        {
            w *= Q -> Eval (Trans, ip);
        }

        r[i] = 0.0;
        for (int k=0; k<spaceDim; ++k)
        {
            r[i] += u_i[k] * v_i[k];
        }

        r[i] *= w;
    }
}

// Compute p*s at all quadrature points on the given element.
void ComputeElementRowOfG_Source(const IntegrationRule *ir,
                                 Array<int> const& vdofs,
                                 Vector const& s, Vector const& p,
                                 FiniteElement const& fe, ElementTransformation & Trans, Vector & r)
{
    MFEM_VERIFY(r.Size() == ir->GetNPoints(), "");
    int dof = fe.GetDof();

    Vector shape(dof);

    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);

        Trans.SetIntPoint(&ip);

        fe.CalcPhysShape(Trans, shape);

        double w = Trans.Weight();

        double u_i = 0.0;
        double v_i = 0.0;

        for (int j=0; j<dof; ++j)
        {
            const int dofj = (vdofs[j] >= 0) ? vdofs[j] : -1 - vdofs[j];
            const double sd = (vdofs[j] >= 0) ? 1.0 : -1.0;
            u_i += sd * p[dofj] * shape(j);
            v_i += sd * s[dofj] * shape(j);
        }

        r[i] = w * u_i * v_i;
    }
}

// Precondition Gt by (V^T M V)^{-1} (of size NB x NB).
void PreconditionNNLS(ParFiniteElementSpace *fespace,
                      BilinearFormIntegrator *massInteg,
                      const CAROM::Matrix* B,
                      const int snapshot,
                      CAROM::Matrix & Gt)
{
    const int NB = B->numColumns();
    const int NQ = Gt.numRows();
    const int nsnap = Gt.numColumns() / NB;

    ParBilinearForm M(fespace);

    M.AddDomainIntegrator(massInteg);
    M.Assemble();
    M.Finalize();
    HypreParMatrix *Mmat = M.ParallelAssemble();

    CAROM::Matrix Mhat(NB, NB, false);
    ComputeCtAB(*Mmat, *B, *B, Mhat);
    Mhat.inverse();

    CAROM::Vector PG(NB, false);

    for (int i = (snapshot >= 0 ? snapshot : 0);
            i < (snapshot >= 0 ? snapshot+1 : nsnap); ++i)
    {
        for (int m=0; m<NQ; ++m)
        {
            for (int j=0; j<NB; ++j)
            {
                PG(j) = 0.0;
                for (int k=0; k<NB; ++k)
                {
                    PG(j) += Mhat(j,k) * Gt(m, k + (i*NB));
                }
            }

            for (int j=0; j<NB; ++j)
                Gt(m, j + (i*NB)) = PG(j);
        }
    }
}

void SolveNNLS(const int rank, const double nnls_tol, const int maxNNLSnnz,
               CAROM::Vector const& w, CAROM::Matrix & Gt,
               CAROM::Vector & sol)
{
    CAROM::NNLSSolver nnls(nnls_tol, 0, maxNNLSnnz, 2);

    CAROM::Vector rhs_ub(Gt.numColumns(), false);
    // G.mult(w, rhs_ub);  // rhs = Gw
    // rhs = Gw. Note that by using Gt and multTranspose, we do parallel communication.
    Gt.transposeMult(w, rhs_ub);

    CAROM::Vector rhs_lb(rhs_ub);
    CAROM::Vector rhs_Gw(rhs_ub);

    const double delta = 1.0e-11;
    for (int i=0; i<rhs_ub.dim(); ++i)
    {
        rhs_lb(i) -= delta;
        rhs_ub(i) += delta;
    }

    nnls.normalize_constraints(Gt, rhs_lb, rhs_ub);
    nnls.solve_parallel_with_scalapack(Gt, rhs_lb, rhs_ub, sol);

    int nnz = 0;
    for (int i=0; i<sol.dim(); ++i)
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
    cout << rank << ": relative residual norm for NNLS solution of Gs = Gw: " <<
         relNorm << endl;
}

// Compute EQP solution from constraints on snapshots.
void SetupEQP_snapshots(const IntegrationRule *ir0, const int rank,
                        ParFiniteElementSpace *fespace_R,
                        ParFiniteElementSpace *fespace_W,
                        const int nsets, const CAROM::Matrix* BR,
                        const CAROM::Matrix* BR_snapshots,
                        const CAROM::Matrix* BW_snapshots,
                        const bool precondition, const double nnls_tol,
                        const int maxNNLSnnz,
                        CAROM::Vector & sol)
{
    const int nqe = ir0->GetNPoints();
    const int ne = fespace_R->GetNE();
    const int NB = BR->numColumns();
    const int NQ = ne * nqe;
    const int nsnap = BR_snapshots->numColumns();

    MFEM_VERIFY(nsnap == BW_snapshots->numColumns() ||
                nsnap + nsets == BW_snapshots->numColumns(), "");
    MFEM_VERIFY(BR->numRows() == BR_snapshots->numRows(), "");
    MFEM_VERIFY(BR->numRows() == fespace_R->GetTrueVSize(), "");
    MFEM_VERIFY(BW_snapshots->numRows() == fespace_W->GetTrueVSize(), "");

    const bool skipFirstW = (nsnap + nsets == BW_snapshots->numColumns());

    // Compute G of size (NB * nsnap) x NQ, but only store its transpose Gt.
    CAROM::Matrix Gt(NQ, NB * nsnap, true);

    // For 0 <= j < NB, 0 <= i < nsnap, 0 <= e < ne, 0 <= m < nqe,
    // G(j + (i*NB), (e*nqe) + m)
    // is the coefficient of v_j^T M(p_i) V v_i at point m of element e,
    // with respect to the integration rule weight at that point,
    // where the "exact" quadrature solution is ir0->GetWeights().

    Vector p_i(BW_snapshots->numRows());
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

    for (int i=0; i<nsnap; ++i)
    {
        for (int j = 0; j < BW_snapshots->numRows(); ++j)
            p_i[j] = (*BW_snapshots)(j, i + skip);

        for (int j = 0; j < BR_snapshots->numRows(); ++j)
            v_i[j] = (*BR_snapshots)(j, i);

        if (skipFirstW && i > 0 && i % nsnapPerSet == 0)
            skip++;

        // Set grid function for a(p)
        ParGridFunction p_gf(fespace_W);

        p_gf.SetFromTrueDofs(p_i);

        GridFunctionCoefficient p_coeff(&p_gf);
        TransformedCoefficient a_coeff(&p_coeff, NonlinearCoefficient);

        ParGridFunction vi_gf(fespace_R);
        vi_gf.SetFromTrueDofs(v_i);

        for (int j=0; j<NB; ++j)
        {
            for (int k = 0; k < BR->numRows(); ++k)
                v_j[k] = (*BR)(k, j);

            ParGridFunction vj_gf(fespace_R);
            vj_gf.SetFromTrueDofs(v_j);

            // TODO: is it better to make the element loop the outer loop?
            for (int e=0; e<ne; ++e)
            {
                Array<int> vdofs;
                DofTransformation *doftrans = fespace_R->GetElementVDofs(e, vdofs);
                const FiniteElement &fe = *fespace_R->GetFE(e);
                ElementTransformation *eltrans = fespace_R->GetElementTransformation(e);

                ComputeElementRowOfG(ir0, vdofs, &a_coeff, vi_gf, vj_gf, fe, *eltrans, r);

                for (int m=0; m<nqe; ++m)
                    Gt((e*nqe) + m, j + (i*NB)) = r[m];
            }
        }

        if (precondition)
        {
            // Preconditioning is done by (V^T M(p_i) V)^{-1} (of size NB x NB).
            PreconditionNNLS(fespace_R, new VectorFEMassIntegrator(a_coeff), BR, i, Gt);
        }
    } // Loop (i) over snapshots

    Array<double> const& w_el = ir0->GetWeights();
    MFEM_VERIFY(w_el.Size() == nqe, "");

    CAROM::Vector w(ne * nqe, true);

    for (int i=0; i<ne; ++i)
    {
        for (int j=0; j<nqe; ++j)
            w((i*nqe) + j) = w_el[j];
    }

    SolveNNLS(rank, nnls_tol, maxNNLSnnz, w, Gt, sol);
}

// Compute EQP solution from constraints on source snapshots.
void SetupEQP_S_snapshots(const IntegrationRule *ir0, const int rank,
                          ParFiniteElementSpace *fespace_W,
                          const CAROM::Matrix* BW,
                          const CAROM::Matrix* BS_snapshots,
                          const bool precondition, const double nnls_tol,
                          const int maxNNLSnnz,
                          CAROM::Vector & sol)
{
    const int nqe = ir0->GetNPoints();
    const int ne = fespace_W->GetNE();
    const int NB = BW->numColumns();
    const int NQ = ne * nqe;
    const int nsnap = BS_snapshots->numColumns();
    const int nrows = BW->numRows();

    MFEM_VERIFY(BS_snapshots->numRows() == BW->numRows(), "");

    // Compute G of size (NB * nsnap) x NQ, but only store its transpose Gt.
    CAROM::Matrix Gt(NQ, NB * nsnap, true);

    Vector s_i(nrows);
    Vector p_j(nrows);

    Vector r(nqe);

    for (int i=0; i<nsnap; ++i)
    {
        for (int j = 0; j < BS_snapshots->numRows(); ++j)
            s_i[j] = (*BS_snapshots)(j, i);

        for (int j=0; j<NB; ++j)
        {
            for (int k = 0; k < BW->numRows(); ++k)
                p_j[k] = (*BW)(k, j);

            // TODO: is it better to make the element loop the outer loop?
            for (int e=0; e<ne; ++e)
            {
                Array<int> vdofs;
                DofTransformation *doftrans = fespace_W->GetElementVDofs(e, vdofs);
                const FiniteElement &fe = *fespace_W->GetFE(e);
                ElementTransformation *eltrans = fespace_W->GetElementTransformation(e);

                ComputeElementRowOfG_Source(ir0, vdofs, s_i, p_j, fe, *eltrans, r);

                for (int m=0; m<nqe; ++m)
                    Gt((e*nqe) + m, j + (i*NB)) = r[m];
            }
        }
    }

    if (precondition)
    {
        // Preconditioning is done by (V^T M V)^{-1} (of size NB x NB).
        PreconditionNNLS(fespace_W, new MassIntegrator(), BW, -1, Gt);
    }

    Array<double> const& w_el = ir0->GetWeights();
    MFEM_VERIFY(w_el.Size() == nqe, "");

    CAROM::Vector w(ne * nqe, true);

    for (int i=0; i<ne; ++i)
    {
        for (int j=0; j<nqe; ++j)
            w((i*nqe) + j) = w_el[j];
    }

    SolveNNLS(rank, nnls_tol, maxNNLSnnz, w, Gt, sol);
}

void WriteMeshEQP(ParMesh *pmesh, const int myid, const int nqe,
                  CAROM::Vector const& eqpSol)
{
    // Find the elements with quadrature points in eqpSol.
    std::set<int> elements;
    Vector elemCount(pmesh->GetNE());
    elemCount = 0.0;

    for (int i=0; i<eqpSol.dim(); ++i)
    {
        if (eqpSol(i) > 1.0e-12)
        {
            const int e = i / nqe;  // Element index
            elements.insert(e);
            elemCount[e] += 1.0;
        }
    }

    // Empty sets, since EQP on samples inside elements.
    std::set<int> faces, edges, vertices;
    CAROM::SampleVisualization(pmesh, elements, elements, faces, edges,
                               vertices, "EQPvis", &elemCount);
}
