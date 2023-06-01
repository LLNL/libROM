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
   
}

// Compute coefficients of the reduced integrator with respect to inputs Q and x
// in HyperelasticNLFIntegrator_ComputeReducedEQP.
void GetEQPCoefficients_HyperelasticNLFIntegrator(ParFiniteElementSpace *fesR,
        std::vector<double> const& rw, std::vector<int> const& qp,
        const IntegrationRule *ir,
        CAROM::Matrix const& V, Vector & res)
{
   
}

void HyperelasticNLFIntegrator_ComputeReducedEQP(ParFiniteElementSpace *fesR,
        std::vector<double> const& rw, std::vector<int> const& qp,
        const IntegrationRule *ir, Coefficient *Q,
        CAROM::Matrix const& V, CAROM::Vector const& x, const int rank, Vector & res)
{
   
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
            PreconditionNNLS(fespace_R, new HyperelasticNLFIntegrator(a_coeff), BR, i, Gt);
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
