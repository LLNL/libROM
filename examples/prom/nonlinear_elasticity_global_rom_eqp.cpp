#include "mfem/Utilities.hpp"
// #include "fem/nonlininteg.hpp"
#include <cmath>
using namespace mfem;
using namespace std;

#include "nonlinear_elasticity_global_rom_eqp.hpp"

// Compute coefficients of the reduced integrator with respect to inputs Q and x
// in HyperelasticNLFIntegrator_ComputeReducedEQP.
void GetEQPCoefficients_HyperelasticNLFIntegrator(ParFiniteElementSpace *fesR,
                                                  std::vector<double> const &rw, std::vector<int> const &qp,
                                                  const IntegrationRule *ir, NeoHookeanModel *model,
                                                  CAROM::Matrix const &V_v, const int rank, Vector &coef, Vector &DS_coef)
{
    const int rvdim = V_v.numColumns();
    const int fomdim = V_v.numRows();
    MFEM_VERIFY(rw.size() == qp.size(), "");
    MFEM_VERIFY(fomdim == fesR->GetTrueVSize(), "");

    MFEM_VERIFY(rank == 0,
                "TODO: generalize to parallel. This uses full dofs in V, which has true dofs");

    const int nqe = ir->GetWeights().Size();

    DofTransformation *doftrans;
    ElementTransformation *eltrans;
    Array<int> vdofs;

    int eprev = -1;
    int dof = 0;
    int dim = 0;

    // Get prolongation matrix
    const Operator *P = fesR->GetProlongationMatrix();

    // Vector to be prolongated
    Vector vj(fomdim);

    // Prolongated vector
    Vector p_vj(P->Height());

    // Element vector
    Vector vj_e;
    // Get the element vector size
    doftrans = fesR->GetElementVDofs(0, vdofs);
    int elvect_size = vdofs.Size();

    // Coefficient vector
    coef.SetSize(elvect_size * rw.size() * rvdim);
    coef = 0.0;

    // Vector for storing DS
    const FiniteElement *fe = fesR->GetFE(0);
    dof = fe->GetDof();
    dim = fe->GetDim();
    DenseMatrix DSh(dof, dim);
    DenseMatrix DS(dof, dim);
    DenseMatrix Jrt(dim);
    DS_coef.SetSize(dof * dim * rw.size() * rvdim);
    DS_coef = 0.0;
    int index = 0;

    // For every basis vector
    for (int j = 0; j < rvdim; ++j)
    {
        // Get basis vector and prolongate
        for (int k = 0; k < V_v.numRows(); ++k)
            vj[k] = V_v(k, j);
        P->Mult(vj, p_vj);

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
                eltrans = fesR->GetElementTransformation(e);

                if (doftrans)
                {
                    MFEM_ABORT("TODO");
                }

                // Get element vectors
                p_vj.GetSubVector(vdofs, vj_e);
                eprev = e;
            }

            // Set integration point in the element transformation
            eltrans->SetIntPoint(&ip);
            model->SetTransformation(*eltrans);
            // Get the transformation weight
            double t = eltrans->Weight();
            // Calculate r[i] = ve_j^T * elvect
            for (int k = 0; k < elvect_size; k++)
            {
                coef[k + (i * elvect_size) + (j * rw.size() * elvect_size)] = vj_e[k] * rw[i] * t;
            }

            // Calculate DS and store
            CalcInverse(eltrans->Jacobian(), Jrt);
            fe->CalcDShape(ip, DSh);
            Mult(DSh, Jrt, DS);
            for (int ii = 0; ii < dof; ++ii)
            {
                for (int jj = 0; jj < dim; ++jj)
                {
                    index = jj + ii * dim;
                    DS_coef[index + (i * dof*dim) + (j * rw.size() * dof*dim)] = DS.Elem(ii,jj);
                }
            }
        }
    }
}
