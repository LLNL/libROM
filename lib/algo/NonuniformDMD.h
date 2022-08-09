/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the modified DMD algorithm on the given snapshot matrix
//              with non-uniform sampling time steps.
//              The implemented dynamic mode decomposition algorithm is modified from
//              Tu et. al's paper "On Dynamic Mode Decomposition: Theory and
//              Applications": https://arxiv.org/abs/1312.0041
//              Instead of approximating the discrete dynamics, this algorithm
//              approximates the continuous dynamics linearly.
//              This algorithm also works in the case that the first sample does
//              not start from t = 0.0 by incorporating a time offset.

#ifndef included_NonuniformDMD_h
#define included_NonuniformDMD_h

#include "DMD.h"
#include "ParametricDMD.h"

namespace CAROM {

/**
 * Class NonuniformDMD implements the modified DMD algorithm on the
 * given snapshot matrix with non-uniform sampling time steps.
 * Instead of linearly approximating the discrete dynamics
 * x(t+dt) = Ax(t) in the original DMD, this algorithm approximates
 * the continuous dynamics linearly by dx/dt = Ax.
 */
class NonuniformDMD : public DMD
{
public:

    /**
     * @brief Constructor.
     *
     * @param[in] dim               The full-order state dimension.
     * @param[in] in_offset         Whether to use state offset
     * @param[in] out_offset        Whether to use derivative offset.
     * @param[in] state_offset       The state offset.
     * @param[in] derivative_offset The derivative offset.
     */
    NonuniformDMD(int dim, bool in_offset = false, bool out_offset = false,  
                  Vector* state_offset = NULL, Vector* derivative_offset = NULL); 

    /**
     * @brief Constructor.
     *
     * @param[in] base_file_name The base part of the filename of the
     *                           database to load when restarting from a save.
     */
    NonuniformDMD(std::string base_file_name);

protected:
    friend void getParametricDMD<NonuniformDMD>(NonuniformDMD*& parametric_dmd,
            std::vector<Vector*>& parameter_points,
            std::vector<NonuniformDMD*>& dmds,
            Vector* desired_point,
            std::string rbf,
            std::string interp_method,
            double closest_rbf_val,
            bool reorthogonalize_W);

    /**
     * @brief Constructor.
     *
     * @param[in] eigs d_eigs
     * @param[in] phi_real d_phi_real
     * @param[in] phi_imaginary d_phi_imaginary
     * @param[in] k d_k
     * @param[in] in_offset d_in_offset
     * @param[in] out_offset d_out_offset
     * @param[in] state_offset d_state_offset
     * @param[in] derivative_offset d_derivative_offset
     * @param[in] dt d_dt
     * @param[in] t_offset d_t_offset
     */
    NonuniformDMD(std::vector<std::complex<double>> eigs, Matrix* phi_real,
                  Matrix* phi_imaginary, int k,
                  bool in_offset, bool out_offset,
                  Vector* state_offset, Vector* derivative_offset,
                  double dt, double t_offset);

private:

    /**
     * @brief Unimplemented default constructor.
     */
    NonuniformDMD();

    /**
     * @brief Unimplemented copy constructor.
     */
    NonuniformDMD(
        const NonuniformDMD& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    NonuniformDMD&
    operator = (
        const NonuniformDMD& rhs);

    /**
     * @brief Returns a pair of pointers to the state and derivative snapshot matrices
     */
    std::pair<Matrix*, Matrix*> computeDMDSnapshotPair(const Matrix* snapshots);

    /**
     * @brief Compute phi.
     */
    void computePhi(struct DMDInternal dmd_internal_obj);

    /**
     * @brief Compute the appropriate exponential function when predicting the solution.
     */
    std::complex<double> computeEigExp(std::complex<double> eig, double t);

};

}

#endif
