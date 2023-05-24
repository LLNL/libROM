/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
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
     * @param[in] alt_output_basis  Whether to use the alternative basis for
     *                              output, i.e. phi = U^(+)*V*Omega^(-1)*X.
     * @param[in] state_offset      The state offset.
     * @param[in] derivative_offset The derivative offset.
     */
    NonuniformDMD(int dim,
                  bool alt_output_basis = false,
                  Vector* state_offset = NULL,
                  Vector* derivative_offset = NULL);

    /**
     * @brief Constructor.
     *
     * @param[in] base_file_name The base part of the filename of the
     *                           database to load when restarting from a save.
     */
    NonuniformDMD(std::string base_file_name);

    /**
     * @brief Destroy the NonuniformDMD object
     */
    ~NonuniformDMD();

    /**
     * @brief Set the offset of a certain order.
     */
    void setOffset(Vector* offset_vector, int order) override;

    /**
     * @brief Load the object state to a file.
     *
     * @param[in] base_file_name The base part of the filename to load the
     *                           database to.
     */
    void load(std::string base_file_name) override;

    /**
     * @brief Save the object state to a file.
     *
     * @param[in] base_file_name The base part of the filename to save the
     *                           database to.
     */
    void save(std::string base_file_name) override;

protected:
    /**
     * @brief Obtain DMD model interpolant at desired parameter point by
     *        interpolation of DMD models from training parameter points.
     *
     * @param[in] parametric_dmd    The interpolant DMD model at the desired point.
     * @param[in] parameter_points  The training parameter points.
     * @param[in] dmds              The DMD objects associated with
     *                              each training parameter point.
     * @param[in] desired_point     The desired point at which to create a parametric DMD.
     * @param[in] rbf               The RBF type ("G" == gaussian,
     *                              "IQ" == inverse quadratic,
     *                              "IMQ" == inverse multiquadric)
     * @param[in] interp_method     The interpolation method type
     *                              ("LS" == linear solve,
     *                              "IDW" == inverse distance weighting,
     *                              "LP" == lagrangian polynomials)
     * @param[in] closest_rbf_val   The RBF parameter determines the width of influence.
     *                              Set the RBF value of the nearest two parameter points to a value between 0.0 to 1.0
     * @param[in] reorthogonalize_W Whether to reorthogonalize the interpolated W (basis) matrix.
     */
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
     * @param[in] dt d_dt
     * @param[in] t_offset d_t_offset
     * @param[in] state_offset d_state_offset
     * @param[in] derivative_offset d_derivative_offset
     */
    NonuniformDMD(std::vector<std::complex<double>> eigs, Matrix* phi_real,
                  Matrix* phi_imaginary, int k,
                  double dt, double t_offset,
                  Vector* state_offset, Vector* derivative_offset);

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
    std::pair<Matrix*, Matrix*> computeDMDSnapshotPair(const Matrix* snapshots)
    override;

    /**
     * @brief Compute phi.
     */
    void computePhi(DMDInternal dmd_internal_obj) override;

    /**
     * @brief Compute the appropriate exponential function when predicting the solution.
     */
    std::complex<double> computeEigExp(std::complex<double> eig, double t) override;

    /**
     * @brief Add the appropriate offset when predicting the solution.
     */
    void addOffset(Vector*& result, double t, int deg) override;

    /**
     * @brief Derivative offset in snapshot.
     */
    Vector* d_derivative_offset = NULL;

};

}

#endif
