/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the DMD algorithm on the given snapshot matrix. The
//              implemented dynamic mode decomposition algorithm is derived from
//              Tu et. al's paper "On Dynamic Mode Decomposition: Theory and
//              Applications": https://arxiv.org/abs/1312.0041
//              This algorithm also works in the case that the first sample does
//              not start from t = 0.0 by incorporating a time offset.

#ifndef included_DMD_h
#define included_DMD_h

#include "ParametricDMD.h"
#include <vector>
#include <complex>

namespace CAROM {

class Matrix;
class Vector;
class ComplexEigenPair;

/**
 * Struct DMDInternal is a struct containing the necessary matrices to compute phi.
 */
struct DMDInternal
{
    /**
     * @brief The snapshot vectors of input.
     */
    Matrix* snapshots_in;

    /**
     * @brief The snapshot vectors of output.
     */
    Matrix* snapshots_out;

    /**
     * @brief The left singular vector basis.
     */
    Matrix* basis;

    /**
     * @brief The right singular vector basis.
     */
    Matrix* basis_right;

    /**
     * @brief The inverse of singular values.
     */
    Matrix* S_inv;

    /**
     * @brief The resultant DMD eigenvalues and eigenvectors.
     */
    ComplexEigenPair* eigenpair;
};

/**
 * Class DMD implements the DMD algorithm on a given snapshot matrix.
 */
class DMD
{
public:

    /**
     * @brief Constructor. Basic DMD with uniform time step size.
     *
     * @param[in] dim              The full-order state dimension.
     * @param[in] dt               The dt between samples.
     * @param[in] alt_output_basis Whether to use the alternative basis for
     *                             output, i.e. phi = U^(+)*V*Omega^(-1)*X.
     * @param[in] state_offset     The state offset.
     */
    DMD(int dim, double dt, bool alt_output_basis = false,
        Vector* state_offset = NULL);

    /**
     * @brief Constructor. DMD from saved models.
     *
     * @param[in] base_file_name The base part of the filename of the
     *                           database to load when restarting from a save.
     */
    DMD(std::string base_file_name);

    /**
     * @brief Destroy the DMD object
     */
    virtual ~DMD();

    /**
     * @brief Set the offset of a certain order.
     */
    virtual void setOffset(Vector* offset_vector, int order);

    /**
     * @brief Sample the new state, u_in. Any samples in d_snapshots
     *        taken at the same or later time will be erased.
     *
     * @pre u_in != 0
     * @pre t >= 0.0
     *
     * @param[in] u_in The new state.
     * @param[in] t    The time of the newly sampled state.
     */
    virtual void takeSample(double* u_in, double t);

    /**
     * @param[in] energy_fraction The energy fraction to keep after doing SVD.
     * @param[in] W0              The initial basis to prepend to W.
     * @param[in] linearity_tol   The tolerance for determining whether a column
                                  of W is linearly independent with W0.
     */
    virtual void train(double energy_fraction, const Matrix* W0 = NULL,
                       double linearity_tol = 0.0);

    /**
     * @param[in] k               The number of modes to keep after doing SVD.
     * @param[in] W0              The initial basis to prepend to W.
     * @param[in] linearity_tol   The tolerance for determining whether a column
                                  of W is linearly independent with W0.
     */
    virtual void train(int k, const Matrix* W0 = NULL, double linearity_tol = 0.0);

    /**
     * @brief Project new initial condition using d_phi.
     *        Calculate pinv(phi) x init, or more precisely,
     *        (phi* x phi)^{-1} x phi* x init, where phi* is the conjugate transpose.
     *
     * @param[in] init     The initial condition.
     * @param[in] t_offset The initial time offset.
     */
    void projectInitialCondition(const Vector* init, double t_offset = -1.0);

    /**
     * @brief Predict state given a time. Uses the projected initial condition of the
     *        training dataset (the first column).
     *
     * @param[in] t   The time of the output state
     * @param[in] deg The derivative degree of the output state
     */
    Vector* predict(double t, int deg = 0);

    /**
     * @brief Get the time offset contained within d_t_offset.
     */
    double getTimeOffset() const;

    /**
     * @brief Returns the number of samples taken.
     *
     * @return The number of samples taken.
     */
    int getNumSamples() const
    {
        return d_snapshots.size();
    }

    int getDimension() const
    {
        return d_k;
    }

    /**
     * @brief Get the snapshot matrix contained within d_snapshots.
     */
    const Matrix* getSnapshotMatrix();

    /**
     * @brief Load the object state from a file.
     *
     * @param[in] base_file_name The base part of the filename to load the
     *                           database from.
     */
    virtual void load(std::string base_file_name);

    /**
     * @brief Load the object state from a file.
     *
     * @param[in] base_file_name The base part of the filename to load the
     *                           database from.
     */
    void load(const char* base_file_name);

    /**
     * @brief Save the object state to a file.
     *
     * @param[in] base_file_name The base part of the filename to save the
     *                           database to.
     */
    virtual void save(std::string base_file_name);

    /**
     * @brief Save the object state to a file.
     *
     * @param[in] base_file_name The base part of the filename to save the
     *                           database to.
     */
    void save(const char* base_file_name);

    /**
     * @brief Output the DMD record in CSV files.
     */
    void summary(std::string base_file_name);

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
    friend void getParametricDMD<DMD>(DMD*& parametric_dmd,
                                      std::vector<Vector*>& parameter_points,
                                      std::vector<DMD*>& dmds,
                                      Vector* desired_point,
                                      std::string rbf,
                                      std::string interp_method,
                                      double closest_rbf_val,
                                      bool reorthogonalize_W);

    /**
     * @brief Constructor. Variant of DMD with non-uniform time step size.
     *
     * @param[in] dim               The full-order state dimension.
     * @param[in] alt_output_basis  Whether to use the alternative basis for
     *                              output, i.e. phi = U^(+)*V*Omega^(-1)*X.
     * @param[in] state_offset      The state offset.
     */
    DMD(int dim, bool alt_output_basis = false, Vector* state_offset = NULL);

    /**
     * @brief Constructor. Specified from DMD components.
     *
     * @param[in] eigs d_eigs
     * @param[in] phi_real d_phi_real
     * @param[in] phi_imaginary d_phi_imaginary
     * @param[in] k d_k
     * @param[in] dt d_dt
     * @param[in] t_offset d_t_offset
     * @param[in] state_offset d_state_offset
     */
    DMD(std::vector<std::complex<double>> eigs, Matrix* phi_real,
        Matrix* phi_imaginary, int k,
        double dt, double t_offset, Vector* state_offset);

    /**
     * @brief Unimplemented default constructor.
     */
    DMD();

    /**
     * @brief Unimplemented copy constructor.
     */
    DMD(const DMD& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    DMD&
    operator = (
        const DMD& rhs);

    /**
     * @brief Internal function to multiply d_phi with the eigenvalues.
     */
    std::pair<Matrix*, Matrix*> phiMultEigs(double t, int deg = 0);

    /**
     * @brief Construct the DMD object.
     */
    void constructDMD(const Matrix* f_snapshots,
                      int rank,
                      int num_procs,
                      const Matrix* W0,
                      double linearity_tol);

    /**
     * @brief Returns a pair of pointers to the minus and plus snapshot matrices
     */
    virtual std::pair<Matrix*, Matrix*> computeDMDSnapshotPair(
        const Matrix* snapshots);

    /**
     * @brief Compute phi.
     */
    virtual void computePhi(DMDInternal dmd_internal_obj);

    /**
     * @brief Compute the appropriate exponential function when predicting the solution.
     */
    virtual std::complex<double> computeEigExp(std::complex<double> eig, double t);

    /**
     * @brief Add the appropriate offset when predicting the solution.
     */
    virtual void addOffset(Vector*& result, double t = 0.0, int deg = 0);

    /**
     * @brief Get the snapshot matrix contained within d_snapshots.
     */
    const Matrix* createSnapshotMatrix(std::vector<Vector*> snapshots);

    /**
     * @brief The rank of the process this object belongs to.
     */
    int d_rank;

    /**
     * @brief The number of processors being run on.
     */
    int d_num_procs;

    /**
     * @brief The total dimension of the sample vector.
     */
    int d_dim;

    /**
     * @brief The time step size between samples.
     */
    double d_dt = -1.0;

    /**
     * @brief The time offset of the first sample.
     */
    double d_t_offset;

    /**
     * @brief std::vector holding the snapshots.
     */
    std::vector<Vector*> d_snapshots;

    /**
     * @brief The stored times of each sample.
     */
    std::vector<Vector*> d_sampled_times;

    /**
     * @brief State offset in snapshot.
     */
    Vector* d_state_offset = NULL;

    /**
     * @brief Whether the DMD has been trained or not.
     */
    bool d_trained;

    /**
     * @brief Whether the initial condition has been projected.
     */
    bool d_init_projected;

    /**
     * @brief The maximum number of singular vectors.
     */
    int d_num_singular_vectors;

    /**
     * @brief std::vector holding the signular values.
     */
    std::vector<double> d_sv;

    /**
     * @brief The energy fraction used to obtain the DMD modes.
     */
    double d_energy_fraction;

    /**
     * @brief The number of columns used after obtaining the SVD decomposition.
     */
    int d_k;

    /**
     * @brief The left singular vector basis.
     */
    Matrix* d_basis = NULL;

    /**
     * @brief Whether to use the alternative basis for output, i.e. phi = U^(+)*V*Omega^(-1)*X.
     */
    bool d_alt_output_basis = false;

    /**
     * @brief A_tilde
     */
    Matrix* d_A_tilde = NULL;

    /**
     * @brief The real part of d_phi.
     */
    Matrix* d_phi_real = NULL;

    /**
     * @brief The imaginary part of d_phi.
     */
    Matrix* d_phi_imaginary = NULL;

    /**
     * @brief The real part of d_phi_squared_inverse.
     */
    Matrix* d_phi_real_squared_inverse = NULL;

    /**
     * @brief The imaginary part of d_phi_squared_inverse.
     */
    Matrix* d_phi_imaginary_squared_inverse = NULL;

    /**
     * @brief The real part of the projected initial condition.
     */
    Vector* d_projected_init_real = NULL;

    /**
     * @brief The imaginary part of the projected initial condition.
     */
    Vector* d_projected_init_imaginary = NULL;

    /**
     * @brief A vector holding the complex eigenvalues of the eigenmodes.
     */
    std::vector<std::complex<double>> d_eigs;

};

}

#endif
