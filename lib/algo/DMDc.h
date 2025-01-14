/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the DMDc algorithm on the given snapshot matrix and control matrix. The
//              implemented dynamic mode decomposition with control algorithm is derived from
//              Proctor et. al's paper "Dynamic mode decomposition with control":
//              https://arxiv.org/abs/1409.6358

#ifndef included_DMDc_h
#define included_DMDc_h

#include "ParametricDMDc.h"
#include <vector>
#include <complex>

namespace CAROM {

class Matrix;
class Vector;
class ComplexEigenPair;

/**
 * Class DMDc implements the DMDc algorithm on a given snapshot matrix.
 */
class DMDc
{
public:

    /**
     * @brief Constructor. Basic DMDc with uniform time step size.
     *
     * @param[in] dim              The full-order state dimension.
     * @param[in] dim_c            The control dimension.
     * @param[in] dt               The dt between samples.
     * @param[in] state_offset     The state offset.
     */
    DMDc(int dim, int dim_c, double dt,
         std::shared_ptr<Vector> state_offset = nullptr);

    /**
     * @brief Constructor. DMDc from saved models.
     *
     * @param[in] base_file_name The base part of the filename of the
     *                           database to load when restarting from a save.
     */
    DMDc(std::string base_file_name);

    /**
     * @brief Destroy the DMDc object
     */
    virtual ~DMDc() {};

    /**
     * @brief Set the state offset.
     */
    virtual void setOffset(std::shared_ptr<Vector> & offset_vector);

    /**
     * @brief Sample the new state, u_in. Any samples in d_snapshots
     *        taken at the same or later time will be erased.
     *
     * @pre u_in != 0
     * @pre t >= 0.0
     *
     * @param[in] u_in      The new state.
     * @param[in] t         The time of the newly sampled state.
     * @param[in] f_in      The control.
     * @param[in] last_step Whether it is the last step.
     */
    virtual void takeSample(double* u_in, double t, double* f_in,
                            bool last_step = false);

    /**
     * @brief Train the DMDc model with energy fraction criterion.
     *        The control matrix B may be available and used in training.
     *        It is yet to be implemented and requires
     *        consideration in sparse matrix multiplication in general.
     *        (See Section III.B. in https://arxiv.org/pdf/1409.6358.pdf)
     *        In default, the control matrix is unknown, with B = NULL.
     *
     * @param[in] energy_fraction The energy fraction to keep after doing SVD.
     * @param[in] B               (Optional) The control matrix B.
     */
    virtual void train(double energy_fraction, const Matrix* B = NULL);

    /**
     * @brief Train the DMDc model with specified reduced dimension.
     *        The control matrix B may be available and used in training.
     *        It is yet to be implemented and requires
     *        consideration in sparse matrix multiplication in general.
     *        (See Section III.B. in https://arxiv.org/pdf/1409.6358.pdf)
     *        In default, the control matrix is unknown, with B = NULL.
     *
     * @param[in] k               The number of modes to keep after doing SVD.
     * @param[in] B               (Optional) The control matrix B.
     */
    virtual void train(int k, const Matrix* B = NULL);

    /**
     * @brief Project U using d_phi, where U is the initial condition and the controls.
     *        Calculate pinv(phi) x U, or more precisely,
     *        (phi* x phi)^{-1} x phi* x U, where phi* is the conjugate transpose.
     *
     * @param[in] init     The initial condition.
     * @param[in] controls The controls.
     * @param[in] t_offset The initial time offset.
     */
    void project(const Vector & init, const Matrix & controls,
                 double t_offset = -1.0);

    /**
     * @brief Predict state given a time. Uses the projected initial condition of the
     *        training dataset (the first column).
     *
     * @param[in] t   The time of the output state
     */
    std::unique_ptr<Vector> predict(double t);

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
    std::unique_ptr<const Matrix> getSnapshotMatrix();

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
     * @brief Output the DMDc record in CSV files.
     */
    void summary(std::string base_file_name);

protected:
    /**
     * @brief Obtain DMD model interpolant at desired parameter point by
     *        interpolation of DMD models from training parameter points.
     *
     * @param[in] parametric_dmdc       The interpolant DMD model at the desired point.
     * @param[in] parameter_points      The training parameter points.
     * @param[in] dmdcs                 The DMD objects associated with
     *                                  each training parameter point.
     * @param[in] controls              The matrices of controls from previous
     *                                  runs which we use to interpolate.
     * @param[in] controls_interpolated The interpolated controls.
     * @param[in] desired_point         The desired point at which to create a parametric DMD.
     * @param[in] rbf                   The RBF type ("G" == gaussian,
     *                                  "IQ" == inverse quadratic,
     *                                  "IMQ" == inverse multiquadric)
     * @param[in] interp_method         The interpolation method type
     *                                  ("LS" == linear solve,
     *                                  "IDW" == inverse distance weighting,
     *                                  "LP" == lagrangian polynomials)
     * @param[in] closest_rbf_val       The RBF parameter determines the width of influence.
     *                                  Set the RBF value of the nearest two parameter points
     *                                  to a value between 0.0 to 1.0
     * @param[in] reorthogonalize_W     Whether to reorthogonalize the interpolated W (basis) matrix.
     */
    friend void getParametricDMDc<DMDc>(std::unique_ptr<DMDc>& parametric_dmdc,
                                        const std::vector<Vector>& parameter_points,
                                        std::vector<DMDc*>& dmdcs,
                                        std::vector<std::shared_ptr<Matrix>> & controls,
                                        std::shared_ptr<Matrix> & controls_interpolated,
                                        const Vector & desired_point,
                                        std::string rbf,
                                        std::string interp_method,
                                        double closest_rbf_val,
                                        bool reorthogonalize_W);

    /**
     * @brief Constructor. Variant of DMDc with non-uniform time step size.
     *
     * @param[in] dim               The full-order state dimension.
     * @param[in] dim_c             The control dimension.
     * @param[in] state_offset      The state offset.
     */
    DMDc(int dim, int dim_c, std::shared_ptr<Vector> state_offset = nullptr);

    /**
     * @brief Constructor. Specified from DMDc components.
     *
     * @param[in] eigs d_eigs
     * @param[in] phi_real d_phi_real
     * @param[in] phi_imaginary d_phi_imaginary
     * @param[in] B_tilde d_B_tilde
     * @param[in] k d_k
     * @param[in] dt d_dt
     * @param[in] t_offset d_t_offset
     * @param[in] state_offset d_state_offset
     * @param[in] basis d_basis, set by DMDc class for offline stages. When
     *            interpolating a new DMDc, we enter the interpolated basis
     *            explicitly.
     */
    DMDc(std::vector<std::complex<double>> & eigs,
         std::shared_ptr<Matrix> & phi_real,
         std::shared_ptr<Matrix> & phi_imaginary,
         std::shared_ptr<Matrix> & B_tilde,
         int k, double dt, double t_offset,
         std::shared_ptr<Vector> & state_offset,
         std::shared_ptr<Matrix> & basis);

    /**
       * @brief Unimplemented default constructor.
       */
    DMDc();

    /**
     * @brief Unimplemented copy constructor.
     */
    DMDc(const DMDc& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    DMDc&
    operator = (
        const DMDc& rhs);

    /**
     * @brief Internal function to multiply d_phi with the eigenvalues.
     */
    std::pair<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>>
            phiMultEigs(double t);

    /**
     * @brief Construct the DMDc object.
     */
    void constructDMDc(const Matrix & f_snapshots,
                       const Matrix & f_controls,
                       int rank,
                       int num_procs,
                       const Matrix* B);

    /**
     * @brief Returns a pair of pointers to the minus and plus snapshot matrices
     */
    virtual std::pair<Matrix*, Matrix*> computeDMDcSnapshotPair(
        const Matrix & snapshots, const Matrix & controls, const Matrix* B);

    /**
     * @brief Compute the appropriate exponential function when predicting the solution.
     */
    virtual std::complex<double> computeEigExp(std::complex<double> eig, double t);

    /**
     * @brief Add the state offset when predicting the solution.
     */
    virtual void addOffset(Vector & result);

    /**
     * @brief Get the snapshot matrix contained within d_snapshots.
     */
    std::unique_ptr<const Matrix>
    createSnapshotMatrix(const std::vector<std::shared_ptr<Vector>> & snapshots);

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
     * @brief The total dimension of the control vector.
     */
    int d_dim_c;

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
    std::vector<std::shared_ptr<Vector>> d_snapshots;

    /**
     * @brief std::vector holding the controls.
     */
    std::vector<std::shared_ptr<Vector>> d_controls;

    /**
     * @brief The stored times of each sample.
     */
    std::vector<double> d_sampled_times;

    /**
     * @brief State offset in snapshot.
     */
    std::shared_ptr<Vector> d_state_offset;

    /**
     * @brief Whether the DMDc has been trained or not.
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
     * @brief The energy fraction used to obtain the DMDc modes.
     */
    double d_energy_fraction;

    /**
     * @brief The number of columns used after obtaining the SVD decomposition.
     */
    int d_k;

    /**
     * @brief The left singular vector basis.
     */
    std::shared_ptr<Matrix> d_basis;

    /**
     * @brief A_tilde
     */
    std::shared_ptr<Matrix> d_A_tilde;

    /**
     * @brief B_tilde
     */
    std::shared_ptr<Matrix> d_B_tilde;

    /**
     * @brief The real part of d_phi.
     */
    std::shared_ptr<Matrix> d_phi_real;

    /**
     * @brief The imaginary part of d_phi.
     */
    std::shared_ptr<Matrix> d_phi_imaginary;

    /**
     * @brief The real part of d_phi_squared_inverse.
     */
    std::shared_ptr<Matrix> d_phi_real_squared_inverse;

    /**
     * @brief The imaginary part of d_phi_squared_inverse.
     */
    std::shared_ptr<Matrix> d_phi_imaginary_squared_inverse;

    /**
     * @brief The real part of the projected initial condition.
     */
    std::shared_ptr<Vector> d_projected_init_real;

    /**
     * @brief The imaginary part of the projected initial condition.
     */
    std::shared_ptr<Vector> d_projected_init_imaginary;

    /**
     * @brief The real part of the projected controls.
     */
    std::unique_ptr<Matrix> d_projected_controls_real;

    /**
     * @brief The imaginary part of the projected controls.
     */
    std::unique_ptr<Matrix> d_projected_controls_imaginary;

    /**
     * @brief A vector holding the complex eigenvalues of the eigenmodes.
     */
    std::vector<std::complex<double>> d_eigs;

};

}

#endif
