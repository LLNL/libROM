/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the DMD algorithm on the given snapshot matrix.

#ifndef included_DMD_h
#define included_DMD_h

#include <vector>
#include <complex>

namespace CAROM {

class Matrix;
class Vector;

class DMD
{
public:

    /**
     * @brief Constructor.
     *
     * @param[in] dim The full-order state dimension.
     */
    DMD(int dim);

    /**
     * @brief Sample the new state, u_in.
     *
     * @pre u_in != 0
     *
     * @param[in] u_in The state at the specified time.
     */
    void takeSample(double* u_in);

    /**
     * @param[in] energy_fraction The energy fraction to keep after doing SVD.
     */
    void train(double energy_fraction);

    /**
     * @param[in] k The number of modes (eigenvalues) to keep after doing SVD.
     */
    void train(int k);

    /**
     * @brief Predict new initial condition using d_phi.
     *
     * @param[in] init The initial condition.
     */
    void projectInitialCondition(const Vector* init);

    /**
     * @brief Predict state given a time. Uses the projected initial condition of the
     *        training dataset (the first column).
     *
     * @param[in] n The time of the outputted state (t/dt)
     */
    Vector* predict(double n);

    /**
     * @brief Predict state given a new initial condition and time.
     *        The initial condition must be projected using projectInitialCondition
     *        for correct results.
     *
     * @param[in] init The initial condition.
     * @param[in] n The time of the outputted state (t/dt)
     */
    Vector* predict(const std::pair<Vector*, Vector*> init, double n);

    /**
     * @brief Get the snapshot matrix contained within d_snapshots.
     */
    const Matrix* getSnapshotMatrix();

private:

    /**
     * @brief Unimplemented default constructor.
     */
    DMD();

    /**
     * @brief Unimplemented copy constructor.
     */
    DMD(
        const DMD& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    DMD&
    operator = (
        const DMD& rhs);

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
     * @brief std::vector holding the snapshots.
     */
    std::vector<Vector> d_snapshots;

    /**
     * @brief Internal function to multiply d_phi with the eigenvalues.
     */
    std::pair<Matrix*, Matrix*> phiMultEigs(double n);

    /**
     * @brief Internal function to obtain the DMD modes.
     */
    void constructDMD(const Matrix* f_snapshots,
                      int rank,
                      int num_procs);

    /**
     * @brief The real part of d_phi.
     */
    Matrix* d_phi_real;

    /**
     * @brief The imaginary part of d_phi.
     */
    Matrix* d_phi_imaginary;

    /**
     * @brief The real part of the projected initial condition.
     */
    Vector* d_projected_init_real;

    /**
     * @brief The imaginary part of the projected initial condition.
     */
    Vector* d_projected_init_imaginary;

    /**
     * @brief A vector holding the complex eigenvalues of the eigenmodes.
     */
    std::vector<std::complex<double>> d_eigs;

    /**
     * @brief The energy fraction used to obtain the DMD modes.
     */
    double d_energy_fraction;

    /**
     * @brief The number of columns used after obtaining the SVD decomposition.
     */
    int d_k;

};

}

#endif
