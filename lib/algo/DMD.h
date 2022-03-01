/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
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

#include <vector>
#include <complex>

namespace CAROM {

class Matrix;
class Vector;

/**
 * Class DMD implements the DMD algorithm on a given snapshot matrix.
 */
class DMD
{
public:

    /**
     * @brief Constructor.
     *
     * @param[in] dim        The full-order state dimension.
     * @param[in] dt         The dt between samples.
     */
    DMD(int dim, double dt);

    /**
     * @brief Constructor.
     *
     * @param[in] base_file_name The base part of the filename of the
     *                           database to load when restarting from a save.
     */
    DMD(std::string base_file_name);

    /**
     * @brief Sample the new state, u_in.
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
     */
    virtual void train(double energy_fraction);

    /**
     * @param[in] k The number of modes (eigenvalues) to keep after doing SVD.
     */
    virtual void train(int k);

    /**
     * @brief Output the DMD record in CSV files.
     */
    void summary(std::string output_path);

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
     * @param[in] t The time of the outputted state
     */
    Vector* predict(double t);

    /**
     * @brief Predict state given a new initial condition and time.
     *        The initial condition must be projected using projectInitialCondition
     *        for correct results.
     *
     * @param[in] init The initial condition.
     * @param[in] t The time of the outputted state
     */
    Vector* predict(const std::pair<Vector*, Vector*> init, double t);

    /**
     * @brief Get the time offset contained within d_t_offset.
     */
    double getTimeOffset() const;

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
     * @brief Save the object state to a file.
     *
     * @param[in] base_file_name The base part of the filename to save the
     *                           database to.
     */
    virtual void save(std::string base_file_name);

protected:

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
     * @brief Internal function to multiply d_phi with the eigenvalues.
     */
    std::pair<Matrix*, Matrix*> phiMultEigs(double t);

    /**
     * @brief Internal function to obtain the DMD modes.
     */
    void constructDMD(const Matrix* f_snapshots,
                      int rank,
                      int num_procs);

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
    double d_dt;

    /**
     * @brief The time offset of the first sample.
     */
    double d_t_offset;

    /**
     * @brief std::vector holding the snapshots.
     */
    std::vector<Vector*> d_snapshots;

    /**
     * @brief Whether the DMD has been trained or not.
     */
    bool d_trained;

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

};

}

#endif
