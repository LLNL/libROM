/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description:Computes the DMD algorithm on the given snapshot matrix.

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
     * @param[in] f_snapshots The snapshot vectors for the RHS.
     * @param[in] energy_fraction The energy fraction to keep after doing SVD.
     * @param[in] rank The rank of this process.
     * @param[in] num_procs The total number of processes.
     */
    DMD(const Matrix* f_snapshots,
        double energy_fraction,
        int rank,
        int num_procs);

    /**
     * @brief Constructor.
     *
     * @param[in] f_snapshots The snapshot vectors for the RHS.
     * @param[in] k The number of modes (eigenvalues) to keep after doing SVD.
     * @param[in] rank The rank of this process.
     * @param[in] num_procs The total number of processes.
     */
    DMD(const Matrix* f_snapshots,
        int k,
        int rank,
        int num_procs);

    /**
     * @brief Predict new initial condition using d_phi.
     *
     * @param[in] init The initial condition.
     * @param[in] t The time of the outputted state.
     */
    std::pair<Vector*, Vector*> projectInitialCondition(const Vector* init);

    /**
     * @brief Predict state given a time. Uses the projected initial condition of the
     *        training dataset (the first column).
     *
     * @param[in] t The time of the outputted state.
     * @param[in] dt The delta time of the simulation.
     */
    Vector* predict(int t, int dt);

    /**
     * @brief Predict state given a new initial condition and time.
     *        The initial condition must be projected using projectInitialCondition
     *        for correct results.
     *
     * @param[in] init The initial condition.
     * @param[in] t The time of the outputted state.
     * @param[in] dt The delta time of the simulation.
     */
    Vector* predict(const std::pair<Vector*, Vector*> init, int t, int dt);

private:

    std::pair<Matrix*, Matrix*> phiMultEigs(int t, int dt);

    void constructDMD(const Matrix* f_snapshots,
                      int rank,
                      int num_procs);

    Matrix* d_phi_real;
    Matrix* d_phi_imaginary;
    Vector* d_projected_init_real;
    Vector* d_projected_init_imaginary;

    std::vector<std::complex<double>> d_eigs;

    double d_energy_fraction;

    int d_k;

};

/**
 * @brief Create a snapshot matrix from a vector of CAROM::Vector snapshots.
 *        The snapshot matrix will remain distributed if the snapshots are
          distributed.
 *
 * @param[in] snapshots The vector of CAROM::Vector snapshots
 */
const Matrix* createSnapshotMatrix(std::vector<Vector> snapshots);

}

#endif
