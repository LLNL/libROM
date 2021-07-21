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
     * @brief Predict state given a time.
     *
     * @param[in] t The time of the outputted state.
     */
    Vector predict(int t);

private:

    void constructDMD(const Matrix* f_snapshots,
                      int rank,
                      int num_procs);

    Matrix* d_phi;

    double d_energy_fraction;

    int d_k;

};

}

#endif
