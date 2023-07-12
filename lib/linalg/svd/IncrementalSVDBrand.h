/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The concrete implementation of the incremental SVD algorithm
//              using Matthew Brand's (exact) "fast update" method.

#ifndef included_IncrementalSVDBrand_h
#define included_IncrementalSVDBrand_h

#include "IncrementalSVD.h"
#include "linalg/Options.h"

namespace CAROM {

/**
 * Class IncrementalSVDBrand implements Brand's fast update incremental SVD
 * algorithm by implementing the pure virtual methods of the IncrementalSVD
 * base class.
 */
class IncrementalSVDBrand : public IncrementalSVD
{
public:
    /**
     * @brief Destructor.
     */
    ~IncrementalSVDBrand();

    /**
     * @brief Returns the basis vectors for the current time interval as a
     *        Matrix.
     *
     * @return The basis vectors for the current time interval.
     */
    const Matrix*
    getSpatialBasis() override;

    /**
     * @brief Returns the temporal basis vectors for the current time interval
     *        as a Matrix.
     *
     * @return The temporal basis vectors for the current time interval.
     */
    const Matrix*
    getTemporalBasis() override;

private:
    friend class BasisGenerator;

    /**
     * @brief Constructor.
     *
     * @param[in] options The struct containing the options for this SVD
     *                    implementation.
     * @param[in] basis_file_name The base part of the name of the file
     *                            containing the basis vectors.  Each process
     *                            will append its process ID to this base
     *                            name.
     * @see Options
     */
    IncrementalSVDBrand(
        Options options,
        const std::string& basis_file_name);

    /**
     * @brief Unimplemented default constructor.
     */
    IncrementalSVDBrand();

    /**
     * @brief Unimplemented copy constructor.
     */
    IncrementalSVDBrand(
        const IncrementalSVDBrand& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    IncrementalSVDBrand&
    operator = (
        const IncrementalSVDBrand& rhs);

    /**
     * @brief Constructs the first SVD.
     *
     * @pre u != 0
     * @pre time >= 0.0
     *
     * @param[in] u The first state.
     * @param[in] time The simulation time for the first state.
     */
    virtual
    void
    buildInitialSVD(
        double* u,
        double time);

    /**
     * @brief Adds the new sampled the state vector, u, to the system.
     *
     * @pre u != 0
     *
     * @param[in] u The new state.
     * @param[in] add_without_increase If true, addLinearlyDependent is invoked.
     *
     * @return True if building the incremental SVD was successful.
     */
    bool
    buildIncrementalSVD(
        double* u, bool add_without_increase = false) override;

    /**
     * @brief Update the current spatial basis vectors.
     */
    void
    updateSpatialBasis();

    /**
     * @brief Update the current temporal basis vectors.
     */
    void
    updateTemporalBasis();

    /**
     * @brief Computes the current basis vectors.
     */
    virtual
    void
    computeBasis();

    /**
     * @brief Add a linearly dependent sample to the SVD.
     *
     * @pre A != 0
     * @pre sigma != 0
     *
     * @param[in] A The left singular vectors.
     * @param[in] W The right singular vectors.
     * @param[in] sigma The singular values.
     */
    void
    addLinearlyDependentSample(
        const Matrix* A,
        const Matrix* W,
        const Matrix* sigma);

    /**
     * @brief Add a new, unique sample to the SVD.
     *
     * @pre j != 0
     * @pre A != 0
     * @pre W != 0
     * @pre sigma != 0
     *
     * @param[in] j The new column of d_U.
     * @param[in] A The left singular vectors.
     * @param[in] W The right singular vectors.
     * @param[in] sigma The singular values.
     */
    void
    addNewSample(
        const Vector* j,
        const Matrix* A,
        const Matrix* W,
        Matrix* sigma);

    /**
     * @brief The matrix U'. U' is not distributed and the entire matrix
     *        exists on each processor.
     */
    Matrix* d_Up;

    /**
     * @brief The tolerance value used to remove small singular values.
     */
    double d_singular_value_tol;

};

}

#endif
