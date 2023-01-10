/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: An abstract class defining the interface to the generic SVD
//              algorithm.

#ifndef included_SVD_h
#define included_SVD_h

#include "linalg/Matrix.h"
#include "linalg/Options.h"
#include <vector>

namespace CAROM {

/**
 * Class SVD defines the interface to the generic SVD algorithm.  The API is
 * intentionally small.  One may collect the samples, compute the SVD, and get
 * the dimension of the system.
 */
class SVD
{
public:
    /**
     * @brief Constructor.
     *
     * @param[in] options The struct containing the options for this abstract
     * SVD class.
     * @see Options
     */
    SVD(
        Options options);

    /**
     * Destructor.
     */
    ~SVD();

    /**
     * @brief Collect the new sample, u_in at supplied time.
     *
     * @pre u_in != 0
     * @pre time >= 0.0
     *
     * @param[in] u_in The new sample.
     * @param[in] time The simulation time of the new sample.
     * @param[in] add_without_increase If true, the addLinearlyDependent is invoked.
     *                                 This only applies to incremental SVD.
     *
     * @return True if the sampling was successful.
     */
    virtual
    bool
    takeSample(
        double* u_in,
        double time,
        bool add_without_increase) = 0;

    /**
     * @brief Returns the dimension of the system on this processor.
     *
     * @return The dimension of the system on this processor.
     */
    int
    getDim() const
    {
        return d_dim;
    }

    /**
     * @brief Returns the basis vectors for the current time interval.
     *
     * @return The basis vectors for the current time interval.
     */
    virtual
    const Matrix*
    getSpatialBasis() = 0;

    /**
     * @brief Returns the temporal basis vectors for the current time interval.
     *
     * @return The temporal basis vectors for the current time interval.
     */
    virtual
    const Matrix*
    getTemporalBasis() = 0;

    /**
     * @brief Returns the singular values for the current time interval.
     *
     * @return The singular values for the current time interval.
     */
    virtual
    const Vector*
    getSingularValues() = 0;

    /**
     * @brief Returns the singular values for the current time interval.
     *
     * @return The singular values for the current time interval.
     */
    virtual
    const Matrix*
    getSnapshotMatrix() = 0;

    /**
     * @brief Returns the number of time intervals on which different sets
     * of basis vectors are defined.
     *
     * @return The number of time intervals on which there are basis vectors.
     */
    int
    getNumBasisTimeIntervals() const
    {
        return static_cast<int>(d_time_interval_start_times.size());
    }

    /**
     * @brief Returns the start time for the requested time interval.
     *
     * @pre 0 <= which_interval
     * @pre which_interval < getNumBasisTimeIntervals()
     *
     * @param[in] which_interval The time interval of interest.
     *
     * @return The start time for the requested time interval.
     */
    double
    getBasisIntervalStartTime(
        int which_interval) const
    {
        CAROM_VERIFY(0 <= which_interval);
        CAROM_VERIFY(which_interval < getNumBasisTimeIntervals());

        std::size_t i = static_cast<std::size_t>(which_interval);
        return d_time_interval_start_times[i];
    }

    /**
     * @brief Returns true if the next sample will result in a new time
     * interval.
     *
     * @return True if the next sample results in the creation of a new time
     *         interval.
     */
    bool
    isNewTimeInterval() const
    {
        return (d_num_samples == 0) ||
               (d_num_samples >= d_samples_per_time_interval);
    }

    /**
     * @brief Increase the number of time intervals by one
     *
     */
    void
    increaseTimeInterval()
    {
        int num_time_intervals =
            static_cast<int>(d_time_interval_start_times.size());
        CAROM_VERIFY(d_max_time_intervals == -1 ||
                     num_time_intervals < d_max_time_intervals);
        d_time_interval_start_times.resize(
            static_cast<unsigned>(num_time_intervals) + 1);
    }

    /**
     * @brief Get the number of samples taken.
     *
     */
    int getNumSamples() const
    {
        return d_num_samples;
    }

protected:
    /**
     * @brief Dimension of the system.
     */
    int d_dim;

    /**
     * @brief Number of samples stored for the current time interval.
     */
    int d_num_samples;

    /**
     * @brief Number of rows in right singular matrix.
     */
    int d_num_rows_of_W;

    /**
     * @brief The maximum number of samples to be collected for a time
     * interval.
     */
    const int d_samples_per_time_interval;

    /**
     * @brief The maximum number of time intervals.
     */
    const int d_max_time_intervals;

    /**
     * @brief The globalized basis vectors for the current time interval.
     *
     * The basis vectors are large and each process owns all of the basis
     * vectors.
     */
    Matrix* d_basis;

    /**
     * @brief The globalized right basis vectors for the current time interval.
     *
     * Depending on the SVD algorithm, it may be  distributed across all
     * processors or each processor may own all of U.
     */
    Matrix* d_basis_right;

    /**
     * @brief The matrix U which is large.
     *
     * Depending on the SVD algorithm, d_U may be  distributed across all
     * processors or each processor may own all of U.
     */
    Matrix* d_U;

    /**
     * @brief The matrix U which is large.
     *
     * Depending on the SVD algorithm, d_W may be  distributed across all
     * processors or each processor may own all of U.
     */
    Matrix* d_W;

    /**
     * @brief The vector S which is small.
     *
     * For all SVD algorithms, S is not distributed and the entire vector
     * exists on each processor.
     */
    Vector* d_S;

    /**
     * @brief The globalized snapshot vectors for the current time interval.
     *
     * The snapshot vectors are large and each process owns all of the snapshot
     * vectors.
     */
    Matrix* d_snapshots;

    /**
     * @brief The simulation time at which each time interval starts.
     */
    std::vector<double> d_time_interval_start_times;

    /**
     * @brief Flag to indicate if results of algorithm should be printed for
     * debugging purposes.
     */
    bool d_debug_algorithm;

private:
    /**
     * @brief Unimplemented default constructor.
     */
    SVD();

    /**
     * @brief Unimplemented copy constructor.
     */
    SVD(
        const SVD& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    SVD&
    operator = (
        const SVD& rhs);
};

}

#endif
