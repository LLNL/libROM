/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The wrapper class for an SVD algorithm and
//              sampler.  This class controls all aspects of basis
//              vector generation.

#ifndef included_BasisGenerator_h
#define included_BasisGenerator_h

#include "BasisWriter.h"
#include "BasisReader.h"
#include "Options.h"
#include "svd/SVD.h"

#include "mpi.h"

#include <cmath>

/* Use C++11 built-in shared pointers if available; else fallback to Boost. */
#if __cplusplus >= 201103L
#include <memory>
#else
#include <boost/shared_ptr.hpp>
#endif

#include <string.h>

namespace CAROM {

class BasisWriter;
class BasisReader;
class Matrix;

/**
 * Class BasisGenerator defines the interface for the generation of basis
 * vectors via the SVD method.  This class wraps the SVD algorithm and sampler
 * and controls all aspects of basis vector generation.
 */
class BasisGenerator
{
public:
    /**
     * @brief Constructor.
     *
     * @param[in] options The struct containing the options for this basis
     *                    generator.
     * @param[in] incremental Whether to conduct static or incremental SVD
     * @param[in] basis_file_name The base part of the name of the file
     *                            containing the basis vectors.  Each process
     *                            will append its process ID to this base
     *                            name.
     * @param[in] file_format The format of the file containing the basis
     *                        vectors.
     */
    BasisGenerator(
        Options options,
        bool incremental,
        const std::string& basis_file_name = "",
        Database::formats file_format = Database::HDF5);

    /**
     * @brief Destructor.
     */
    ~BasisGenerator();

    /**
     * @brief Returns true if it is time for the next svd sample.
     *
     * @pre time >= 0.0
     *
     * @param[in] time Time of interest.
     *
     * @return True if it is time for the next sample to be taken.
     */
    bool
    isNextSample(
        double time);

    /**
     * @brief Check whether right basis vectors will be updated.
     *
     * * @return True if the right basis vectors will be updated.
     */
    bool
    updateRightSV()
    {
        return d_update_right_SV;
    }

    /**
     * @brief Sample the new state, u_in, at the given time.
     *
     * @pre u_in != 0
     * @pre time >= 0.0
     *
     * @param[in] u_in The state at the specified time.
     * @param[in] time The simulation time for the state.
     * @param[in] dt The current simulation dt.
     * @param[in] add_without_increase If true, the addLinearlyDependent is
     *                                 invoked. This only applies to incremental
     *                                 SVD.
     *
     * @return True if the sampling was successful.
     */
    bool
    takeSample(
        double* u_in,
        double time,
        double dt,
        bool add_without_increase = false);

    /**
     * @brief Signal that the final sample has been taken.
     *
     * @param[in] kind A string equal to "basis" or "snapshot", representing
     *                 which one will be written.
     */
    void
    endSamples(const std::string& kind = "basis")
    {
        if (d_basis_writer) {
            d_basis_writer->writeBasis(kind);
        }
    }

    /**
     * @brief Write current snapshot matrix.
     */
    void
    writeSnapshot()
    {
        if (d_basis_writer) {
            d_basis_writer->writeBasis("snapshot");
        }
    }

    /**
     * @brief Load previously saved sample (basis or state).
     *
     * @param[in] base_file_name The base part of the name of the files
     *                           holding the basis / snapshot vectors.
     * @param[in] kind A string equal to "basis" or "snapshot", representing
     *                 which kind of data to load.
     * @param[in] cut_off The maximum number of bases or snapshots to read.
     * @param[in] db_format Format of the file to read.
     */
    void
    loadSamples(const std::string& base_file_name,
                const std::string& kind  = "basis",
                int cut_off = 1e9,
                Database::formats db_format = Database::HDF5);

    /**
     * @brief Computes next time an svd sample is needed.
     *
     * @pre u_in != 0
     * @pre rhs_in != 0
     * @pre time >= 0.0
     *
     * @param[in] u_in The state at the specified time.
     * @param[in] rhs_in The right hand side at the specified time.
     * @param[in] time The simulation time for the state.
     */
    double
    computeNextSampleTime(
        double* u_in,
        double* rhs_in,
        double time);

    /**
     * @brief Returns the basis vectors for the current time interval as a
     * Matrix.
     *
     * @return The basis vectors for the current time interval.
     */
    const Matrix*
    getSpatialBasis()
    {
        return d_svd->getSpatialBasis();
    }

    /**
     * @brief Returns the temporal basis vectors for the current time interval as a
     * Matrix.
     *
     * @return The temporal basis vectors for the current time interval.
     */
    const Matrix*
    getTemporalBasis()
    {
        return d_svd->getTemporalBasis();
    }

    /**
     * @brief Returns the singular values for the current time interval as a
     * Vector.
     *
     * @return The singular values for the current time interval.
     */
    const Vector*
    getSingularValues()
    {
        return d_svd->getSingularValues();
    }

    /**
     * @brief Returns the snapshot matrix for the current time interval.
     *
     * @return The snapshot matrix for the current time interval.
     */
    const Matrix*
    getSnapshotMatrix()
    {
        return d_svd->getSnapshotMatrix();
    }

    /**
     * @brief Returns the number of time intervals on which different sets of
     * basis vectors are defined.
     *
     * @return The number of time intervals on which there are basis vectors.
     */
    int
    getNumBasisTimeIntervals() const
    {
        return d_svd->getNumBasisTimeIntervals();
    }

    /**
     * @brief Returns the start time for the requested time interval.
     *
     * @pre 0 <= which_interval
     * @pre which_interval < getNumBasisTimeIntervals()
     *
     * @param[in] which_interval Time interval whose start time is needed.
     *
     * @return The start time for the requested time interval.
     */
    double
    getBasisIntervalStartTime(
        int which_interval) const
    {
        CAROM_VERIFY(0 <= which_interval);
        CAROM_VERIFY(which_interval < getNumBasisTimeIntervals());
        return d_svd->getBasisIntervalStartTime(which_interval);
    }

    /**
     * @brief Returns the number of samples taken.
     *
     * @return The number of samples taken.
     */
    int getNumSamples() const
    {
        return d_svd->getNumSamples();
    }

protected:
    /**
     * @brief Writer of basis vectors.
     */
    BasisWriter* d_basis_writer;

    /**
     * @brief Reader of basis vectors.
     */
    BasisReader* d_basis_reader;

    /**
     * @brief Whether to write snapshots instead of bases.
     */
    bool d_write_snapshots;

    /**
     * @brief Pointer to the abstract SVD algorithm object.
     */
#if __cplusplus >= 201103L
    std::shared_ptr<SVD> d_svd;
#else
    boost::shared_ptr<SVD> d_svd;
#endif

private:
    /**
     * @brief Unimplemented default constructor.
     */
    BasisGenerator();
    /**
     * @brief Unimplemented copy constructor.
     */
    BasisGenerator(
        const BasisGenerator& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    BasisGenerator&
    operator = (
        const BasisGenerator& rhs);

    /**
     * @brief Resets sample time step.
     *
     * @param[in] new_dt New value of sample time step.
     */
    void
    resetDt(
        double new_dt);

    /**
     * @brief Returns the dimension of the system on this processor.
     *
     * @return The dimension of the system on this processor.
     */
    int
    getDim()
    {
        return d_svd->getDim();
    }

    /**
     * @brief If using incremental or static SVD
     */
    bool d_incremental;

    /**
     * @brief if true, isNextSample returns always true
     */
    bool d_update_right_SV;

    /**
     * @brief Sampling control tolerance.
     *
     * Limits error in projection of solution into the reduced order space.
     */
    double d_tol;

    /**
     * @brief Maximum time between samples.
     */
    double d_max_time_between_samples;

    /**
     * @brief Minimum sampling time step scale factor.
     */
    double d_min_sampling_time_step_scale;

    /**
     * @brief Sampling time step scale factor to apply to algorithm.
     */
    double d_sampling_time_step_scale;

    /**
     * @brief Maximum sampling time step scale factor.
     */
    double d_max_sampling_time_step_scale;

    /**
     * @brief Current time step.
     */
    double d_dt;

    /**
     * @brief Next time at which a sample should be taken.
     */
    double d_next_sample_time;

    /**
     * @brief The number of processors being run on.
     */
    int d_num_procs;
};

}

#endif
