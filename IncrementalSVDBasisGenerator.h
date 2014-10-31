/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The concrete wrapper class for a specific incremental SVD
 *              algorithm and time stepper.  Implements interface of
 *              SVDBasisGenerator.
 *
 *****************************************************************************/

#ifndef included_IncrementalSVDBasisGenerator_h
#define included_IncrementalSVDBasisGenerator_h

#include "SVDBasisGenerator.h"
#include "IncrementalSVDTimeStepper.h"

namespace CAROM {

/**
 * Class IncrementalSVDBasisGenerator implements the interface of base class
 * SVDBasisGenerator for the incremental svd algorithm.  Either the fast update
 * or the naive incremental algorithm may be specified through the constructor.
 */
class IncrementalSVDBasisGenerator : public SVDBasisGenerator
{
   public:
      /**
       * @brief Constructor.
       *
       * @pre dim > 0
       * @pre redundancy_tol > 0.0
       * @pre increments_per_time_interval > 0
       * @pre sampling_tol > 0.0
       * @pre max_time_between_snapshots > 0.0
       *
       * @param[in] dim The dimension of the system on this processor.
       * @param[in] redundancy_tol Tolerance to determine if a sample is
       *                           redundant or not.
       * @param[in] skip_redundant If true skip redundant samples.
       * @param[in] increments_per_time_interval The maximum number of samples
       *                                         in each time interval.
       * @param[in] sampling_tol Time step control tolerance.  Limits error in
       *                         projection of solution into reduced order
       *                         space.
       * @param[in] max_time_between_snapshots Hard upper bound on time step.
       * @param[in] fast_update If true use the fast update algorithm.
       * @param[in] basis_file_name The base part of the name of the file
       *                            containing the basis vectors.  Each process
       *                            will append its process ID to this base
       *                            name.
       * @param[in] file_format The format of the file containing the basis
       *                        vectors.
       */
      IncrementalSVDBasisGenerator(
         int dim,
         double redundancy_tol,
         bool skip_redundant,
         int increments_per_time_interval,
         double sampling_tol,
         double max_time_between_snapshots,
         bool fast_update,
         const std::string& basis_file_name,
         Database::formats file_format = Database::HDF5);

      /**
       * @brief Destructor.
       */
      virtual
      ~IncrementalSVDBasisGenerator();

      /**
       * @brief Returns true if it is time for the next svd snapshot.
       *
       * @pre time >= 0
       *
       * @param[in] time Time of interest.
       *
       * @return True if it is time for the next snapshot to be taken.
       */
      virtual
      bool
      isNextSnapshot(
         double time);

      /**
       * @brief Add a snapshot to the incremental svd at the given time.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] time The simulation time for the state.
       */
      virtual
      void
      takeSnapshot(
         double* u_in,
         double time);

      /**
       * @brief Computes next time an svd snapshot is needed.
       *
       * @pre u_in != 0
       * @pre rhs_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] rhs_in The right hand side at the specified time.
       * @param[in] time The simulation time for the state.
       */
      virtual
      double
      computeNextSnapshotTime(
         double* u_in,
         double* rhs_in,
         double time);

      /**
       * @brief Returns the basis vectors for the current time interval as a
       * Matrix.
       *
       * @return The basis vectors for the current time interval.
       */
      virtual
      const Matrix*
      getBasis();

      /**
       * @brief Returns the number of time intervals on which different sets of
       * basis vectors are defined.
       *
       * @return The number of time intervals on which there are basis vectors.
       */
      virtual
      int
      getNumBasisTimeIntervals() const;

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
      virtual
      double
      getBasisIntervalStartTime(
         int which_interval) const;

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      IncrementalSVDBasisGenerator();

      /**
       * @brief Unimplemented copy constructor.
       */
      IncrementalSVDBasisGenerator(
         const IncrementalSVDBasisGenerator& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      IncrementalSVDBasisGenerator&
      operator = (
         const IncrementalSVDBasisGenerator& rhs);

      /**
       * @brief Pointer to the underlying time step control object.
       */
      boost::shared_ptr<IncrementalSVDTimeStepper> d_isvdts;
};

}

#endif
