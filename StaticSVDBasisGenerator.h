/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The concrete wrapper class for static SVD algorithm and time
 *              stepper.  Implements interface of SVDBasisGenerator.
 *
 *****************************************************************************/

#ifndef included_StaticSVDBasisGenerator_h
#define included_StaticSVDBasisGenerator_h

#include "SVDBasisGenerator.h"
#include "StaticSVDTimeStepper.h"

namespace CAROM {

/**
 * Class StaticSVDBasisGenerator implements the interface of base class
 * SVDBasisGenerator for the static svd algorithm.
 */
class StaticSVDBasisGenerator : public SVDBasisGenerator
{
   public:
      /**
       * @brief Constructor.
       *
       * @pre dim > 0
       * @pre increments_per_time_interval > 0
       *
       * @param[in] dim The dimension of the system on this processor.
       * @param[in] increments_per_time_interval The maximum number of samples
       *                                         in each time interval.
       * @param[in] basis_file_name The base part of the name of the file
       *                            containing the basis vectors.  Each process
       *                            will append its process ID to this base
       *                            name.
       * @param[in] file_format The format of the file containing the basis
       *                        vectors.
       */
      StaticSVDBasisGenerator(
         int dim,
         int increments_per_time_interval,
         const std::string& basis_file_name,
         Database::formats file_format = Database::HDF5);

      /**
       * @brief Destructor.
       */
      virtual
      ~StaticSVDBasisGenerator();

      /**
       * @brief Returns true if it is time for the next svd snapshot.
       *
       * @pre time >= 0
       *
       * @parm[in] time Time of interest.
       *
       * @return True if it is time for the next snapshot to be taken.
       */
      virtual
      bool
      isNextSnapshot(
         double time);

      /**
       * @brief Add a snapshot to the static svd at the given time.
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
      StaticSVDBasisGenerator();

      /**
       * @brief Unimplemented copy constructor.
       */
      StaticSVDBasisGenerator(
         const StaticSVDBasisGenerator& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      StaticSVDBasisGenerator&
      operator = (
         const StaticSVDBasisGenerator& rhs);

      /**
       * @brief Pointer to the time step control object.
       */
      boost::shared_ptr<StaticSVDTimeStepper> d_svdts;
};

}

#endif
