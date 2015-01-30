/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2015 Lawrence Livermore National Security, LLC
 * Description: An abstract class defining the interface to the generic SVD
 *              algorithm.
 *
 *****************************************************************************/

#ifndef included_SVD_h
#define included_SVD_h

#include "Matrix.h"
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
       * @pre dim > 0
       * @pre sample_per_time_interval > 0
       *
       * @param[in] dim The dimension of the system distributed to this
       *                processor.
       * @param[in] samples_per_time_interval The maximium number of samples
       *                                      collected in a time interval.
       * @param[in] debug_algorithm If true results of the algorithm will be
       *                            printed to facilitate debugging.
       */
      SVD(
         int dim,
         int samples_per_time_interval,
         bool debug_algorithm = false);

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
       */
      virtual
      void
      takeSample(
         const double* u_in,
         double time) = 0;

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
      getBasis() = 0;

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
         CAROM_ASSERT(0 <= which_interval);
         CAROM_ASSERT(which_interval < getNumBasisTimeIntervals());
         return d_time_interval_start_times[which_interval];
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
       * @brief The maximum number of samples to be collected for a time
       * interval.
       */
      int d_samples_per_time_interval;

      /**
       * @brief The globalized basis vectors for the current time interval.
       *
       * The basis vectors are large and each process owns all of the basis
       * vectors.
       */
      Matrix* d_basis;

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
