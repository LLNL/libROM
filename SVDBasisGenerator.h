/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The abstract wrapper class for an abstract SVD algorithm and
 *              sampler.  This class provides interfaces to each so that
 *              an application only needs to instantiate one concrete
 *              implementation of this class to control all aspects of basis
 *              vector generation.
 *
 *****************************************************************************/

#ifndef included_SVDBasisGenerator_h
#define included_SVDBasisGenerator_h

#include "BasisWriter.h"

#include <string.h>

namespace CAROM {

class BasisWriter;
class Matrix;

/**
 * Class SVDBasisGenerator is an abstract base class defining the interface for
 * the generation of basis vectors via the svd method.  This class wraps the
 * abstract SVD algorithm and sampler and provides interfaces to each so
 * that an application only needs to instantiate one concrete implementation of
 * this class to control all aspects of basis vector generation.
 */
class SVDBasisGenerator
{
   public:
      /**
       * @brief Constructor.
       *
       * @param[in] basis_file_name The base part of the name of the file
       *                            containing the basis vectors.  Each process
       *                            will append its process ID to this base
       *                            name.
       * @param[in] file_format The format of the file containing the basis
       *                        vectors.
       */
      SVDBasisGenerator(
         const std::string& basis_file_name,
         Database::formats file_format = Database::HDF5);

      /**
       * @brief Destructor.
       */
      virtual
      ~SVDBasisGenerator();

      /**
       * @brief Returns true if it is time for the next svd sample.
       *
       * @pre time >= 0.0
       *
       * @param[in] time Time of interest.
       *
       * @return True if it is time for the next sample to be taken.
       */
      virtual
      bool
      isNextSample(
         double time) = 0;

      /**
       * @brief Sample the new state, u_in, at the given time.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] time The simulation time for the state.
       */
      virtual
      void
      takeSample(
         double* u_in,
         double time) = 0;

      /**
       * @brief Signal that the final sample has been taken.
       */
      void
      endSamples()
      {
         if (d_basis_writer) {
            d_basis_writer->writeBasis();
         }
      }

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
      virtual
      double
      computeNextSampleTime(
         double* u_in,
         double* rhs_in,
         double time) = 0;

      /**
       * @brief Returns the basis vectors for the current time interval as a
       * Matrix.
       *
       * @return The basis vectors for the current time interval.
       */
      virtual
      const Matrix*
      getBasis() = 0;

      /**
       * @brief Returns the number of time intervals on which different sets of
       * basis vectors are defined.
       *
       * @return The number of time intervals on which there are basis vectors.
       */
      virtual
      int
      getNumBasisTimeIntervals() const = 0;

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
         int which_interval) const = 0;

   protected:
      /**
       * @brief Writer of basis vectors.
       */
      BasisWriter* d_basis_writer;

   private:
      /**
       * @brief Unimplemented copy constructor.
       */
      SVDBasisGenerator(
         const SVDBasisGenerator& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      SVDBasisGenerator&
      operator = (
         const SVDBasisGenerator& rhs);
};

}

#endif
