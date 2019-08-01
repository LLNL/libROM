/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The abstract wrapper class for an abstract SVD algorithm and
//              sampler.  This class provides interfaces to each so that an
//              application only needs to instantiate one concrete
//              implementation of this class to control all aspects of basis
//              vector generation.

#ifndef included_SVDBasisGenerator_h
#define included_SVDBasisGenerator_h

#include "BasisWriter.h"
#include "BasisReader.h"
#include "SVDSampler.h"

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
      bool
      isNextSample(
         double time)
      {
         CAROM_ASSERT(time >= 0.0);
         return d_svdsampler->isNextSample(time);
      }

      /**
       * @brief Returns true if it needs to update right basis vectors.
       */
      bool
      updateRightSV() {return d_svdsampler->isUpdateRightSV(); };

      /**
       * @brief Sample the new state, u_in, at the given time.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] time The simulation time for the state.
       * @param[in] dt The current simulation dt.
       *
       * @return True if the sampling was successful.
       */
      bool
      takeSample(
         double* u_in,
         double time,
         double dt,
         bool add_without_increase = false)
      {
         CAROM_ASSERT(u_in != 0);
         CAROM_ASSERT(time >= 0);

         // Check that u_in is not non-zero.
         Vector u_vec(u_in, d_svdsampler->getDim(), true);
         if (u_vec.norm() == 0.0) {
            return false;
         }


         if (getNumBasisTimeIntervals() > 0 &&
             d_svdsampler->isNewTimeInterval()) {
            d_svdsampler->resetDt(dt);
            // YC: commenting this out unless someone think this is necessary
            //     we call writeBasis in "endSamples()" function below.
            //     I think that one call is enough.
            //     I will remove completely if no one has other opinion.
            //if (d_basis_writer) {
            //   d_basis_writer->writeBasis();
            //}
         }

         return d_svdsampler->takeSample(u_in, time, add_without_increase);
      }

      /**
       * @brief Signal that the final sample has been taken.
       *
       * @param[in] kind "basis" or "snapshot" to write one or the other.
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
       * @param[in] db_format Format of the file to read.
       */
      void
      loadSamples(const std::string& base_file_name,
                  const std::string& kind = "basis",
                  Database::formats db_format = Database::HDF5)
      { 
	      CAROM_ASSERT(!base_file_name.empty());
         CAROM_ASSERT(kind == "basis" || kind == "snapshot");
         
         if (d_basis_reader) delete d_basis_reader;
         
	      d_basis_reader = new BasisReader(base_file_name, db_format);
         d_basis_reader->readBasis(base_file_name, db_format);
         double time = 0.0;
         const Matrix* mat;
         
         if (kind == "basis") {
            mat = d_basis_reader->getSpatialBasis(time);
         }
         else if (kind == "snapshot") {
            mat = d_basis_reader->getSnapshotMatrix(time);
         }
         
         int num_rows = mat->numRows();
         int num_cols = mat->numColumns();
         double* u_in = new double[num_rows*num_cols];
         for (int j = 0; j < num_cols; j++) {
            for (int i = 0; i < num_rows; i++) {
               u_in[i+j*num_rows] = mat->item(i,j);
            }
            d_svdsampler->takeSample(u_in+j*num_rows, time, false);
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
      double
      computeNextSampleTime(
         double* u_in,
         double* rhs_in,
         double time)
      {
         CAROM_ASSERT(u_in != 0);
         CAROM_ASSERT(rhs_in != 0);
         CAROM_ASSERT(time >= 0);

         return d_svdsampler->computeNextSampleTime(u_in, rhs_in, time);
      }

      /**
       * @brief Returns the basis vectors for the current time interval as a
       * Matrix.
       *
       * @return The basis vectors for the current time interval.
       */
      const Matrix*
      getSpatialBasis()
      {
         return d_svdsampler->getSpatialBasis();
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
         return d_svdsampler->getTemporalBasis();
      }

      /**
       * @brief Returns the singular values for the current time interval as a
       * Matrix.
       *
       * @return The singular values for the current time interval.
       */
      const Matrix*
      getSingularValues()
      {
         return d_svdsampler->getSingularValues();
      }

      /**
       * @brief Returns the snapshot matrix for the current time interval.
       *
       * @return The snapshot matrix for the current time interval.
       */
      const Matrix*
      getSnapshotMatrix()
      {
         return d_svdsampler->getSnapshotMatrix();
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
         return d_svdsampler->getNumBasisTimeIntervals();
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
         CAROM_ASSERT(0 <= which_interval);
         CAROM_ASSERT(which_interval < getNumBasisTimeIntervals());
         return d_svdsampler->getBasisIntervalStartTime(which_interval);
      }

   protected:
      /**
       * @brief Constructor.
       *
       * Although all member functions are implemented by delegation to either
       * d_basis_writer or d_svdsampler, this class is still abstract.  In this
       * context it is not yet known which SVDSampler to instantiate.  Hence an
       * instance of this class may not be constructed.
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
       * @brief Writer of basis vectors.
       */
      BasisWriter* d_basis_writer;

      /**
       * @brief Reader of basis vectors.
       */
      BasisReader* d_basis_reader;


      /**
       * @brief Pointer to the underlying sampling control object.
       */
#if __cplusplus >= 201103L
      std::shared_ptr<SVDSampler> d_svdsampler;
#else
      boost::shared_ptr<SVDSampler> d_svdsampler;
#endif

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
