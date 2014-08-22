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

// A class which implements a Reduced Order Model based on the static svd
// algorithm.
class StaticSVDBasisGenerator : public SVDBasisGenerator
{
   public:
      // Constructor.
      StaticSVDBasisGenerator(
         int dim,
         int increments_per_time_interval,
         const std::string& basis_file_name,
         Database::formats file_format = Database::HDF5);

      // Destructor.
      virtual
      ~StaticSVDBasisGenerator();

      // Returns true if it is time for the next svd snapshot.
      virtual
      bool
      isNextSnapshot(
         double time);

      // Add a snapshot to the static svd.
      virtual
      void
      takeSnapshot(
         double* u_in,
         double time);

      // Computes next time an svd snapshot is needed.
      virtual
      double
      computeNextSnapshotTime(
         double* u_in,
         double* rhs_in,
         double time);

      // Returns the basis vectors for the current time interval as a Matrix.
      virtual
      const Matrix*
      getBasis();

      // Returns the number of time intervals on which different sets of basis
      // vectors are defined.
      virtual
      int
      getNumBasisTimeIntervals() const;

      // Returns the start time for the requested time interval.
      virtual
      double
      getBasisIntervalStartTime(
         int which_interval) const;

   private:
      // Unimplemented default constructor.
      StaticSVDBasisGenerator();

      // Unimplemented copy constructor.
      StaticSVDBasisGenerator(
         const StaticSVDBasisGenerator& other);

      // Unimplemented assignment operator.
      StaticSVDBasisGenerator&
      operator = (
         const StaticSVDBasisGenerator& rhs);

      boost::shared_ptr<StaticSVDTimeStepper> d_svdts;
};

}

#endif
