/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: An abstract class defining the interface to the generic SVD
 *              algorithm.
 *
 *****************************************************************************/

#include "SVD.h"

namespace CAROM {

SVD::SVD(
   int dim,
   int samples_per_time_interval,
   bool debug_rom) :
   d_dim(dim),
   d_num_samples(0),
   d_samples_per_time_interval(samples_per_time_interval),
   d_basis(0),
   d_time_interval_start_times(0),
   d_debug_rom(debug_rom)
{
   CAROM_ASSERT(dim > 0);
   CAROM_ASSERT(samples_per_time_interval > 0);
}

SVD::~SVD()
{
   if (d_basis) {
      delete d_basis;
   }
}

}
