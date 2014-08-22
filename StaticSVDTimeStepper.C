/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The class that determines the next time at which a sample
 *              should be taken for basis generation using the static SVD
 *              approach.
 *
 *****************************************************************************/

#include "StaticSVDTimeStepper.h"

namespace CAROM {

StaticSVDTimeStepper::StaticSVDTimeStepper(
   int dim,
   int increments_per_time_interval) :
   d_svd(new StaticSVD(dim, increments_per_time_interval))
{
}

StaticSVDTimeStepper::~StaticSVDTimeStepper()
{
}

}
