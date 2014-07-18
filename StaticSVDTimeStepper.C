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
