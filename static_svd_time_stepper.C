#include "static_svd_time_stepper.h"

namespace CAROM {

static_svd_time_stepper::static_svd_time_stepper(
   int dim,
   int increments_per_time_interval) :
   d_svd(new static_svd(dim, increments_per_time_interval))
{
}

static_svd_time_stepper::~static_svd_time_stepper()
{
}

}
