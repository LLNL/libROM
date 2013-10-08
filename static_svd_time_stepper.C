#include "static_svd_time_stepper.h"

namespace CAROM {

static_svd_time_stepper::static_svd_time_stepper(
   int dim) :
   d_svd(new static_svd(dim))
{
}

static_svd_time_stepper::~static_svd_time_stepper()
{
}

}
