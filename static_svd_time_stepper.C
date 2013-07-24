#include "static_svd_time_stepper.h"

namespace CAROM {

static_svd_time_stepper::static_svd_time_stepper(
   int* argc,
   char*** argv,
   int dim) :
   d_svd(new static_svd(argc, argv, dim))
{
}

static_svd_time_stepper::~static_svd_time_stepper()
{
}

}
