#include "static_svd_rom.h"

namespace CAROM {

static_svd_rom::static_svd_rom(
   int dim,
   int increments_per_time_interval) :
   d_svdts(new static_svd_time_stepper(dim, increments_per_time_interval))
{
}

static_svd_rom::~static_svd_rom()
{
}

}
