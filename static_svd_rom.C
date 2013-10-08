#include "static_svd_rom.h"

namespace CAROM {

static_svd_rom::static_svd_rom(
   int dim) :
   d_svdts(new static_svd_time_stepper(dim))
{
}

static_svd_rom::~static_svd_rom()
{
}

}
