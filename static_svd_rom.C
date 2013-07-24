#include "static_svd_rom.h"

namespace CAROM {

static_svd_rom::static_svd_rom(
   int* argc,
   char*** argv,
   int dim) :
   d_svdts(new static_svd_time_stepper(argc, argv, dim))
{
}

static_svd_rom::~static_svd_rom()
{
}

}
