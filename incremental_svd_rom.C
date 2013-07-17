#include "incremental_svd_rom.h"

namespace CAROM {

incremental_svd_rom::incremental_svd_rom(
   int* argc,
   char*** argv,
   int dim,
   double epsilon,
   bool skip_redundant,
   int max_time_steps_between_snapshots) :
   d_isvdts(new incremental_svd_time_stepper(argc,
                                             argv,
                                             dim,
                                             epsilon,
                                             skip_redundant,
                                             max_time_steps_between_snapshots))
{
}

incremental_svd_rom::~incremental_svd_rom()
{
}

}
