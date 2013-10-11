#include "incremental_svd_rom.h"

namespace CAROM {

incremental_svd_rom::incremental_svd_rom(
   int dim,
   double epsilon,
   bool skip_redundant,
   int rom_size,
   int max_time_steps_between_snapshots) :
   d_isvdts(new incremental_svd_time_stepper(dim,
                                             epsilon,
                                             skip_redundant,
                                             rom_size,
                                             max_time_steps_between_snapshots))
{
}

incremental_svd_rom::~incremental_svd_rom()
{
}

}
