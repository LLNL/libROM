#include "incremental_svd_rom.h"

namespace CAROM {

incremental_svd_rom::incremental_svd_rom(
   int dim,
   double epsilon,
   bool skip_redundant,
   int rom_size,
   double tolerance,
   double max_time_between_snapshots,
   bool fast_update) :
   d_isvdts(new incremental_svd_time_stepper(dim,
                                             epsilon,
                                             skip_redundant,
                                             rom_size,
                                             tolerance,
                                             max_time_between_snapshots,
                                             fast_update))
{
}

incremental_svd_rom::~incremental_svd_rom()
{
}

}
