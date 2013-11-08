#include "Utilities.h"

#include <mpi.h>
#include <stdlib.h>

namespace CAROM {

void
Utilities::abort(
   const std::string& message,
   const std::string& filename,
   int line)
{
   std::cerr << "Program abort called in file ``" << filename
             << "'' at line " << line << std::endl;
   std::cerr << "ERROR MESSAGE: " << std::endl << message.c_str() << std::endl;
   std::cerr << std::flush;;

   int size;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   if (size > 1) {
      MPI_Abort(MPI_COMM_WORLD, -1);
   }
   else {
      exit(-1);
   }
}

}
