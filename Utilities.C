#include "Utilities.h"
#include "ParallelBuffer.h"

#include "mpi.h"

#include <iomanip>
#include <stdlib.h>

namespace CAROM {

void
Utilities::abort(
   const std::string& message,
   const std::string& filename,
   int line)
{
   static ParallelBuffer perr_buffer;
   std::ostringstream os;
   os << "Program abort called in file ``" << filename
      << "'' at line " << line << std::endl;
   perr_buffer.outputString(os.str());
   os.str("");
   os << "ERROR MESSAGE: " << std::endl << message.c_str() << std::endl;
   perr_buffer.outputString(os.str());

   int size;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   if (size > 1) {
      MPI_Abort(MPI_COMM_WORLD, -1);
   }
   else {
      exit(-1);
   }
}

std::string
Utilities::processorToString(
   int processorID)
{
   std::ostringstream os;
   if (processorID < 0) {
      os << '-' << std::setw(6) << std::setfill('0') << -processorID;
   } else {
      os << std::setw(7) << std::setfill('0') << processorID;
   }
   os << std::flush;

   return os.str();
}

}
