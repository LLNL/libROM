#ifndef included_utilities_h
#define included_utilities_h

#include "CAROM_config.h"
#include <iostream>
#include <sstream>
#include <string>

namespace CAROM {

/*!
 * A null use of a variable, use to avoid GNU compiler
 * warnings about unused variables.
 */
#define NULL_USE(variable)                               \
   do {                                                  \
      if (0) { char* temp = (char *)&variable; temp++; } \
   } while (0)

#ifdef DEBUG_CHECK_ASSERTIONS

#define CAROM_ASSERT(EXP)                                       \
   do {                                                         \
      if (!(EXP)) {                                             \
         std::ostringstream os;                                 \
         os << "Failed assertion: " << # EXP << std::ends;      \
         CAROM::Utilities::abort(os.str(), __FILE__, __LINE__); \
      }                                                         \
   } while (0)
#else

/*
 * No assertion checking
 */
#define CAROM_ASSERT(EXP)

#endif

/* Holder of basic, static, utility routines. */
struct Utilities
{
   // Cleanly end the program when something horrible happend and print a
   // message about what took place.  Take into account whether MPI is or isn't
   // running to decide how to die.
   static void
   abort(
      const std::string& message,
      const std::string& filename,
      int line);

   // Convert a processor ID to a string.
   static std::string
   processorToString(
      int processorID);
};

}

#endif
