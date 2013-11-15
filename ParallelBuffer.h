#ifndef included_ParallelBuffer_h
#define included_ParallelBuffer_h

#include "CAROM_config.h"

#include <iostream>
#include <string>

namespace CAROM {

/**
 * Class ParallelBuffer is a simple I/O stream utility that
 * intercepts output from an ostream and redirects the output as necessary
 * for parallel I/O.  This class defines a stream buffer class for an
 * ostream class.
 */
class ParallelBuffer:public std::streambuf
{
   public:
      /**
       * Create a parallel buffer class.  The object will require further
       * initialization to set up the I/O streams and prefix string.
       */
      ParallelBuffer();

      /**
       * The destructor simply deallocates any internal data
       * buffers.  It does not modify the output streams.
       */
      ~ParallelBuffer();

      /**
       * Write a text string to the output stream.  Note that the string is
       * not actually written until an end-of-line is detected.
       */
      void
      outputString(
         const std::string& text)
      {
         outputString(text, static_cast<int>(text.length()));
      }

      /**
       * Write a text string of the specified length to the output file.  Note
       * that the string is not actually written until an end-of-line is
       * detected.
       */
      void
      outputString(
         const std::string& text,
         const int length);

      /**
       * Synchronize the parallel buffer (called from streambuf).
       */
      int
      sync();

#if !defined(__INTEL_COMPILER) && (defined(__GNUG__))
      /**
       * Write the specified number of characters into the output stream
       * (called from streambuf).
       */
      std::streamsize
      xsputn(
         const std::string& text,
         std::streamsize n);
#endif

      /**
       * Write an overflow character into the parallel buffer (called from
       * streambuf).
       */
      int
      overflow(
         int ch);

#ifdef _MSC_VER
      /**
       * Read an overflow character from the parallel buffer (called from
       * streambuf).  This is not implemented.  It is needed by the
       * MSVC++ stream implementation.
       */
      int
      underflow()
      {
         return EOF;
      }
#endif

   private:
      void
      copyToBuffer(
         const std::string& text,
         const int length);
      void
      outputBuffer();           // output internal buffer data to streams

      std::string d_prefix;     // string prefix to prepend output strings
      std::ostream* d_ostream;  // output stream for buffer
      char* d_buffer;           // internal buffer to store accumulated string
      int d_buffer_size;        // size of the internal output buffer
      int d_buffer_ptr;         // number of charcters in the output buffer

      static const int DEFAULT_BUFFER_SIZE;
};

}

#endif
