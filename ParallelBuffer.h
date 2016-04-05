/******************************************************************************
 *
 * Copyright (c) 2013-2016, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 * Written by William Arrighi wjarrighi@llnl.gov
 * CODE-686965
 * All rights reserved.
 *
 * This file is part of libROM.
 * For details, see https://computation.llnl.gov/librom
 * Please also read README_BSD_NOTICE.
 *
 * Redistribution and use in source and binary forms, with or without
 * modifications, are permitted provided that the following conditions are met:
 *
 *    o Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the disclaimer below.
 *    o Redistribution in binary form must reproduce the above copyright
 *      notice, this list of conditions and the disclaimer (as noted below) in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *    o Neither the name of the LLNS/LLNL nor the names of its contributors may
 *      be used to endorse or promote products derived from this software
 *      without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
 * LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OR SUCH DAMAGE.
 *
 *****************************************************************************/

// Description: A simple I/O stream class that intercepts output from an
//              ostream and redirects the output as necessary for parallel
//              I/O.

#ifndef included_ParallelBuffer_h
#define included_ParallelBuffer_h

#include "CAROM_config.h"

#include <iostream>
#include <string>

namespace CAROM {

/**
 * Class ParallelBuffer is a simple I/O stream utility that intercepts output
 * from an ostream and redirects the output as necessary for parallel I/O.
 * This class defines a stream buffer class for an ostream class.
 */
class ParallelBuffer : public std::streambuf
{
   public:
      /**
       * @brief Default constructor.
       *
       * The object will require further initialization to set up the I/O
       * streams and prefix string.
       */
      ParallelBuffer();

      /**
       * @brief Destructor.
       *
       * Simply deallocates any internal data buffers.  It does not modify the
       * output streams.
       */
      ~ParallelBuffer();

      /**
       * @brief Write a text string to the output stream.
       *
       * Note that the string is not actually written until an end-of-line is
       * detected.
       *
       * @param[in] text The string to be written.
       */
      void
      outputString(
         const std::string& text)
      {
         outputString(text, static_cast<int>(text.length()));
      }

      /**
       * @brief Write a text string of the specified length to the output file.
       *
       * Note that the string is not actually written until an end-of-line is
       * detected.
       *
       * @param[in] text The string to be written.
       * @param[in] length The length of the string.
       */
      void
      outputString(
         const std::string& text,
         int length);

      /**
       * @brief Synchronize the parallel buffer (called from streambuf).
       *
       * @return 0
       */
      int
      sync();

#if !defined(__INTEL_COMPILER) && (defined(__GNUG__))
      /**
       * @brief Write the specified number of characters into the output stream
       * (called from streambuf).
       *
       * @param[in] text The string to write.
       * @param[in] n The number of characters of the string to write.
       *
       * @return n
       */
      std::streamsize
      xsputn(
         const std::string& text,
         std::streamsize n);
#endif

      /**
       * @brief Write an overflow character into the parallel buffer (called
       * from streambuf).
       *
       * @param[in] ch The character to write.
       *
       * @return 0
       */
      int
      overflow(
         int ch);

#ifdef _MSC_VER
      /**
       * @brief Read an overflow character from the parallel buffer (called
       * from streambuf).
       *
       * This is not implemented.  It is needed by the MSVC++ stream
       * implementation.
       *
       * @return EOF
       */
      int
      underflow()
      {
         return EOF;
      }
#endif

   private:
      /**
       * @brief Copy data from the text string into the internal output buffer.
       *
       * If the internal buffer is not large enough to hold all of the string
       * data, then allocate a new internal buffer.
       *
       * @param[in] text String to copy to the internal output buffer.
       * @param[in] length Length of string to copy.
       */
      void
      copyToBuffer(
         const std::string& text,
         int length);

      /**
       * @brief Output internal buffer data to streams/
       */
      void
      outputBuffer();

      /**
       * @brief Prefix to prepend to output strings.
       */
      std::string d_prefix;

      /**
       * @brief Output stream for buffer.
       */
      std::ostream* d_ostream;

      /**
       * @brief Internal buffer to store accumulated string.
       */
      char* d_buffer;

      /**
       * @brief Size of the internal output buffer.
       */
      int d_buffer_size;

      /**
       * @brief Number of charcters in the output buffer.
       */
      int d_buffer_ptr;

      static const int DEFAULT_BUFFER_SIZE;
};

}

#endif
