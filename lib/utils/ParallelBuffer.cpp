/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: A simple I/O stream class that intercepts output from an
//              ostream and redirects the output as necessary for parallel
//              I/O.

#include "ParallelBuffer.h"
#include "Utilities.h"

#include "mpi.h"

#include <string>
#include <cstring>
#include <cstdio>

namespace CAROM {

const int ParallelBuffer::DEFAULT_BUFFER_SIZE = 128;

/*
 *************************************************************************
 *
 * Construct a parallel buffer object.  The object will require further
 * initialization to set up I/O streams and the prefix string.
 *
 *************************************************************************
 */
ParallelBuffer::ParallelBuffer()
{
    int mpi_init;
    MPI_Initialized(&mpi_init);
    int rank;
    if (mpi_init) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
    else {
        rank = 0;
    }
    d_prefix = "P=" + Utilities::processorToString(rank) + ":";
    d_ostream = &std::cerr;
    d_buffer = 0;
    d_buffer_size = 0;
    d_buffer_ptr = 0;
}

/*
 *************************************************************************
 *
 * The destructor deallocates internal data buffer.  It does not modify
 * the output streams.
 *
 *************************************************************************
 */
ParallelBuffer::~ParallelBuffer()
{
    if (d_buffer) {
        delete[] d_buffer;
    }
}

/*
 *************************************************************************
 *
 * Write a text string of the specified length to the output stream.
 * Note that the string data is accumulated into the internal output
 * buffer until an end-of-line is detected.
 *
 *************************************************************************
 */
void
ParallelBuffer::outputString(
    const std::string& text,
    int length)
{
    if (length > 0) {

        /*
         * If we need to allocate the internal buffer, then do so.
         */
        if (!d_buffer) {
            d_buffer = new char[DEFAULT_BUFFER_SIZE];
            d_buffer_size = DEFAULT_BUFFER_SIZE;
            d_buffer_ptr = 0;
        }

        /*
         * If the buffer pointer is zero, then prepend the prefix.
         */
        if (d_buffer_ptr == 0) {
            copyToBuffer(d_prefix, static_cast<int>(d_prefix.length()));
        }

        /*
         * Search for an end-of-line in the string.
         */
//      const int eol_ptr = static_cast<int>(text.find('\n'));
        int eol_ptr = 0;
        for ( ; (eol_ptr < length) && (text[eol_ptr] != '\n'); eol_ptr++) {}

        /*
         * If no end-of-line is found, copy the entire text string but do not
         * output.  Otherwise copy the text string through the end-of-line,
         * output, and recurse with the remainder of the text string if there are
         * more characters in it.
         */
        if (eol_ptr == length) {
            copyToBuffer(text, length);
        } else {
            const int ncopy = eol_ptr + 1;
            copyToBuffer(text, ncopy);
            outputBuffer();
            if (ncopy < length) {
                outputString(text.substr(ncopy), length - ncopy);
            }
        }
    }
}

/*
 *************************************************************************
 *
 * Synchronize the parallel buffer and write string data.  This routine
 * is called from streambuf.
 *
 *************************************************************************
 */
int
ParallelBuffer::sync()
{
    const int n = static_cast<int>(pptr() - pbase());
    if (n > 0) {
        outputString(pbase(), n);
    }
    return 0;
}

/*
 *************************************************************************
 *
 * Write the specified number of characters into the output stream.
 * This routine is called from streambuf.  If this routine is not
 * provided, then overflow() is called instead for each character.
 *
 * Note that this routine is not required; it only
 * offers some efficiency over overflow().
 *
 *************************************************************************
 */
#if !defined(__INTEL_COMPILER) && (defined(__GNUG__))
std::streamsize
ParallelBuffer::xsputn(
    const std::string& text,
    std::streamsize n)
{
    sync();
    if (n > 0) outputString(text, static_cast<int>(n));
    return n;
}
#endif

/*
 *************************************************************************
 *
 * Write a single character into the parallel buffer.  This routine is
 * called from streambuf.
 *
 *************************************************************************
 */
int
ParallelBuffer::overflow(
    int ch)
{
    const int n = static_cast<int>(pptr() - pbase());
    if (n && sync()) {
        return EOF;
    }
    if (ch != EOF) {
        char character[2];
        character[0] = (char)ch;
        character[1] = 0;
        outputString(character, 1);
    }
    pbump(-n);
    return 0;
}

/*
 *************************************************************************
 *
 * Copy data from the text string into the internal output buffer.
 * If the internal buffer is not large enough to hold all of the string
 * data, then allocate a new internal buffer.
 *
 *************************************************************************
 */
void
ParallelBuffer::copyToBuffer(
    const std::string& text,
    int length)
{
    /*
     * First check whether we need to increase the size of the buffer
     */
    if (d_buffer_ptr + length > d_buffer_size) {
        int new_size;
        if (d_buffer_ptr + length > 2 * d_buffer_size) {
            new_size = d_buffer_ptr + length;
        }
        else {
            new_size = 2 * d_buffer_size;
        }
        char* new_buffer = new char[new_size];

        if (d_buffer_ptr > 0) {
            (void)strncpy(new_buffer, d_buffer, d_buffer_ptr);
        }
        delete[] d_buffer;

        d_buffer = new_buffer;
        d_buffer_size = new_size;
    }
    CAROM_ASSERT(d_buffer_ptr + length <= d_buffer_size);

    /*
     * Copy data from the input into the internal buffer and increment pointer
     */
    strncpy(d_buffer + d_buffer_ptr, text.c_str(), length);
    d_buffer_ptr += length;
}

/*
 *************************************************************************
 *
 * Output buffered stream data to the active output streams and reset
 * the buffer pointer to its empty state.
 *
 *************************************************************************
 */
void
ParallelBuffer::outputBuffer()
{
    if (d_buffer_ptr > 0) {
        if (d_ostream) {
            d_ostream->write(d_buffer, d_buffer_ptr);
            d_ostream->flush();
        }
        d_buffer_ptr = 0;
    }
}

}
