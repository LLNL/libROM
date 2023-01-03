/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Holds common utilities used among the different sampling algorithms.

#ifndef included_Utilities_h
#define included_Utilities_h

#include "mpi.h"

namespace CAROM {

#ifndef DOXYGEN_IGNORE

/**
 * @brief Struct to hold the local maximum absolute value of a basis vector,
 *        the row it is in, and the processor that owns it. We will reduce this
 *        to find the global row containing the maximum of a basis vector.
 */
typedef struct
{
    double row_val;
    int row;
    int proc;
} RowInfo;

/**
 * @brief The function to use as an MPI_Op in the reduction to determine the row
 *        and processor owning the row of the absolute maximum of a basis vector.
 */
void RowInfoMax(RowInfo* a, RowInfo* b, int* len, MPI_Datatype* type);

#endif /* DOXYGEN_IGNORE */

}

#endif
