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

#include "Utilities.h"

namespace CAROM {

void
RowInfoMax(RowInfo* a, RowInfo* b, int* len, MPI_Datatype* type)
{
    for (int i = 0; i < *len; ++i) {
        if (a[i].row_val > b[i].row_val || (a[i].row_val == b[i].row_val
                                            && a[i].proc < b[i].proc)) {
            b[i].row_val = a[i].row_val;
            b[i].row = a[i].row;
            b[i].proc = a[i].proc;
        }
    }
}

}
