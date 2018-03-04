/******************************************************************************
 *
 * Copyright (c) 2013-2017, Lawrence Livermore National Security, LLC.
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

// Description: This source file is a test runner that uses the Google Test
// Framework to run unit tests.

#include<gtest/gtest.h>
#include <mpi.h>
#include "../Matrix.h"

// Test fixture for testing
TEST(MatrixSerialQRCPTest, SecondDifferenceMatrix) {
  // Allocate space for second_difference matrix
  const int size = 4;
  const bool is_distributed = false;
  CAROM::Matrix second_difference(size, size, is_distributed);

  /* A "second difference" matrix is a tridiagonal matrix that looks like:

     [ 2 -1         ]
     [-1  2 -1      ]
     [   .. .. ..   ]
     [      -1  2 -1]
     [         -1  2]

     which looks like the stencil from second-order central
     differencing (or the stiffness matrix for finite elements on a
     subspace of H1 with linear basis functions) for Poisson's
     equation, modulo boundary conditions (that would be applied to
     the first and last rows).
   */

  // The code below won't work if the size is less than 2.
  assert(size >= 2);

  // Zero out the entire matrix
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      second_difference.item(i, j) = 0;
    }
  }

  // Build first row of second difference matrix
  second_difference.item(0, 0) =  2;
  second_difference.item(0, 1) = -1;

  // Build second through penultimate rows of
  // second difference matrix
  for (int row = 1; row < size - 1; row++) {
    second_difference.item(row, row - 1) = -1;
    second_difference.item(row, row    ) =  2;
    second_difference.item(row, row + 1) = -1;
  }

  // Build last row of second difference matrix
  second_difference.item(size - 1, size - 2) = -1;
  second_difference.item(size - 1, size - 1) =  2;

  // Allocate memory for QRCP arguments
  int* row_pivot                 =       new int[size];
  int* row_pivot_owner           =       new int[size];
  const int row_pivots_requested =               size ;

  // Compute pivots
  second_difference.qrcp_pivots_transpose(row_pivot,
					  row_pivot_owner,
					  row_pivots_requested);

  // All row pivot owners should be on the current rank;
  // all row pivots shouold be less than the size of the matrix
  int is_mpi_initialized, is_mpi_finalized;
  CAROM_ASSERT(MPI_Initialized(&is_mpi_initialized) == MPI_SUCCESS);
  CAROM_ASSERT(MPI_Finalized(&is_mpi_finalized) == MPI_SUCCESS);
  int my_rank = 0;
  if(is_mpi_initialized && !is_mpi_finalized) {
    const MPI_Comm my_comm = MPI_COMM_WORLD;
    CAROM_ASSERT(MPI_Comm_rank(my_comm, &my_rank) == MPI_SUCCESS);
  }

  for (int i = 0; i < row_pivots_requested; i++) {
    EXPECT_EQ(row_pivot_owner[i], my_rank);
    EXPECT_TRUE(row_pivot[i] < size);
  }

  // Row pivots should be known values
  EXPECT_EQ(row_pivot[0], 1);
  EXPECT_EQ(row_pivot[1], 3);
  EXPECT_EQ(row_pivot[2], 0);
  EXPECT_EQ(row_pivot[3], 2);

  // Free allocated arrays
  delete [] row_pivot;
  delete [] row_pivot_owner;
}

// Simple smoke test to make sure runner is properly linked
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
  SUCCEED();
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
