/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

#ifdef CAROM_HAS_GTEST

#include<gtest/gtest.h>
#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include "utils/HDFDatabase.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring> // for memcpy
#include <random>
#include "mpi.h"
#include "utils/mpi_utils.h"

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
    SUCCEED();
}

TEST(HDFDatabase, Test_file_access_property_list)
{
    hid_t access_plist;
    hid_t file_id;
    
    access_plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpi(access_plist, MPI_COMM_WORLD, MPI_INFO_NULL);
    
    // H5Fopen must be called collectively
    file_id = H5Fopen("test.h5", H5F_ACC_RDWR, access_plist);
    
    // H5Fclose must be called collectively
    H5Fclose(file_id);
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}

#else // #ifndef CAROM_HAS_GTEST

int main()
{
    std::cout << "libROM was compiled without Google Test support, so unit "
              << "tests have been disabled. To enable unit tests, compile "
              << "libROM with Google Test support." << std::endl;
}

#endif // #endif CAROM_HAS_GTEST