/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/
#include <string>
#include <cmath>
#include <cassert>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {
    double offlineTime = stod(argv[1]);
    double onlineTime = stod(argv[2]);
    double speedupTol = stod(argv[3]);
    if (offlineTime / onlineTime < speedupTol) {
        cerr << "offlineTime / onlineTime = " << offlineTime / onlineTime << endl;
        cerr << "speedupTol: " << speedupTol << " was not surpassed." << endl;
        abort();
    }

    return 0;
}