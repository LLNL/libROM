/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
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
    double relError= stod(argv[1]);
    double errorBound = stod(argv[2]);
    if (abs(relError) > errorBound) {
        cerr << "errorBound = " << errorBound << endl;
        cerr << "abs(relError) = " << abs(relError) << endl;
        cerr << "Error bound: " << errorBound << " was surpassed." << endl;
        abort();
    }

    return 0;
}