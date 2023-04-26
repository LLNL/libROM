/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <iostream>
#include <sstream>

using namespace std;

void compareFiles(ifstream &baselineFile, ifstream &targetFile, double errorBound) {
    string baselineLine, targetLine;
    double baselineNum, targetNum;
    int fileLine = 1;
    bool outOfRange;

    while(!baselineFile.eof()) {
        outOfRange = false;
        getline(baselineFile, baselineLine);
        getline(targetFile, targetLine);
        if (baselineLine == "" || targetLine == "") {
            assert(baselineLine == targetLine || !(cerr << "The files are not the same length." << endl));
            break;
        }

        auto posOfData = baselineLine.find_last_of(' ');
        auto stripped = baselineLine.substr(posOfData != string::npos ? posOfData : 0);
        baselineNum = stod(stripped);

        posOfData = targetLine.find_last_of(' ');
        stripped = targetLine.substr(posOfData != string::npos ? posOfData : 0);
        targetNum = stod(stripped);
        double diff = baselineNum - targetNum;
        double error;
        if (baselineNum == 0.0) {
            error = abs(baselineNum - targetNum);
        }
        else {
            error = abs(baselineNum - targetNum) / baselineNum;
        }

        if (error > errorBound) {
            cerr << "baseline = " << baselineNum << ", diff = " << diff << endl;
            cerr << "error = " << error << endl;
            cerr << "Error bound: " << errorBound << " was surpassed on line: " << fileLine << endl;
            abort();
        }
        fileLine++;
    }
    assert(targetFile.eof() || !(cerr << "The files are not the same length." << endl));
}

int main(int argc, char *argv[]) {
    ifstream baselineFile, targetFile;
    baselineFile.open((string) argv[1]);
    targetFile.open((string) argv[2]);
    double errorBound = stod(argv[3]);
    compareFiles(baselineFile, targetFile, errorBound);

    return 0;
}