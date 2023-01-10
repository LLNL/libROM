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
#include <iostream>
#include "mfem.hpp"

using namespace std;
using namespace mfem;

int getDimensions(string &filePath) {
    string line;
    int count = 0;

    ifstream file(filePath);
    while (getline(file, line)) {
        count++;
    }
    return count;
}

void compareSolutions(string &baselineFile, string &targetFile, double errorBound, int numProcessors) {
    int* baselineDim = new int[numProcessors];
    istream** baselineFiles = new istream*[numProcessors];
    int* targetDim = new int[numProcessors];
    istream** targetFiles = new istream*[numProcessors];
    std::filebuf* baselinefb = new filebuf[numProcessors];
    std::filebuf* targetfb = new filebuf[numProcessors];

    for (int i = 0; i < numProcessors; i++) {
        if (i > 0) {
            baselineFile.back() = '0' + i;
            targetFile.back() = '0' + i;
        }
        if (baselinefb[i].open(baselineFile, ios::in)) {
            baselineFiles[i] = new istream(&baselinefb[i]);
            baselineDim[i] = getDimensions(baselineFile);
        }
        else {
            cerr << "Something went wrong with opening the following file. Most likely it doesn't exist: " << baselineFile << endl;
            abort();
        }
        cout << "Solution Comparator is Opening file: " << targetFile << endl;
        if (targetfb[i].open(targetFile, ios::in)) {
            targetFiles[i] = new istream(&targetfb[i]);
            targetDim[i] = getDimensions(targetFile);
        }
        else {
            cerr << "Something went wrong with opening the following file. Most likely it doesn't exist: " << targetFile << endl;
            abort();
        }
    }
    Vector baseline = Vector();
    Vector target = Vector();
    baseline.Load(baselineFiles, numProcessors, baselineDim);
    target.Load(targetFiles, numProcessors, targetDim);

    if (baseline.Size() != target.Size()) {
        cerr << "The solution vectors are different dimensions." << endl;
        cerr << "Baseline dim: " << baseline.Size() << ", Target dim: " << target.Size() << endl;
        abort();
    }

    Vector diff = Vector(baseline.Size());
    try {
        subtract(baseline, target, diff);
    }
    catch (const exception& e) {
        cerr << "Something went wrong when calculating the difference \
between the solution vectors." << endl;
        abort();
    }
    double baselineNormL2 = baseline.Norml2();
    double diffNormL2 = diff.Norml2();
    // Check for NaN
    if(std::isnan(baselineNormL2)){
        std::cerr << "baselineNormL2 is NaN" << std::endl;
        if(std::isnan(diffNormL2)){
            std::cerr << "diffNormL2 is NaN" << std::endl;
        }
        abort();
    }
    if(std::isnan(diffNormL2)){
        std::cerr << "diffNormL2 is NaN" << std::endl;
        abort();
    }
    double error;
    if (baselineNormL2 == 0.0) {
        error = diffNormL2;
    }
    else {
        error = diffNormL2 / baselineNormL2;
    }


    if(std::isnan(baselineNormL2)){
        std::cout << "baselineNormL2 is NaN" << std::endl;
        abort();
    }

    if(std::isnan(diffNormL2)){
        std::cout << "diffNormL2 is NaN" << std::endl;
        abort();
    }


    // Test whether l2 norm is smaller than error bound
    if (error > errorBound) {
        cerr << "baselineNormL2 = " << baselineNormL2 << ", diffNormL2 = " << diffNormL2 << endl;
        cerr << "error = " << error << endl;
        cerr << "Error bound: " << errorBound << " was surpassed for the l2 norm of the difference of the solutions." << endl;
        abort();
    }
}

int main(int argc, char *argv[]) {
    string baselinePath((string) argv[1]);
    string targetPath((string) argv[2]);
    double errorBound = stod(argv[3]);
    int numProcessors = stoi(argv[4]);

    compareSolutions(baselinePath, targetPath, errorBound, numProcessors);
    return 0;
}