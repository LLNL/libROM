/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

#include "linalg/BasisReader.h"
#include "linalg/Matrix.h"
#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <mpi.h>

using namespace std;

void compareBasis(string &baselineFile, string &targetFile, double errorBound, int numProcessors) {

    MPI_Init(NULL, NULL);
    // Get the number of processes
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    vector<double> vecNormL2;
    vector<double> reducedVecNormL2;
    vector<double> diffVecNormL2;
    vector<double> reducedDiffVecNormL2;

    CAROM::BasisReader baselineReader(baselineFile);
    CAROM::Matrix *baselineBasis = (CAROM::Matrix*) baselineReader.getSpatialBasis(0.0);
    CAROM::Vector *baselineSV = (CAROM::Vector*) baselineReader.getSingularValues(0.0);
    CAROM::BasisReader targetReader(targetFile);
    CAROM::Matrix *targetBasis = (CAROM::Matrix*) targetReader.getSpatialBasis(0.0);
    CAROM::BasisReader diffReader(baselineFile);
    CAROM::Matrix *diffBasis = (CAROM::Matrix*) diffReader.getSpatialBasis(0.0);

    // Get basis dimensions
    int baselineNumRows = baselineBasis->numRows();
    int baselineNumColumns = baselineBasis->numColumns();
    int targetNumRows = targetBasis->numRows();
    int targetNumColumns = targetBasis->numColumns();

    // Test basis matrices have the same dimensions
    if (baselineNumRows != targetNumRows) {
        cerr << "The number of rows of the two basis matrices \
are not equal in the following files: " << baselineFile << " and " << targetFile << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (baselineNumColumns != targetNumColumns) {
        cerr << "The number of columns of the two basis matrices \
are not equal in the following file: " << baselineFile << " and " << targetFile << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (baselineSV->dim() != baselineNumColumns)
    {
        cerr << "The number of singular values does not equal the \
number of basis vectors in the following file: " << baselineFile << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    vecNormL2.resize(baselineNumColumns, 0.0);
    reducedVecNormL2.resize(baselineNumColumns, 0.0);
    diffVecNormL2.resize(baselineNumColumns, 0.0);
    reducedDiffVecNormL2.resize(baselineNumColumns, 0.0);

    try {
        *diffBasis -= (*targetBasis);
    }
    catch (const exception& e) {
        cerr << "Something went wrong when calculating the difference \
between the basis matrices in the following files: " << baselineFile << " and " << targetFile << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Compute l2-norm
    for (unsigned int i = 0; i < baselineNumColumns; i++) {
        for (unsigned int j = 0; j < baselineNumRows; j++) {
            vecNormL2[i] += pow(baselineBasis->item(j,i), 2);
            diffVecNormL2[i] += pow(diffBasis->item(j,i), 2);
        }
    }

    MPI_Reduce(vecNormL2.data(), reducedVecNormL2.data(), baselineNumColumns, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(diffVecNormL2.data(), reducedDiffVecNormL2.data(), baselineNumColumns, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double baselineNormL2 = 0.0;
        double diffNormL2 = 0.0;
        for (int i = 0; i < baselineNumColumns; i++) {
            const double weight = (*baselineSV)(i) / (*baselineSV)(0);
            baselineNormL2 += reducedVecNormL2[i];
            diffNormL2 += weight * weight * reducedDiffVecNormL2[i];
        }

        baselineNormL2 = sqrt(baselineNormL2);
        diffNormL2 = sqrt(diffNormL2);

        double error;
        if (baselineNormL2 == 0.0) {
            error = diffNormL2;
        }
        else {
            error = diffNormL2 / baselineNormL2;
        }

        // Test whether l2 norm is smaller than error bound
        if (error > errorBound) {
            cerr << "baselineNormL2 = " << baselineNormL2 << ", diffNormL2 = " << diffNormL2 << endl;
            cerr << "error = " << error << endl;
            cerr << "Error bound: " << errorBound << " was surpassed for the l2 norm of the difference of the basis matrices." << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Finalize the MPI environment.
    MPI_Finalize();
}

int main (int argc, char *argv[]) {
    string baselinePath((string) argv[1]);
    string targetPath((string) argv[2]);
    double errorBound = stod(argv[3]);
    int numProcessors = stoi(argv[4]);

    compareBasis(baselinePath, targetPath, errorBound, numProcessors);
    return 0;

}
