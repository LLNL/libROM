/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: This function reads snapshot or basis files and computes the
//		SVD over the group. Also provides the option to subtract the
//		mean or a supplied offset vector before computing the SVD.
//
// Assumptions: You are running this file with the same number of processors
//		as when the snapshots/bases were saved.
//
// Inputs:
//	[1] List of file names without extension.
//	[2] Optional: "-b" or "basis" to compute SVD over bases instead
//	    of snapshots. (Snapshots are default.)
// 	[3] Optional: "-m" or "mean" to subtract the mean before
//	    computing the SVD. (No subtraction is default.)
//      [4] Optional: "-o" or "offset" to subtract an arbitrary offset vector
//          before computing the SVD. Offset vector file supplied as next argument.
//
// Outputs:
//	[1] total_snapshot.*:           single hdf5 file (per rank) of all loaded data.
//	[2] mean.*:                     single hdf5 file (per rank) with subtracted mean.
//	[3] total.* OR total_meansub.*: single hdf5 file (per rank) with the new POD basis.
//
// Examples: mpirun -n 2 ./combine_samples file1 file2 file3 -b -m
//           mpirun -n 2 ./combine_samples -o offsetfile file1 file2 file3

#include "linalg/BasisGenerator.h"
#include "linalg/scalapack_wrapper.h"
#include "linalg/BasisReader.h"
#include "linalg/BasisWriter.h"
#include <vector>
#include "mpi.h"
#include <string>
#include <stdio.h>

int main(int argc, char* argv[])
{
    // Initialize MPI
    MPI_Init(&argc, &argv);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<std::string> sample_names;
    int snaps = 0;
    int dim   = 0;
    bool subtract_mean = false;
    bool subtract_offset = false;
    std::string kind = "snapshot";
    std::string offset_file;

    bool offset_arg = false;

    if (argc >= 2) {
        for (int i = 1; i < argc; i++) {
            if (!strcmp(argv[i], "basis") || !strcmp(argv[i], "-b")) {
                if (rank==0) std::cout << "Argument " << i << " identified as basis or -b" <<
                                           std::endl;
                kind = "basis";
            }
            else if (offset_arg) {
                offset_file = argv[i];
                offset_arg = false;
            }
            else if (!strcmp(argv[i], "offset") || !strcmp(argv[i], "-o")) {
                if (rank==0) std::cout << "Will subtract offset" << std::endl;
                subtract_offset = true;
                offset_arg = true;
            }
            else if (!strcmp(argv[i], "mean") || !strcmp(argv[i], "-m")) {
                if (rank==0) std::cout << "Will subtract mean" << std::endl;
                subtract_mean = true;
            }
            else {
                sample_names.push_back(argv[i]);
            }
        }
    }
    else {
        if (rank==0) std::cout << "No arguments passed." << std::endl;
        return 1;
    }

    CAROM_VERIFY(!(subtract_mean
                   && subtract_offset)); // both should not be selected

    /*-- Read dimension and count number of snapshots/bases --*/
    if (rank==0) std::cout <<
                               "Opening files to read dimension and count number of snapshots/bases" <<
                               std::endl;
    int dimFirst = 0;
    for (const auto& sample_name: sample_names) {
        CAROM::BasisReader reader(sample_name);

        dim    = reader.getDim(kind, 0);
        snaps += reader.getNumSamples(kind, 0);
        if (dimFirst == 0) dimFirst = dim;

        CAROM_VERIFY(dim ==
                     dimFirst); // ensures all files have the same dimension distribution
    }

    CAROM_VERIFY((snaps > 0) && (dim > 0));

    /*-- Load data from input files --*/
    std::string generator_filename = "total";
    std::unique_ptr<CAROM::BasisGenerator> static_basis_generator;
    static_basis_generator.reset(new CAROM::BasisGenerator(
                                     CAROM::Options(dim, snaps).setMaxBasisDimension(snaps), false,
                                     generator_filename));

    if (rank==0) std::cout << "Loading data from " << kind << std::endl;
    for(const auto& sample_name: sample_names) {
        static_basis_generator->loadSamples(sample_name, kind);
    }

    if (rank==0) std::cout << "Saving data uploaded as a snapshot matrix" <<
                               std::endl;
    static_basis_generator->writeSnapshot();

    if (!subtract_mean && !subtract_offset) {
        /*-- Compute SVD and save file --*/
        if (rank==0) std::cout << "Computing SVD" << std::endl;
        int rom_dim = static_basis_generator->getSpatialBasis()->numColumns();
        if (rank==0) std::cout << "U ROM Dimension: " << rom_dim << std::endl;
        static_basis_generator->endSamples();
    }
    else {
        /*-- load data from hdf5 file to find the mean and subtract it --*/
        if (rank==0) std::cout << "Reading snapshots" << std::endl;
        CAROM::Matrix* snapshots = (CAROM::Matrix*)
                                   static_basis_generator->getSnapshotMatrix();

        int num_rows = snapshots->numRows();
        int num_cols = snapshots->numColumns();
        CAROM::Vector offset(num_rows, true);

        if (subtract_mean) {
            /*-- Find the mean per row and write to hdf5 --*/
            if (rank==0) std::cout << "Subtracting mean" << std::endl;
            for (int row = 0; row < num_rows; ++row) {
                offset(row) = 0;
                for (int col = 0; col < num_cols; ++col) {
                    offset(row) += snapshots->item(row,col);
                }
            }
            offset.mult(1.0 / (double)num_cols, offset);
            offset.write("mean");
        }
        else if (subtract_offset) {
            if (rank==0) std::cout << "Subtracting offset from file: " << offset_file <<
                                       std::endl;
            offset.read(offset_file);
        }

        /*-- Subtract mean or offset from snapshot/bases matrix --*/
        for (int row = 0; row < num_rows; ++row) {
            for (int col = 0; col < num_cols; ++col) {
                snapshots->item(row,col) = snapshots->item(row,col) - offset(row);
            }
        }

        /*-- Load new samples with mean/offset subtracted --*/
        // TODO: Add capability to librom for editing an existing snapshot matrix,
        // so second generator is not needed.
        std::unique_ptr<CAROM::BasisGenerator> static_basis_generator2;
        static_basis_generator2.reset(new CAROM::BasisGenerator(
                                          CAROM::Options(dim, snaps).setMaxBasisDimension(snaps), false,
                                          generator_filename+"_sub"));

        CAROM::Vector snap_cur(num_rows, true);
        for (int col = 0; col < num_cols; col++) {
            snap_cur = *snapshots->getColumn(col);
            static_basis_generator2->takeSample(snap_cur.getData(), 0.0, false);
        }

        /*-- Compute SVD and save file --*/
        if (rank==0) std::cout << "Computing SVD" << std::endl;
        int rom_dim = static_basis_generator2->getSpatialBasis()->numColumns();
        if (rank==0) std::cout << "U ROM Dimension: " << rom_dim << std::endl;
        static_basis_generator2->endSamples();

        static_basis_generator2 = nullptr;
    }
    static_basis_generator = nullptr;

    MPI_Finalize();
    return 0;
}
