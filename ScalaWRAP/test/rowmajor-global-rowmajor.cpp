/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

#include "catch.hpp"
#include "test-common.hpp"
#include "../src/LocalMatrix.hpp"
#include "../src/ScalaMat.hpp"
#include "../src/MPI_utils.hpp"
#include "../src/C_interface.h"

#include <iostream>
#include <memory>
#include <random>
#include <string>

#define NOISY_TESTS false

using namespace ScalaWRAP;
using std::cout;

TEST_CASE("Copying row major array to distributed and back to row major",
          "[data-movement]")
{
    MPI_Barrier(MPI_COMM_WORLD);

    // Setup here gets executed for each SECTION.
    constexpr int global_m = 20, global_n = 8;
    std::unique_ptr<double[]> local_data(new double[global_m * global_n]);
    for (int j = 0; j < global_n; ++j)
    {
        for (int i = 0; i < global_m; ++i)
        {
            local_data[j * global_m + i] = (mpi_rank() + 1) * (j * global_m + i);
        }
    }

    for (int rank = 0; rank < mpi_size(); ++rank)
    {
        SECTION(std::string("default_copy_from") + std::to_string(rank))
        {
            MPI_Barrier(MPI_COMM_WORLD);
            ScalaMat A(global_m, global_n);
            if (mpi_rank() == rank && NOISY_TESTS)
            {
                cout << "Testing row major copy with default blocks from rank "
                     << rank << '\n';
                cout << "Local matrix on rank " << rank << ":\n";
                print_local_mat(local_data.get(), global_m, global_n, ROW_MAJOR);
                cout << "\n\n";
            }
            LocalMatrix localmat(local_data.get(), global_m, global_n, rank, false, ROW_MAJOR);

            A.submatrix(rowrange(1, global_m), colrange(1, global_n)) = std::move(localmat);
            if (mpi_rank() == 0 && NOISY_TESTS)
                cout << "Distributed matrix:\n";
            mpi_foreach(
                [&A]() { if (NOISY_TESTS) {
                            A.dump_debug_info(cout, true); cout.flush();
                        } });
            if (mpi_rank() == 0 && NOISY_TESTS)
                cout << "\n\n";

            std::unique_ptr<double[]> diff(new double[global_m * global_n]);
            LocalMatrix(diff.get(), global_m, global_n, rank, false, ROW_MAJOR) =
                A(rowrange(1, global_m), colrange(1, global_n));
            if (mpi_rank() == rank && NOISY_TESTS)
            {
                cout << "Local matrix as retrieved on rank " << rank << ":\n";
                print_local_mat(diff.get(), global_m, global_n, ROW_MAJOR);
                cout << "\n\n";
            }

            double squared_diff = 0.0;
            if (mpi_rank() == rank)
            {
                for (int j = 0; j < global_m * global_n; ++j)
                {
                    diff[j] = local_data[j] - diff[j];
                    squared_diff += diff[j] * diff[j];
                }
            }
            MPI_Bcast(&squared_diff, 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);

            REQUIRE(squared_diff < 1e-12);
        }
    }

    for (int rank = 0; rank < mpi_size(); ++rank)
    {
        SECTION(std::string("smallblock_copy_from" + std::to_string(rank)))
        {
            MPI_Barrier(MPI_COMM_WORLD);
            ScalaMat A(global_m, global_n, 2, 2);
            if (mpi_rank() == rank && NOISY_TESTS)
            {
                cout << "Testing copy with 2 x 2 blocks from " << rank << '\n';
                cout << "Testing copy from row major to global back to row major\n";
                cout << "Local matrix on rank " << rank << ":\n";
                print_local_mat(local_data.get(), global_m, global_n, ROW_MAJOR);
                cout << "\n\n";
            }
            LocalMatrix localmat(local_data.get(), global_m, global_n, rank, false, ROW_MAJOR);

            A.submatrix(rowrange(1, global_m), colrange(1, global_n)) = std::move(localmat);
            if (mpi_rank() == 0 && NOISY_TESTS)
                cout << "Distributed matrix:\n";
            mpi_foreach(
                [&A]() {
                    if (NOISY_TESTS) {
                        A.dump_debug_info(cout, true); cout.flush();
                    } });
            if (mpi_rank() == 0 && NOISY_TESTS)
                cout << "\n\n";

            std::unique_ptr<double[]> diff(new double[global_m * global_n]);
            LocalMatrix(diff.get(), global_m, global_n, rank, false, ROW_MAJOR) =
                A(rowrange(1, global_m), colrange(1, global_n));
            if (mpi_rank() == rank && NOISY_TESTS)
            {
                cout << "Local matrix as retrieved on rank " << rank << ":\n";
                print_local_mat(local_data.get(), global_m, global_n, ROW_MAJOR);
                cout << "\n\n";
            }

            double squared_diff = 0.0;
            if (mpi_rank() == rank)
            {
                for (int j = 0; j < global_m * global_n; ++j)
                {
                    diff[j] = local_data[j] - diff[j];
                    squared_diff += diff[j] * diff[j];
                }
            }
            MPI_Bcast(&squared_diff, 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);

            REQUIRE(squared_diff < 1e-12);
        }
    }

    for (int rank = 0; rank < mpi_size(); ++rank)
    {
        SECTION(std::string("smallblock_copy_w_offsets_from" + std::to_string(rank)))
        {
            MPI_Barrier(MPI_COMM_WORLD);
            int pi, pj;
            std::tie(pi, pj) = get_random_coord(Context::default_context());
            ScalaMat A(global_m, global_n, 2, 2, Context::default_context(), pi, pj);
            if (mpi_rank() == rank && NOISY_TESTS)
            {
                cout << "Testing copy with 2 x 2 blocks from " << rank << '\n';
                cout << "Testing copy from row major to global back to row major\n";
                cout << "Testing with random rowsrc and colsrc\n";
                cout << "Local matrix on rank " << rank << ":\n";
                print_local_mat(local_data.get(), global_m, global_n, ROW_MAJOR);
                cout << "\n\n";
            }
            LocalMatrix localmat(local_data.get(), global_m, global_n, rank, false, ROW_MAJOR);

            A.submatrix(rowrange(1, global_m), colrange(1, global_n)) = std::move(localmat);
            if (mpi_rank() == 0 && NOISY_TESTS)
                cout << "Distributed matrix:\n";
            mpi_foreach(
                [&A]() {
                    if (NOISY_TESTS) {
                        A.dump_debug_info(cout, true); cout.flush();
                    } });
            if (mpi_rank() == 0 && NOISY_TESTS)
                cout << "\n\n";

            std::unique_ptr<double[]> diff(new double[global_m * global_n]);
            LocalMatrix(diff.get(), global_m, global_n, rank, false, ROW_MAJOR) =
                A(rowrange(1, global_m), colrange(1, global_n));
            if (mpi_rank() == rank && NOISY_TESTS)
            {
                cout << "Local matrix as retrieved on rank " << rank << ":\n";
                print_local_mat(local_data.get(), global_m, global_n, ROW_MAJOR);
                cout << "\n\n";
            }

            double squared_diff = 0.0;
            if (mpi_rank() == rank)
            {
                for (int j = 0; j < global_m * global_n; ++j)
                {
                    diff[j] = local_data[j] - diff[j];
                    squared_diff += diff[j] * diff[j];
                }
            }
            MPI_Bcast(&squared_diff, 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);

            REQUIRE(squared_diff < 1e-12);
        }
    }

    for (int rank = 0; rank < mpi_size(); ++rank)
    {
        SECTION(std::string("randblock_copy_w_offset_from_" + std::to_string(rank)))
        {
            MPI_Barrier(MPI_COMM_WORLD);
            int pi, pj;
            std::tie(pi, pj) = get_random_coord(Context::default_context());
            int bs[2] = { 0 };
            if (mpi_rank() == 0) {
                bs[0] = random_between(1, global_m);
                bs[1] = random_between(1, global_n);
            }
            MPI_Bcast(bs, 2, MPI_INT, 0, MPI_COMM_WORLD);

            ScalaMat A(global_m, global_n, bs[0], bs[1], Context::default_context(), pi, pj);
            if (mpi_rank() == rank && NOISY_TESTS)
            {
                cout << "Testing copy with random blocksize from " << rank << '\n';
                cout << "Testing copy from row major to global back to row major\n";
                cout << "Testing with random rowsrc and colsrc\n";
                cout << "Local matrix on rank " << rank << ":\n";
                print_local_mat(local_data.get(), global_m, global_n, ROW_MAJOR);
                cout << "\n\n";
            }
            LocalMatrix localmat(local_data.get(), global_m, global_n, rank, false, ROW_MAJOR);

            A.submatrix(rowrange(1, global_m), colrange(1, global_n)) = std::move(localmat);
            if (mpi_rank() == 0 && NOISY_TESTS)
                cout << "Distributed matrix:\n";
            mpi_foreach(
                [&A]() {
                    if (NOISY_TESTS) {
                        A.dump_debug_info(cout, true); cout.flush();
                    } });
            if (mpi_rank() == 0 && NOISY_TESTS)
                cout << "\n\n";

            std::unique_ptr<double[]> diff(new double[global_m * global_n]);
            LocalMatrix(diff.get(), global_m, global_n, rank, false, ROW_MAJOR) =
                A(rowrange(1, global_m), colrange(1, global_n));
            if (mpi_rank() == rank && NOISY_TESTS)
            {
                cout << "Local matrix as retrieved on rank " << rank << ":\n";
                print_local_mat(local_data.get(), global_m, global_n, ROW_MAJOR);
                cout << "\n\n";
            }

            double squared_diff = 0.0;
            if (mpi_rank() == rank)
            {
                for (int j = 0; j < global_m * global_n; ++j)
                {
                    diff[j] = local_data[j] - diff[j];
                    squared_diff += diff[j] * diff[j];
                }
            }
            MPI_Bcast(&squared_diff, 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);

            REQUIRE(squared_diff < 1e-12);
        }
    }

    for (int rank = 0; rank < mpi_size(); ++rank)
    {
        SECTION(std::string("randblock_copy_w_offset_and_random_map_from_" + std::to_string(rank)))
        {
            MPI_Barrier(MPI_COMM_WORLD);
            // Create a random process map with a subset of processes
            int nsub;
            if (mpi_rank() == 0)
                nsub = random_between(1, mpi_size());
            MPI_Bcast(&nsub, 1, MPI_INT, 0, MPI_COMM_WORLD);

            std::vector<int> pids;
            pids.resize(nsub);
            if (mpi_rank() == 0) {
                pids.resize(mpi_size());
                for (unsigned i = 0; i < mpi_size(); ++i)
                    pids[i] = i;
                std::shuffle(pids.begin(), pids.end(), std::random_device());
                pids.resize(nsub);
            }
            MPI_Bcast(pids.data(), pids.size(), MPI_INT, 0, MPI_COMM_WORLD);
            int nprow = best_distribution(nsub).first;

            auto ctxt = Context::make_context(pids, nprow);

            int pi, pj;
            int rank = pids[0];
            std::tie(pi, pj) = get_random_coord(ctxt, rank);
            int bs[2] = { 0 };
            if (mpi_rank() == 0) {
                bs[0] = random_between(1, global_m);
                bs[1] = random_between(1, global_n);
            }
            MPI_Bcast(bs, 2, MPI_INT, 0, MPI_COMM_WORLD);

            ScalaMat A(global_m, global_n, bs[0], bs[1], ctxt, pi, pj);
            if (mpi_rank() == rank && NOISY_TESTS)
            {
                cout << "Testing copy with random blocksize from " << rank << '\n';
                cout << "Testing copy from row major to global back to row major\n";
                cout << "Testing with random rowsrc and colsrc\n";
                cout << "Testing with a random process map of a random subset\n";
                cout << "Local matrix on rank " << rank << ":\n";
                print_local_mat(local_data.get(), global_m, global_n, ROW_MAJOR);
                cout << "\n\n";
            }
            LocalMatrix localmat(local_data.get(), global_m, global_n, rank, false, ROW_MAJOR);

            A.submatrix(rowrange(1, global_m), colrange(1, global_n)) = std::move(localmat);
            if (mpi_rank() == 0 && NOISY_TESTS)
                cout << "Distributed matrix:\n";
            mpi_foreach(
                [&A]() {
                    if (NOISY_TESTS) {
                        A.dump_debug_info(cout, true); cout.flush();
                    } });
            if (mpi_rank() == 0 && NOISY_TESTS)
                cout << "\n\n";

            std::unique_ptr<double[]> diff(new double[global_m * global_n]);
            LocalMatrix(diff.get(), global_m, global_n, rank, false, ROW_MAJOR) =
                A(rowrange(1, global_m), colrange(1, global_n));
            if (mpi_rank() == rank && NOISY_TESTS)
            {
                cout << "Local matrix as retrieved on rank " << rank << ":\n";
                print_local_mat(local_data.get(), global_m, global_n, ROW_MAJOR);
                cout << "\n\n";
            }

            double squared_diff = 0.0;
            if (mpi_rank() == rank)
            {
                for (int j = 0; j < global_m * global_n; ++j)
                {
                    diff[j] = local_data[j] - diff[j];
                    squared_diff += diff[j] * diff[j];
                }
            }
            MPI_Bcast(&squared_diff, 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);

            REQUIRE(squared_diff < 1e-12);
        }
    }
}

#undef NOISY_TESTS