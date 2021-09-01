/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: A class that reads basis vectors from a file.

#include "BasisReader.h"
#include "HDFDatabase.h"
#include "Matrix.h"
#include "Vector.h"
#include "mpi.h"

namespace CAROM {

BasisReader::BasisReader(
    const std::string& base_file_name,
    Database::formats db_format) :
    d_spatial_basis_vectors(NULL),
    d_temporal_basis_vectors(0),
    d_singular_values(0),
    d_last_basis_idx(-1),
    d_snapshots(0),
    full_file_name(""),
    base_file_name_(base_file_name)
{
    CAROM_ASSERT(!base_file_name.empty());

    int mpi_init;
    MPI_Initialized(&mpi_init);
    int rank;
    if (mpi_init) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
    else {
        rank = 0;
    }

    char tmp[100];
    sprintf(tmp, ".%06d", rank);
    full_file_name = base_file_name + tmp;
    if (db_format == Database::HDF5) {
        d_database = new HDFDatabase();
    }

    std::cout << "Opening file: " << full_file_name << std::endl;
    d_database->open(full_file_name, "r");

    int num_time_intervals;
    double foo;
    d_database->getDouble("num_time_intervals", foo);
    num_time_intervals = static_cast<int>(foo);
    d_time_interval_start_times.resize(num_time_intervals);
    for (int i = 0; i < num_time_intervals; ++i) {
        sprintf(tmp, "time_%06d", i);
        d_database->getDouble(tmp, d_time_interval_start_times[i]);
    }
}

BasisReader::~BasisReader()
{
    delete d_spatial_basis_vectors;
    delete d_temporal_basis_vectors;
    delete d_singular_values;
    d_database->close();
    delete d_database;
    delete d_snapshots;
}

void
BasisReader::readBasis(
    const std::string& base_file_name,
    Database::formats db_format)
{
    CAROM_ASSERT(!base_file_name_.empty());

    int mpi_init;
    MPI_Initialized(&mpi_init);
    int rank;
    if (mpi_init) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
    else {
        rank = 0;
    }

    if (db_format == Database::HDF5) {
        delete d_database;
        d_database = new HDFDatabase();
    }

    char tmp[100];
    sprintf(tmp, ".%06d", rank);
    if (base_file_name_ != base_file_name && !base_file_name.empty()) {
        full_file_name = base_file_name + tmp;
    }

    std::cout << "Opening file: " << full_file_name << std::endl;
    d_database->open(full_file_name, "r");

    int num_time_intervals;
    double foo;
    d_database->getDouble("num_time_intervals", foo);
    num_time_intervals = static_cast<int>(foo);
    d_time_interval_start_times.resize(num_time_intervals);
    for (int i = 0; i < num_time_intervals; ++i) {
        sprintf(tmp, "time_%06d", i);
        d_database->getDouble(tmp, d_time_interval_start_times[i]);
    }
}

const Matrix*
BasisReader::getSpatialBasis(
    double time)
{
    CAROM_ASSERT(0 < numTimeIntervals());
    CAROM_ASSERT(0 <= time);
    int num_time_intervals = numTimeIntervals();
    int i;
    for (i = 0; i < num_time_intervals-1; ++i) {
        if (d_time_interval_start_times[i] <= time &&
                time < d_time_interval_start_times[i+1]) {
            break;
        }
    }
    d_last_basis_idx = i;
    char tmp[100];
    int num_rows;
    sprintf(tmp, "spatial_basis_num_rows_%06d", i);
    d_database->getInteger(tmp, num_rows);
    int num_cols;
    sprintf(tmp, "spatial_basis_num_cols_%06d", i);
    d_database->getInteger(tmp, num_cols);
    if (d_spatial_basis_vectors) {
        delete d_spatial_basis_vectors;
    }
    d_spatial_basis_vectors = new Matrix(num_rows, num_cols, true);
    sprintf(tmp, "spatial_basis_%06d", i);
    d_database->getDoubleArray(tmp,
                               &d_spatial_basis_vectors->item(0, 0),
                               num_rows*num_cols);
    return d_spatial_basis_vectors;
}

const Matrix*
BasisReader::getSpatialBasis(
    double time,
    int n)
{
    CAROM_ASSERT(0 < numTimeIntervals());
    CAROM_ASSERT(0 <= time);
    int num_time_intervals = numTimeIntervals();
    int i;
    for (i = 0; i < num_time_intervals-1; ++i) {
        if (d_time_interval_start_times[i] <= time &&
                time < d_time_interval_start_times[i+1]) {
            break;
        }
    }
    d_last_basis_idx = i;
    char tmp[100];
    int num_rows;
    sprintf(tmp, "spatial_basis_num_rows_%06d", i);
    d_database->getInteger(tmp, num_rows);
    int num_cols;
    sprintf(tmp, "spatial_basis_num_cols_%06d", i);
    d_database->getInteger(tmp, num_cols);
    if (d_spatial_basis_vectors) {
        delete d_spatial_basis_vectors;
    }
    CAROM_VERIFY(0 < n <= num_cols);
    d_spatial_basis_vectors = new Matrix(num_rows, n, true);
    sprintf(tmp, "spatial_basis_%06d", i);
    d_database->getDoubleArray(tmp,
                               &d_spatial_basis_vectors->item(0, 0),
                               num_rows*n,
                               n,
                               num_cols - n);
    return d_spatial_basis_vectors;
}

const Matrix*
BasisReader::getSpatialBasis(
    double time,
    double ef)
{
    const Vector* sv = getSingularValues(time);
    double total_energy = 0.0;
    double energy = 0.0;
    for (int i = 0; i < sv->dim(); i++)
    {
        total_energy += sv->item(i);
    }

    int num_used_singular_values = 0;
    for (int i = 0; i < sv->dim(); i++)
    {
        energy += sv->item(i);
        num_used_singular_values++;
        if (energy >= ef)
        {
            break;
        }
    }

    delete sv;
    return getSpatialBasis(time, num_used_singular_values);
}

const Matrix*
BasisReader::getTemporalBasis(
    double time)
{
    CAROM_ASSERT(0 < numTimeIntervals());
    CAROM_ASSERT(0 <= time);
    int num_time_intervals = numTimeIntervals();
    int i;
    for (i = 0; i < num_time_intervals-1; ++i) {
        if (d_time_interval_start_times[i] <= time &&
                time < d_time_interval_start_times[i+1]) {
            break;
        }
    }
    d_last_basis_idx = i;
    char tmp[100];
    int num_rows;
    sprintf(tmp, "temporal_basis_num_rows_%06d", i);
    d_database->getInteger(tmp, num_rows);
    int num_cols;
    sprintf(tmp, "temporal_basis_num_cols_%06d", i);
    d_database->getInteger(tmp, num_cols);
    if (d_temporal_basis_vectors) {
        delete d_temporal_basis_vectors;
    }
    d_temporal_basis_vectors = new Matrix(num_rows, num_cols, true);
    sprintf(tmp, "temporal_basis_%06d", i);
    d_database->getDoubleArray(tmp,
                               &d_temporal_basis_vectors->item(0, 0),
                               num_rows*num_cols);
    return d_temporal_basis_vectors;
}

const Matrix*
BasisReader::getTemporalBasis(
    double time,
    int n)
{
    CAROM_ASSERT(0 < numTimeIntervals());
    CAROM_ASSERT(0 <= time);
    int num_time_intervals = numTimeIntervals();
    int i;
    for (i = 0; i < num_time_intervals-1; ++i) {
        if (d_time_interval_start_times[i] <= time &&
                time < d_time_interval_start_times[i+1]) {
            break;
        }
    }
    d_last_basis_idx = i;
    char tmp[100];
    int num_rows;
    sprintf(tmp, "temporal_basis_num_rows_%06d", i);
    d_database->getInteger(tmp, num_rows);
    int num_cols;
    sprintf(tmp, "temporal_basis_num_cols_%06d", i);
    d_database->getInteger(tmp, num_cols);
    if (d_temporal_basis_vectors) {
        delete d_temporal_basis_vectors;
    }
    CAROM_VERIFY(0 < n <= num_cols);
    d_temporal_basis_vectors = new Matrix(num_rows, n, true);
    sprintf(tmp, "temporal_basis_%06d", i);
    d_database->getDoubleArray(tmp,
                               &d_temporal_basis_vectors->item(0, 0),
                               num_rows*n,
                               n,
                               num_cols - n);
    return d_temporal_basis_vectors;
}

const Matrix*
BasisReader::getTemporalBasis(
    double time,
    double ef)
{
    const Vector* sv = getSingularValues(time);
    double total_energy = 0.0;
    double energy = 0.0;
    for (int i = 0; i < sv->dim(); i++)
    {
        total_energy += sv->item(i);
    }

    int num_used_singular_values = 0;
    for (int i = 0; i < sv->dim(); i++)
    {
        energy += sv->item(i);
        num_used_singular_values++;
        if (energy >= ef)
        {
            break;
        }
    }

    delete sv;
    return getTemporalBasis(time, num_used_singular_values);
}

const Vector*
BasisReader::getSingularValues(
    double time)
{
    CAROM_ASSERT(0 < numTimeIntervals());
    CAROM_ASSERT(0 <= time);
    int num_time_intervals = numTimeIntervals();
    int i;
    for (i = 0; i < num_time_intervals-1; ++i) {
        if (d_time_interval_start_times[i] <= time &&
                time < d_time_interval_start_times[i+1]) {
            break;
        }
    }
    d_last_basis_idx = i;
    char tmp[100];
    int size;
    sprintf(tmp, "singular_value_size_%06d", i);
    d_database->getInteger(tmp, size);
    if (d_singular_values) {
        delete d_singular_values;
    }
    d_singular_values = new Vector(size, false);
    sprintf(tmp, "singular_value_%06d", i);
    d_database->getDoubleArray(tmp,
                               &d_singular_values->item(0),
                               size);
    return d_singular_values;
}

const Vector*
BasisReader::getSingularValues(
    double time,
    double ef)
{
    const Vector* sv = getSingularValues(time);
    double total_energy = 0.0;
    double energy = 0.0;
    for (int i = 0; i < sv->dim(); i++)
    {
        total_energy += sv->item(i);
    }

    int num_used_singular_values = 0;
    for (int i = 0; i < sv->dim(); i++)
    {
        energy += sv->item(i);
        num_used_singular_values++;
        if (energy >= ef)
        {
            break;
        }
    }

    Vector* truncated_sv = new Vector(num_used_singular_values, false);
    for (int i = 0; i < num_used_singular_values; i++)
    {
        truncated_sv->item(i) = sv->item(i);
    }

    delete sv;
    return truncated_sv;
}

const Matrix*
BasisReader::getSnapshotMatrix(
    double time)
{
    CAROM_ASSERT(0 < numTimeIntervals());
    CAROM_ASSERT(0 <= time);
    int num_time_intervals = numTimeIntervals();
    int i;
    for (i = 0; i < num_time_intervals-1; ++i) {
        if (d_time_interval_start_times[i] <= time &&
                time < d_time_interval_start_times[i+1]) {
            break;
        }
    }
    d_last_basis_idx = i;
    char tmp[100];
    int num_rows;
    sprintf(tmp, "snapshot_matrix_num_rows_%06d", i);
    d_database->getInteger(tmp, num_rows);
    int num_cols;
    sprintf(tmp, "snapshot_matrix_num_cols_%06d", i);
    d_database->getInteger(tmp, num_cols);
    if (d_snapshots) {
        delete d_snapshots;
    }
    d_snapshots = new Matrix(num_rows, num_cols, false);
    sprintf(tmp, "snapshot_matrix_%06d", i);
    d_database->getDoubleArray(tmp,
                               &d_snapshots->item(0, 0),
                               num_rows*num_cols);
    return d_snapshots;
}

Matrix
BasisReader::getMatlabBasis(
    double time)
{
    CAROM_ASSERT(0 < numTimeIntervals());
    CAROM_ASSERT(0 <= time);
    int num_time_intervals = numTimeIntervals();
    int i;
    for (i = 0; i < num_time_intervals-1; ++i) {
        if (d_time_interval_start_times[i] <= time &&
                time < d_time_interval_start_times[i+1]) {
            break;
        }
    }
    d_last_basis_idx = i;
    char tmp[100];
    int num_rows;
    double foo;
    sprintf(tmp, "num_rows_%06d", i);
    d_database->getDouble(tmp, foo);
    num_rows = static_cast<int>(foo);
    int num_cols;
    sprintf(tmp, "num_cols_%06d", i);
    d_database->getDouble(tmp, foo);
    num_cols = static_cast<int>(foo);
    if (d_spatial_basis_vectors) {
        delete d_spatial_basis_vectors;
    }
    d_spatial_basis_vectors = new Matrix(num_rows, num_cols, true);
    sprintf(tmp, "matlab_spatial_basis_%06d", i);
    d_database->getDoubleArray(tmp,
                               &d_spatial_basis_vectors->item(0, 0),
                               num_rows*num_cols);
    return *d_spatial_basis_vectors;
}
}
