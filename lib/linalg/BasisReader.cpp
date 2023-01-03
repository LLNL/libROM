/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: A class that reads basis vectors from a file.

#include "BasisReader.h"
#include "utils/HDFDatabase.h"
#include "utils/CSVDatabase.h"
#include "Matrix.h"
#include "Vector.h"
#include "mpi.h"

namespace CAROM {

BasisReader::BasisReader(
    const std::string& base_file_name,
    Database::formats db_format) :
    d_last_basis_idx(-1),
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
    else if (db_format == Database::CSV) {
        d_database = new CSVDatabase();
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
    d_database->close();
    delete d_database;
}

Matrix*
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
    int num_rows = getDim("basis",time);
    int num_cols = getNumSamples("basis",time);

    char tmp[100];
    Matrix* spatial_basis_vectors = new Matrix(num_rows, num_cols, true);
    sprintf(tmp, "spatial_basis_%06d", i);
    d_database->getDoubleArray(tmp,
                               &spatial_basis_vectors->item(0, 0),
                               num_rows*num_cols);
    return spatial_basis_vectors;
}

Matrix*
BasisReader::getSpatialBasis(
    double time,
    int n)
{
    return getSpatialBasis(time, 1, n);
}

Matrix*
BasisReader::getSpatialBasis(
    double time,
    int start_col,
    int end_col)
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
    int num_rows = getDim("basis",time);
    int num_cols = getNumSamples("basis",time);

    char tmp[100];
    CAROM_VERIFY(0 < start_col <= num_cols);
    CAROM_VERIFY(start_col <= end_col && end_col <= num_cols);
    int num_cols_to_read = end_col - start_col + 1;

    Matrix* spatial_basis_vectors = new Matrix(num_rows, num_cols_to_read, true);
    sprintf(tmp, "spatial_basis_%06d", i);
    d_database->getDoubleArray(tmp,
                               &spatial_basis_vectors->item(0, 0),
                               num_rows*num_cols_to_read,
                               start_col - 1,
                               num_cols_to_read,
                               num_cols);
    return spatial_basis_vectors;
}

Matrix*
BasisReader::getSpatialBasis(
    double time,
    double ef)
{
    Vector* sv = getSingularValues(time);
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
        if (energy / total_energy >= ef)
        {
            break;
        }
    }

    delete sv;
    return getSpatialBasis(time, num_used_singular_values);
}

Matrix*
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
    int num_rows = getDim("temporal_basis",time);
    int num_cols = getNumSamples("temporal_basis",time);

    char tmp[100];
    Matrix* temporal_basis_vectors = new Matrix(num_rows, num_cols, true);
    sprintf(tmp, "temporal_basis_%06d", i);
    d_database->getDoubleArray(tmp,
                               &temporal_basis_vectors->item(0, 0),
                               num_rows*num_cols);
    return temporal_basis_vectors;
}

Matrix*
BasisReader::getTemporalBasis(
    double time,
    int n)
{
    return getTemporalBasis(time, 1, n);
}

Matrix*
BasisReader::getTemporalBasis(
    double time,
    int start_col,
    int end_col)
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
    int num_rows = getDim("temporal_basis",time);
    int num_cols = getNumSamples("temporal_basis",time);

    char tmp[100];
    CAROM_VERIFY(0 < start_col <= num_cols);
    CAROM_VERIFY(start_col <= end_col && end_col <= num_cols);
    int num_cols_to_read = end_col - start_col + 1;

    Matrix* temporal_basis_vectors = new Matrix(num_rows, num_cols_to_read, true);
    sprintf(tmp, "temporal_basis_%06d", i);
    d_database->getDoubleArray(tmp,
                               &temporal_basis_vectors->item(0, 0),
                               num_rows*num_cols_to_read,
                               start_col - 1,
                               num_cols_to_read,
                               num_cols);
    return temporal_basis_vectors;
}

Matrix*
BasisReader::getTemporalBasis(
    double time,
    double ef)
{
    Vector* sv = getSingularValues(time);
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
        if (energy / total_energy >= ef)
        {
            break;
        }
    }

    delete sv;
    return getTemporalBasis(time, num_used_singular_values);
}

Vector*
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

    Vector* singular_values = new Vector(size, false);
    sprintf(tmp, "singular_value_%06d", i);
    d_database->getDoubleArray(tmp,
                               &singular_values->item(0),
                               size);
    return singular_values;
}

Vector*
BasisReader::getSingularValues(
    double time,
    double ef)
{
    Vector* sv = getSingularValues(time);
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
        if (energy / total_energy >= ef)
        {
            break;
        }
    }

    Vector* truncated_singular_values = new Vector(num_used_singular_values,
            false);
    for (int i = 0; i < num_used_singular_values; i++)
    {
        truncated_singular_values->item(i) = sv->item(i);
    }

    delete sv;
    return truncated_singular_values;
}

int
BasisReader::getDim(
    const std::string kind,
    double time)
{
    CAROM_ASSERT(0 < numTimeIntervals());
    CAROM_ASSERT(0 <= time);
    CAROM_ASSERT((kind == "basis") ||
                 (kind == "snapshot") ||
                 (kind == "temporal_basis"));

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
    if (kind == "basis") sprintf(tmp, "spatial_basis_num_rows_%06d", i);
    else if (kind == "snapshot") sprintf(tmp, "snapshot_matrix_num_rows_%06d", i);
    else if (kind == "temporal_basis") sprintf(tmp, "temporal_basis_num_rows_%06d",
                i);
    d_database->getInteger(tmp, num_rows);
    return num_rows;
}

int
BasisReader::getNumSamples(
    const std::string kind,
    double time)
{
    CAROM_ASSERT(0 < numTimeIntervals());
    CAROM_ASSERT(0 <= time);
    CAROM_ASSERT((kind == "basis") ||
                 (kind == "snapshot") ||
                 (kind == "temporal_basis"));

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
    int num_cols;
    if (kind == "basis") sprintf(tmp, "spatial_basis_num_cols_%06d", i);
    else if (kind == "snapshot") sprintf(tmp, "snapshot_matrix_num_cols_%06d", i);
    else if (kind == "temporal_basis") sprintf(tmp, "temporal_basis_num_cols_%06d",
                i);
    d_database->getInteger(tmp, num_cols);
    return num_cols;
}

Matrix*
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
    int num_rows = getDim("snapshot",time);
    int num_cols = getNumSamples("snapshot",time);

    char tmp[100];
    Matrix* snapshots = new Matrix(num_rows, num_cols, false);
    sprintf(tmp, "snapshot_matrix_%06d", i);
    d_database->getDoubleArray(tmp,
                               &snapshots->item(0, 0),
                               num_rows*num_cols);
    return snapshots;
}

Matrix*
BasisReader::getSnapshotMatrix(
    double time,
    int n)
{
    return getSnapshotMatrix(time, 1, n);
}

Matrix*
BasisReader::getSnapshotMatrix(
    double time,
    int start_col,
    int end_col)
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
    int num_rows = getDim("snapshot",time);
    int num_cols = getNumSamples("snapshot",time);

    CAROM_VERIFY(0 < start_col <= num_cols);
    CAROM_VERIFY(start_col <= end_col && end_col <= num_cols);
    int num_cols_to_read = end_col - start_col + 1;

    char tmp[100];
    Matrix* snapshots = new Matrix(num_rows, num_cols_to_read, false);
    sprintf(tmp, "snapshot_matrix_%06d", i);
    d_database->getDoubleArray(tmp,
                               &snapshots->item(0, 0),
                               num_rows*num_cols_to_read,
                               start_col - 1,
                               num_cols_to_read,
                               num_cols);
    return snapshots;
}
}
