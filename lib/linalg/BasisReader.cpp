/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
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

    d_database->open(full_file_name, "r");
}

BasisReader::~BasisReader()
{
    d_database->close();
    delete d_database;
}

Matrix*
BasisReader::getSpatialBasis()
{
    d_last_basis_idx = 0;
    int num_rows = getDim("basis");
    int num_cols = getNumSamples("basis");

    Matrix* spatial_basis_vectors = new Matrix(num_rows, num_cols, true);

    // This 0 index is the remainder from time interval concept.
    // This remains only for backward compatibility purpose.
    d_database->getDoubleArray("spatial_basis_000000",
                               &spatial_basis_vectors->item(0, 0),
                               num_rows*num_cols);
    return spatial_basis_vectors;
}

Matrix*
BasisReader::getSpatialBasis(
    int n)
{
    return getSpatialBasis(1, n);
}

Matrix*
BasisReader::getSpatialBasis(
    int start_col,
    int end_col)
{
    d_last_basis_idx = 0;
    int num_rows = getDim("basis");
    int num_cols = getNumSamples("basis");

    char tmp[100];
    CAROM_VERIFY(0 < start_col <= num_cols);
    CAROM_VERIFY(start_col <= end_col && end_col <= num_cols);
    int num_cols_to_read = end_col - start_col + 1;

    Matrix* spatial_basis_vectors = new Matrix(num_rows, num_cols_to_read, true);
    // This 0 index is the remainder from time interval concept.
    // This remains only for backward compatibility purpose.
    sprintf(tmp, "spatial_basis_000000");
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
    double ef)
{
    Vector* sv = getSingularValues();
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
    return getSpatialBasis(num_used_singular_values);
}

Matrix*
BasisReader::getTemporalBasis()
{
    d_last_basis_idx = 0;
    int num_rows = getDim("temporal_basis");
    int num_cols = getNumSamples("temporal_basis");

    char tmp[100];
    Matrix* temporal_basis_vectors = new Matrix(num_rows, num_cols, true);
    // This 0 index is the remainder from time interval concept.
    // This remains only for backward compatibility purpose.
    sprintf(tmp, "temporal_basis_000000");
    d_database->getDoubleArray(tmp,
                               &temporal_basis_vectors->item(0, 0),
                               num_rows*num_cols);
    return temporal_basis_vectors;
}

Matrix*
BasisReader::getTemporalBasis(
    int n)
{
    return getTemporalBasis(1, n);
}

Matrix*
BasisReader::getTemporalBasis(
    int start_col,
    int end_col)
{
    d_last_basis_idx = 0;
    int num_rows = getDim("temporal_basis");
    int num_cols = getNumSamples("temporal_basis");

    char tmp[100];
    CAROM_VERIFY(0 < start_col <= num_cols);
    CAROM_VERIFY(start_col <= end_col && end_col <= num_cols);
    int num_cols_to_read = end_col - start_col + 1;

    Matrix* temporal_basis_vectors = new Matrix(num_rows, num_cols_to_read, true);
    // This 0 index is the remainder from time interval concept.
    // This remains only for backward compatibility purpose.
    sprintf(tmp, "temporal_basis_000000");
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
    double ef)
{
    Vector* sv = getSingularValues();
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
    return getTemporalBasis(num_used_singular_values);
}

Vector*
BasisReader::getSingularValues()
{
    d_last_basis_idx = 0;
    char tmp[100];
    int size;
    // This 0 index is the remainder from time interval concept.
    // This remains only for backward compatibility purpose.
    sprintf(tmp, "singular_value_size_000000");
    d_database->getInteger(tmp, size);

    Vector* singular_values = new Vector(size, false);
    // This 0 index is the remainder from time interval concept.
    // This remains only for backward compatibility purpose.
    sprintf(tmp, "singular_value_000000");
    d_database->getDoubleArray(tmp,
                               &singular_values->item(0),
                               size);
    return singular_values;
}

Vector*
BasisReader::getSingularValues(
    double ef)
{
    Vector* sv = getSingularValues();
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
    const std::string kind)
{
    CAROM_ASSERT((kind == "basis") ||
                 (kind == "snapshot") ||
                 (kind == "temporal_basis"));

    d_last_basis_idx = 0;
    char tmp[100];
    int num_rows;

    // This 0 index is the remainder from time interval concept.
    // This remains only for backward compatibility purpose.
    if (kind == "basis") sprintf(tmp, "spatial_basis_num_rows_000000");
    else if (kind == "snapshot") sprintf(tmp, "snapshot_matrix_num_rows_000000");
    else if (kind == "temporal_basis") sprintf(tmp,
                "temporal_basis_num_rows_000000");

    d_database->getInteger(tmp, num_rows);
    return num_rows;
}

int
BasisReader::getNumSamples(
    const std::string kind)
{
    CAROM_ASSERT((kind == "basis") ||
                 (kind == "snapshot") ||
                 (kind == "temporal_basis"));

    d_last_basis_idx = 0;
    char tmp[100];
    int num_cols;

    // This 0 index is the remainder from time interval concept.
    // This remains only for backward compatibility purpose.
    if (kind == "basis") sprintf(tmp, "spatial_basis_num_cols_000000");
    else if (kind == "snapshot") sprintf(tmp, "snapshot_matrix_num_cols_000000");
    else if (kind == "temporal_basis") sprintf(tmp,
                "temporal_basis_num_cols_000000");

    d_database->getInteger(tmp, num_cols);
    return num_cols;
}

Matrix*
BasisReader::getSnapshotMatrix()
{
    d_last_basis_idx = 0;
    int num_rows = getDim("snapshot");
    int num_cols = getNumSamples("snapshot");

    char tmp[100];
    Matrix* snapshots = new Matrix(num_rows, num_cols, false);
    // This 0 index is the remainder from time interval concept.
    // This remains only for backward compatibility purpose.
    sprintf(tmp, "snapshot_matrix_000000");
    d_database->getDoubleArray(tmp,
                               &snapshots->item(0, 0),
                               num_rows*num_cols);
    return snapshots;
}

Matrix*
BasisReader::getSnapshotMatrix(
    int n)
{
    return getSnapshotMatrix(1, n);
}

Matrix*
BasisReader::getSnapshotMatrix(
    int start_col,
    int end_col)
{
    d_last_basis_idx = 0;
    int num_rows = getDim("snapshot");
    int num_cols = getNumSamples("snapshot");

    CAROM_VERIFY(0 < start_col <= num_cols);
    CAROM_VERIFY(start_col <= end_col && end_col <= num_cols);
    int num_cols_to_read = end_col - start_col + 1;

    char tmp[100];
    Matrix* snapshots = new Matrix(num_rows, num_cols_to_read, false);
    // This 0 index is the remainder from time interval concept.
    // This remains only for backward compatibility purpose.
    sprintf(tmp, "snapshot_matrix_000000");
    d_database->getDoubleArray(tmp,
                               &snapshots->item(0, 0),
                               num_rows*num_cols_to_read,
                               start_col - 1,
                               num_cols_to_read,
                               num_cols);
    return snapshots;
}
}
