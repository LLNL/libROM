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
#include "utils/HDFDatabaseMPIO.h"
#include "Matrix.h"
#include "Vector.h"
#include "mpi.h"
#include "utils/mpi_utils.h"

namespace CAROM {

BasisReader::BasisReader(
    const std::string& base_file_name,
    Database::formats db_format,
    const int dim,
    MPI_Comm comm) :
    d_dim(dim),
    full_file_name(""),
    base_file_name_(base_file_name),
    d_format(db_format)
{
    CAROM_VERIFY(!base_file_name.empty());
    CAROM_VERIFY(comm == MPI_COMM_NULL || comm == MPI_COMM_WORLD);
    d_distributed = comm != MPI_COMM_NULL;

    int mpi_init;
    MPI_Initialized(&mpi_init);
    int rank;
    if (mpi_init && d_distributed) {
        MPI_Comm_rank(comm, &rank);
    }
    else {
        rank = 0;
    }

    full_file_name = base_file_name;

    // Enforce hdf data format.
    if (d_format == Database::formats::HDF5)
    {
        d_database = new HDFDatabase();
    }
    else if (d_format == Database::formats::HDF5_MPIO)
    {
        /*
            For MPIO case, local dimension needs to be specified.
            We allow 0 local dimension. (global dimension still needs to be positive)
        */
        std::vector<int> tmp;
        d_global_dim = get_global_offsets(d_dim, tmp, comm);
        CAROM_VERIFY(d_dim >= 0);
        CAROM_VERIFY(d_global_dim > 0);
        d_database = new HDFDatabaseMPIO();
    }
    else
        CAROM_ERROR("BasisWriter only supports HDF5/HDF5_MPIO data format!\n");

    d_database->open(full_file_name, "r", comm);
}

BasisReader::~BasisReader()
{
    d_database->close();
    delete d_database;
}

std::unique_ptr<Matrix>
BasisReader::getSpatialBasis()
{
    int num_rows = getDim("basis");
    int num_cols = getNumSamples("basis");

    Matrix* spatial_basis_vectors = new Matrix(num_rows, num_cols, d_distributed);

    d_database->getDoubleArray("spatial_basis",
                               &spatial_basis_vectors->item(0, 0),
                               num_rows*num_cols,
                               d_distributed);
    return std::unique_ptr<Matrix>(spatial_basis_vectors);
}

std::unique_ptr<Matrix>
BasisReader::getSpatialBasis(
    int n)
{
    return getSpatialBasis(1, n);
}

std::unique_ptr<Matrix>
BasisReader::getSpatialBasis(
    int start_col,
    int end_col)
{
    int num_rows = getDim("basis");
    int num_cols = getNumSamples("basis");

    char tmp[100];
    CAROM_VERIFY(0 < start_col && start_col <= num_cols);
    CAROM_VERIFY(start_col <= end_col && end_col <= num_cols);
    int num_cols_to_read = end_col - start_col + 1;

    Matrix* spatial_basis_vectors = new Matrix(num_rows, num_cols_to_read,
            d_distributed);
    sprintf(tmp, "spatial_basis");
    d_database->getDoubleArray(tmp,
                               &spatial_basis_vectors->item(0, 0),
                               num_rows*num_cols_to_read,
                               start_col - 1,
                               num_cols_to_read,
                               num_cols,
                               d_distributed);
    return std::unique_ptr<Matrix>(spatial_basis_vectors);
}

std::unique_ptr<Matrix>
BasisReader::getSpatialBasis(
    double ef)
{
    std::unique_ptr<Vector> sv = getSingularValues();
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

    return getSpatialBasis(num_used_singular_values);
}

std::unique_ptr<Matrix>
BasisReader::getTemporalBasis()
{
    int num_rows = getDim("temporal_basis");
    int num_cols = getNumSamples("temporal_basis");

    char tmp[100];
    Matrix* temporal_basis_vectors = new Matrix(num_rows, num_cols, false);
    sprintf(tmp, "temporal_basis");
    d_database->getDoubleArray(tmp,
                               &temporal_basis_vectors->item(0, 0),
                               num_rows*num_cols, d_distributed);
    return std::unique_ptr<Matrix>(temporal_basis_vectors);
}

std::unique_ptr<Matrix>
BasisReader::getTemporalBasis(
    int n)
{
    return getTemporalBasis(1, n);
}

std::unique_ptr<Matrix>
BasisReader::getTemporalBasis(
    int start_col,
    int end_col)
{
    int num_rows = getDim("temporal_basis");
    int num_cols = getNumSamples("temporal_basis");

    char tmp[100];
    CAROM_VERIFY(0 < start_col && start_col <= num_cols);
    CAROM_VERIFY(start_col <= end_col && end_col <= num_cols);
    int num_cols_to_read = end_col - start_col + 1;

    Matrix* temporal_basis_vectors = new Matrix(num_rows, num_cols_to_read, false);
    sprintf(tmp, "temporal_basis");
    d_database->getDoubleArray(tmp,
                               &temporal_basis_vectors->item(0, 0),
                               num_rows*num_cols_to_read,
                               start_col - 1,
                               num_cols_to_read,
                               num_cols,
                               d_distributed);
    return std::unique_ptr<Matrix>(temporal_basis_vectors);
}

std::unique_ptr<Matrix>
BasisReader::getTemporalBasis(
    double ef)
{
    std::unique_ptr<Vector> sv = getSingularValues();
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

    return getTemporalBasis(num_used_singular_values);
}

std::unique_ptr<Vector>
BasisReader::getSingularValues()
{
    char tmp[100];
    int size;
    sprintf(tmp, "singular_value_size");
    d_database->getInteger(tmp, size);

    Vector* singular_values = new Vector(size, false);
    sprintf(tmp, "singular_value");
    d_database->getDoubleArray(tmp,
                               &singular_values->item(0),
                               size,
                               d_distributed);
    return std::unique_ptr<Vector>(singular_values);
}

std::unique_ptr<Vector>
BasisReader::getSingularValues(
    double ef)
{
    std::unique_ptr<Vector> sv = getSingularValues();
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

    return std::unique_ptr<Vector>(truncated_singular_values);
}

int
BasisReader::getDim(
    const std::string kind)
{
    CAROM_ASSERT((kind == "basis") ||
                 (kind == "snapshot") ||
                 (kind == "temporal_basis"));

    char tmp[100];
    int num_rows;

    if (kind == "basis") sprintf(tmp, "spatial_basis_num_rows");
    else if (kind == "snapshot") sprintf(tmp, "snapshot_matrix_num_rows");
    else if (kind == "temporal_basis") sprintf(tmp,
                "temporal_basis_num_rows");

    d_database->getInteger(tmp, num_rows);
    /* only basis and snapshot are stored as distributed matrices */
    if ((kind != "temporal_basis") && (d_format == Database::formats::HDF5_MPIO))
    {
        /*
            for MPIO database, return specified local dimension.
            only checks the global dimension match.
        */
        CAROM_VERIFY(d_global_dim == num_rows);
        return d_dim;
    }
    else
        return num_rows;
}

int
BasisReader::getNumSamples(
    const std::string kind)
{
    CAROM_ASSERT((kind == "basis") ||
                 (kind == "snapshot") ||
                 (kind == "temporal_basis"));

    char tmp[100];
    int num_cols;

    if (kind == "basis") sprintf(tmp, "spatial_basis_num_cols");
    else if (kind == "snapshot") sprintf(tmp, "snapshot_matrix_num_cols");
    else if (kind == "temporal_basis") sprintf(tmp,
                "temporal_basis_num_cols");

    d_database->getInteger(tmp, num_cols);
    return num_cols;
}

std::unique_ptr<Matrix>
BasisReader::getSnapshotMatrix()
{
    int num_rows = getDim("snapshot");
    int num_cols = getNumSamples("snapshot");

    char tmp[100];
    Matrix* snapshots = new Matrix(num_rows, num_cols, d_distributed);
    sprintf(tmp, "snapshot_matrix");
    d_database->getDoubleArray(tmp,
                               &snapshots->item(0, 0),
                               num_rows*num_cols,
                               d_distributed);
    return std::unique_ptr<Matrix>(snapshots);
}

std::unique_ptr<Matrix>
BasisReader::getSnapshotMatrix(
    int n)
{
    return getSnapshotMatrix(1, n);
}

std::unique_ptr<Matrix>
BasisReader::getSnapshotMatrix(
    int start_col,
    int end_col)
{
    int num_rows = getDim("snapshot");
    int num_cols = getNumSamples("snapshot");

    CAROM_VERIFY(0 < start_col && start_col <= num_cols);
    CAROM_VERIFY(start_col <= end_col && end_col <= num_cols);
    int num_cols_to_read = end_col - start_col + 1;

    char tmp[100];
    Matrix* snapshots = new Matrix(num_rows, num_cols_to_read, d_distributed);
    sprintf(tmp, "snapshot_matrix");
    d_database->getDoubleArray(tmp,
                               &snapshots->item(0, 0),
                               num_rows*num_cols_to_read,
                               start_col - 1,
                               num_cols_to_read,
                               num_cols,
                               d_distributed);
    return std::unique_ptr<Matrix>(snapshots);
}
}
