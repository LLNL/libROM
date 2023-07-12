/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The abstract wrapper class for an abstract SVD algorithm and
//              sampler. This class provides interfaces to each so that an
//              application only needs to instantiate one concrete
//              implementation of this class to control all aspects of basis
//              vector generation.

#include "BasisGenerator.h"
#include "svd/StaticSVD.h"
#include "svd/RandomizedSVD.h"
#include "svd/IncrementalSVDStandard.h"
#include "svd/IncrementalSVDFastUpdate.h"
#include "svd/IncrementalSVDBrand.h"

namespace CAROM {

BasisGenerator::BasisGenerator(
    Options options,
    bool incremental,
    const std::string& basis_file_name,
    Database::formats file_format) :
    d_incremental(incremental),
    d_basis_writer(0),
    d_basis_reader(0),
    d_write_snapshots(options.write_snapshots)
{
    CAROM_VERIFY(options.dim > 0);
    CAROM_VERIFY(options.samples_per_time_interval > 0);
    CAROM_VERIFY(options.singular_value_tol >= 0);
    CAROM_VERIFY(options.max_time_intervals == -1
                 || options.max_time_intervals > 0);
    CAROM_VERIFY(options.max_basis_dimension > 0);
    if (incremental)
    {
        CAROM_VERIFY(options.linearity_tol > 0.0);
        CAROM_VERIFY(options.initial_dt > 0.0);
        CAROM_VERIFY(options.sampling_tol > 0.0);
        CAROM_VERIFY(options.max_time_between_samples > 0.0);
        CAROM_VERIFY(options.min_sampling_time_step_scale >= 0.0);
        CAROM_VERIFY(options.sampling_time_step_scale >= 0.0);
        CAROM_VERIFY(options.max_sampling_time_step_scale >= 0.0);
        CAROM_VERIFY(options.min_sampling_time_step_scale <=
                     options.max_sampling_time_step_scale);
    }

    if (!basis_file_name.empty()) {
        d_basis_writer = new BasisWriter(this, basis_file_name, file_format);
    }
    d_update_right_SV = options.update_right_SV;
    if (incremental)
    {
        d_tol = options.sampling_tol;
        d_max_time_between_samples = options.max_time_between_samples;
        d_min_sampling_time_step_scale = options.min_sampling_time_step_scale;
        d_sampling_time_step_scale = options.sampling_time_step_scale;
        d_max_sampling_time_step_scale = options.max_sampling_time_step_scale;
        d_dt = options.initial_dt;
        d_next_sample_time = 0.0;

        if (options.fast_update_brand) {
            d_svd.reset(
                new IncrementalSVDBrand(
                    options,
                    basis_file_name));
        }
        else if (options.fast_update) {
            d_svd.reset(
                new IncrementalSVDFastUpdate(
                    options,
                    basis_file_name));
        }
        else {
            d_svd.reset(
                new IncrementalSVDStandard(
                    options,
                    basis_file_name));
        }

        // Get the number of processors.
        int mpi_init;
        MPI_Initialized(&mpi_init);
        if (mpi_init) {
            MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
        }
        else {
            d_num_procs = 1;
        }
    }
    else
    {
        if (options.randomized) {
            d_svd.reset(
                new RandomizedSVD(
                    options));
        }
        else {
            d_svd.reset(
                new StaticSVD(
                    options));
        }
    }
}

bool
BasisGenerator::isNextSample(
    double time)
{
    CAROM_VERIFY(time >= 0.0);
    if (d_incremental)
    {
        if(d_update_right_SV)
            return true;
        else
            return time >= d_next_sample_time;
    }

    return true;
}

bool
BasisGenerator::takeSample(
    double* u_in,
    double time,
    double dt,
    bool add_without_increase)
{
    CAROM_VERIFY(u_in != 0);
    CAROM_VERIFY(time >= 0);

    // Check that u_in is not non-zero.
    Vector u_vec(u_in, getDim(), true);
    if (u_vec.norm() == 0.0) {
        printf("WARNING: BasisGenerator::takeSample skipped trivial sample at time %.4E\n",
               time);
        return false;
    }

    if (getNumBasisTimeIntervals() > 0 &&
            d_svd->isNewTimeInterval()) {
        resetDt(dt);
        if (d_basis_writer) {
            if (d_write_snapshots) {
                writeSnapshot();
            }
            else {
                d_basis_writer->writeBasis("basis");
            }
        }
    }

    return d_svd->takeSample(u_in, time, add_without_increase);
}

void
BasisGenerator::loadSamples(const std::string& base_file_name,
                            const std::string& kind,
                            int cut_off,
                            Database::formats db_format)
{
    CAROM_ASSERT(!base_file_name.empty());
    CAROM_VERIFY(kind == "basis" || kind == "snapshot");

    if (d_basis_reader) delete d_basis_reader;

    d_basis_reader = new BasisReader(base_file_name, db_format);
    double time = 0.0;
    const Matrix* mat;
    const Vector* singular_vals;

    if (kind == "basis") {
        mat = d_basis_reader->getSpatialBasis(time);
        singular_vals = d_basis_reader->getSingularValues(time);
    }
    else if (kind == "snapshot") {
        mat = d_basis_reader->getSnapshotMatrix(time);
    }

    int num_rows = mat->numRows();
    int num_cols = mat->numColumns();
    int max_cols = num_cols;
    if (cut_off < num_cols) max_cols = cut_off;

    for (int j = 0; j < max_cols; j++) {
        double* u_in = new double[num_rows];
        for (int i = 0; i < num_rows; i++) {
            if (kind == "basis") {
                u_in[i] = mat->item(i,j) * singular_vals->item(j);
            }
            else {
                u_in[i] = mat->item(i,j);
            }
        }
        d_svd->takeSample(u_in, time, false);
        delete[] u_in;
    }
}

double
BasisGenerator::computeNextSampleTime(
    double* u_in,
    double* rhs_in,
    double time)
{
    CAROM_VERIFY(u_in != 0);
    CAROM_VERIFY(rhs_in != 0);
    CAROM_VERIFY(time >= 0);
    if (d_incremental)
    {

        // Check that u_in is not non-zero.
        int dim = getDim();
        Vector u_vec(u_in, dim, true);
        if (u_vec.norm() == 0.0) {
            return d_next_sample_time;
        }

        // Get the current basis vectors.
        const Matrix* basis = getSpatialBasis();

        // Compute l = basis' * u
        Vector* l = basis->transposeMult(u_vec);

        // basisl = basis * l
        Vector* basisl = basis->mult(l);

        // Compute u - basisl.
        Vector* eta = u_vec.minus(basisl);

        delete l;
        delete basisl;

        // Compute l = basis' * rhs
        Vector rhs_vec(rhs_in, dim, true);
        l = basis->transposeMult(rhs_vec);

        // basisl = basis * l
        basisl = basis->mult(l);

        // Compute rhs - basisl.
        Vector* eta_dot = rhs_vec.minus(basisl);

        delete l;
        delete basisl;

        // Compute the l-inf norm of eta + d_dt*eta_dot.
        double global_norm;
        double local_norm = 0.0;
        for (int i = 0; i < dim; ++i) {
            double val = fabs(eta->item(i) + d_dt*eta_dot->item(i));
            if (val > local_norm) {
                local_norm = val;
            }
        }
        delete eta;
        delete eta_dot;
        if (d_num_procs == 1) {
            global_norm = local_norm;
        }
        else {
            MPI_Allreduce(&local_norm,
                          &global_norm,
                          1,
                          MPI_DOUBLE,
                          MPI_MAX,
                          MPI_COMM_WORLD);
        }

        // Compute dt from this norm.
        double tmp = d_sampling_time_step_scale*sqrt(d_tol/global_norm);
        if (tmp < d_min_sampling_time_step_scale) {
            d_dt *= d_min_sampling_time_step_scale;
        }
        else if (tmp > d_max_sampling_time_step_scale) {
            d_dt *= d_max_sampling_time_step_scale;
        }
        else {
            d_dt *= tmp;
        }
        if (d_dt < 0) {
            d_dt = 0.0;
        }
        else if (d_dt > d_max_time_between_samples) {
            d_dt = d_max_time_between_samples;
        }

        // Return next sample time.
        d_next_sample_time = time + d_dt;
        return d_next_sample_time;
    }

    return time;
}

void
BasisGenerator::resetDt(
    double new_dt)
{
    if (d_incremental)
    {
        d_dt = new_dt;
    }
}

BasisGenerator::~BasisGenerator()
{
    if (d_basis_writer) {
        delete d_basis_writer;
    }
    if (d_basis_reader) {
        delete d_basis_reader;
    }
}

}
