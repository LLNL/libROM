/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
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

#include <iomanip>
#include <fstream>

namespace CAROM {

BasisGenerator::BasisGenerator(
    Options options,
    bool incremental,
    const std::string& basis_file_name,
    Database::formats file_format) :
    d_dim(options.dim),
    d_incremental(incremental),
    d_basis_writer(0),
    d_basis_reader(0),
    d_write_snapshots(options.write_snapshots)
{
    CAROM_VERIFY(options.dim > 0);
    CAROM_VERIFY(options.max_num_samples > 0);
    CAROM_VERIFY(options.singular_value_tol >= 0);
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
        d_basis_writer = new BasisWriter(*this, basis_file_name, file_format);
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
                    options, basis_file_name));
        }
        else if (options.fast_update) {
            d_svd.reset(
                new IncrementalSVDFastUpdate(
                    options, basis_file_name));
        }
        else {
            d_svd.reset(
                new IncrementalSVDStandard(
                    options, basis_file_name));
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
            d_svd.reset(new RandomizedSVD(options));
        }
        else {
            d_svd.reset(new StaticSVD(options));
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
    bool add_without_increase)
{
    CAROM_VERIFY(u_in != 0);
    CAROM_VERIFY(d_svd->getNumSamples() < d_svd->getMaxNumSamples());

    // Check that u_in is not non-zero.
    Vector u_vec(u_in, getDim(), true);
    if (u_vec.norm() == 0.0) {
        printf("WARNING: BasisGenerator::takeSample skipped trivial sample.\n");
        return false;
    }

    return d_svd->takeSample(u_in, add_without_increase);
}

void
BasisGenerator::loadSampleRange(const std::string& base_file_name,
                                const std::string& kind,
                                int col_min,
                                int col_max,
                                Database::formats db_format)
{
    CAROM_ASSERT(!base_file_name.empty());
    CAROM_VERIFY(kind == "basis" || kind == "snapshot");

    if (d_basis_reader) delete d_basis_reader;

    d_basis_reader = new BasisReader(base_file_name, db_format, d_dim);
    std::unique_ptr<const Matrix> mat;
    std::unique_ptr<const Vector> singular_vals;

    if (kind == "basis") {
        mat = d_basis_reader->getSpatialBasis();
        singular_vals = d_basis_reader->getSingularValues();
    }
    else if (kind == "snapshot") {
        mat = d_basis_reader->getSnapshotMatrix();
    }

    int num_rows = mat->numRows();
    int num_cols = mat->numColumns();
    if (col_min < 0) col_min = 0;
    if (col_max > num_cols-1) col_max = num_cols-1;

    CAROM_VERIFY(col_max >= col_min);

    for (int j = col_min; j <= col_max; j++) {
        double* u_in = new double[num_rows];
        for (int i = 0; i < num_rows; i++) {
            if (kind == "basis") {
                u_in[i] = mat->item(i,j) * singular_vals->item(j);
            }
            else {
                u_in[i] = mat->item(i,j);
            }
        }
        d_svd->takeSample(u_in, false);
        delete[] u_in;
    }
}

void
BasisGenerator::loadSamples(const std::string& base_file_name,
                            const std::string& kind,
                            int cutoff,
                            Database::formats db_format)
{
    loadSampleRange(base_file_name, kind, 0, cutoff-1, db_format);
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
        std::shared_ptr<const Matrix> basis = getSpatialBasis();

        // Compute l = basis' * u
        Vector l(basis->numColumns(), false);
        basis->transposeMult(u_vec, l);

        // basisl = basis * l
        Vector basisl(basis->numRows(), true);
        basis->mult(l, basisl);

        // Compute u - basisl.
        std::unique_ptr<Vector> eta = u_vec.minus(basisl);

        // Compute l = basis' * rhs
        Vector rhs_vec(rhs_in, dim, true);
        basis->transposeMult(rhs_vec, l);

        // basisl = basis * l
        basis->mult(l, basisl);

        // Compute rhs - basisl.
        std::unique_ptr<Vector> eta_dot = rhs_vec.minus(basisl);

        // Compute the l-inf norm of eta + d_dt*eta_dot.
        double global_norm;
        double local_norm = 0.0;
        for (int i = 0; i < dim; ++i) {
            double val = fabs(eta->item(i) + d_dt*eta_dot->item(i));
            if (val > local_norm) {
                local_norm = val;
            }
        }
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

void
BasisGenerator::finalSummary(
    const double energyFractionThreshold,
    int & cutoff,
    const std::string & cutoffOutputPath,
    const int first_sv, const bool squareSV)
{
    const int rom_dim = getSpatialBasis()->numColumns();
    std::shared_ptr<const Vector> sing_vals = getSingularValues();

    CAROM_VERIFY(rom_dim <= sing_vals->dim());

    double sum = 0.0;
    for (int sv = first_sv; sv < sing_vals->dim(); ++sv) {
        const double s = (*sing_vals)(sv);
        sum += squareSV ? s * s : s;
    }

    int p = std::floor(-std::log10(energyFractionThreshold));
    std::vector<double> energy_fractions(p);

    for (int i = 0; i < p; ++i) {
        energy_fractions[i] = 1 - std::pow(10, -1 - i);
    }

    cutoff = first_sv;
    bool reached_cutoff = false;
    double partialSum = 0.0;
    int count = 0;

    std::ostream* output_stream;

    if (!cutoffOutputPath.empty()) {
        output_stream = new std::ofstream(cutoffOutputPath);
    } else {
        output_stream = &std::cout;
    }

    for (int sv = first_sv; sv < sing_vals->dim(); ++sv) {
        const double s = (*sing_vals)(sv);
        partialSum += squareSV ? s * s : s;
        for (int i = count; i < p; ++i)
        {
            if (partialSum / sum > 1.0 - std::pow(10, -1 - i))
            {
                *output_stream << "For energy fraction: 0.";
                for (int j = 0; j < i+1; ++j) *output_stream << "9";
                *output_stream << ", take first " << sv+1 << " of "
                               << sing_vals->dim() << " basis vectors" << std::endl;
                count += 1;
            }
            else
            {
                break;
            }
        }
        if (!reached_cutoff && partialSum / sum > 1.0 - energyFractionThreshold)
        {
            cutoff = sv+1;
            reached_cutoff = true;
        }
    }

    if (!reached_cutoff) cutoff = sing_vals->dim();
    *output_stream << std::fixed << std::setprecision(p+1);
    *output_stream << "For energy fraction: " << 1.0 - energyFractionThreshold <<
                   ", take first "
                   << cutoff << " of " << sing_vals->dim() << " basis vectors" << std::endl;

    if (!cutoffOutputPath.empty()) {
        static_cast<std::ofstream*>(output_stream)->close();
        delete output_stream;
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
