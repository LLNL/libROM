/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

#include "mpi.h"
#include "NNLS.h"
#include "scalapack_wrapper.h"

#include <set>

using namespace std;

namespace CAROM {

extern "C" {
    void blacs_get_(int*, int*, int*);
    void blacs_pinfo_(int*, int*);
    void blacs_gridinit_(int*, char*, int*, int*);
    void blacs_gridinfo_(int*, int*, int*, int*, int*);
    void blacs_gridexit_(int* icontxt);
    void blacs_freebuff_(int* icontxt, int* wait);
    void blacs_exit_(int* cont);
    void descinit_(int* DESC, int* M, int* N, int* MB, int* NB,
                   int* IRSRC, int* ICSRC, int* ICTXT,
                   int* LLD, int* INFO);
    void pdormqr_(char* SIDE, char* TRANS, int* M, int* N, int* K,
                  double* A, int* IA, int* JA, int* DESCA, double* TAU,
                  double* C, int* IC, int* JC, int* DESCC,
                  double* WORK, int* LWORK, int* INFO );
    void pdgeqrf_(int* M, int* N,
                  double* A, int* IA, int* JA, int* DESCA, double* TAU,
                  double* WORK, int* LWORK, int* INFO );
    void pdtrsm_(char* SIDE, char* UPLO, char* TRANS, char* DIAG,
                 int* M, int* N, double* ALPHA,
                 double* A, int* IA, int* JA, int* DESCA,
                 double* B, int* IB, int* JB, int* DESCB );
    void pdgemr2d_(int* m, int* n,
                   double* A, int* IA, int* JA, int* descA,
                   double* B, int* IB, int* JB, int* descB,
                   int* gcontext);
    void pdlacpy_(char* UPLO, int* M, int* N,
                  double* A, int* IA, int* JA, int* DESCA,
                  double* B, int* IB, int* JB, int* DESCB );
    int numroc_(int* N, int* NB, int* IPROC, int* ISRCPROC, int* NPROCS);
    void pdgemv_( char* TRANS, int * M, int * N, double * ALPHA,
                  double * A, int * IA, int * JA, int * DESCA,
                  double * X, int * IX, int * JX, int * DESCX, int * INCX,
                  double * BETA,
                  double * Y, int * IY, int * JY, int * DESCY, int * INCY );
}

NNLSSolver::NNLSSolver(double const_tol, int min_nnz, int max_nnz,
                       int verbosity,
                       double res_change_termination_tol,
                       double zero_tol, int n_outer, int n_inner)
    : const_tol_(const_tol), min_nnz_(min_nnz), max_nnz_(max_nnz),
      verbosity_(verbosity),
      res_change_termination_tol_(res_change_termination_tol),
      zero_tol_(zero_tol), n_outer_(n_outer), n_inner_(n_inner),
      n_proc_max_for_partial_matrix_(15),
      NNLS_qrres_on_(false),
      qr_residual_mode_(QRresidualMode::hybrid)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
    std::cout << "NNLSSolver init on rank " << d_rank << " out of "
              << d_num_procs << " processes" << std::endl;
}

NNLSSolver::~NNLSSolver()
{}

void NNLSSolver::set_verbosity(const int verbosity_in)
{
    verbosity_ = verbosity_in;
}

void NNLSSolver::set_qrresidual_mode(const QRresidualMode qr_residual_mode)
{
    qr_residual_mode_ = qr_residual_mode;
    if (qr_residual_mode_ == QRresidualMode::on) {
        NNLS_qrres_on_ = true;
    }
}

void NNLSSolver::normalize_constraints(Matrix& matTrans, Vector& rhs_lb,
                                       Vector& rhs_ub)
{
    // We scale everything so that rescaled half gap is the same for all constraints
    const unsigned int n = matTrans.numRows();
    const unsigned int m = matTrans.numColumns();

    CAROM_VERIFY(rhs_lb.dim() == m && rhs_ub.dim() == m);

    Vector rhs_avg = rhs_ub;
    rhs_avg += rhs_lb;
    rhs_avg *= 0.5;

    Vector rhs_halfgap = rhs_ub;
    rhs_halfgap -= rhs_lb;
    rhs_halfgap *= 0.5;

    Vector rhs_avg_glob = rhs_avg;
    Vector rhs_halfgap_glob = rhs_halfgap;
    Vector halfgap_target(m, true);
    halfgap_target = 1.0e3 * const_tol_;

    for (int i=0; i<m; ++i)
    {
        const double s = halfgap_target(i) / rhs_halfgap_glob(i);
        for (int j=0; j<n; ++j)
        {
            matTrans(j,i) *= s;
        }

        rhs_lb(i) = (rhs_avg(i) * s) - halfgap_target(i);
        rhs_ub(i) = (rhs_avg(i) * s) + halfgap_target(i);
    }
}

void NNLSSolver::solve_parallel_with_scalapack(const Matrix& matTrans,
        const Vector& rhs_lb, const Vector& rhs_ub, Vector& soln)
{
    CAROM_VERIFY(matTrans.distributed());

    int n = matTrans.numRows();
    int m = matTrans.numColumns();
    int n_tot = n;
    MPI_Allreduce(MPI_IN_PLACE, &n_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    CAROM_VERIFY(rhs_lb.dim() == m && rhs_lb.dim() == m && soln.dim() == n);
    if (max_nnz_ == 0)
        max_nnz_ = matTrans.numDistributedRows();

    // prepare right hand side
    Vector rhs_avg(rhs_ub);
    rhs_avg += rhs_lb;
    rhs_avg *= 0.5;

    Vector rhs_halfgap(rhs_ub);
    rhs_halfgap -= rhs_lb;
    rhs_halfgap *= 0.5;

    Vector rhs_avg_glob(rhs_avg);
    Vector rhs_halfgap_glob(rhs_halfgap);

    int izero = 0;
    int ione = 1;
    double fone = 1.0;
    int ictxt;
    char lside = 'L';
    char trans = 'T';
    char notrans = 'N';
    char layout = 'R';

    int n_proc = std::min(n_proc_max_for_partial_matrix_,d_num_procs);
    int nb = 3; // block column size

    // initialize blacs process grid
    blacs_get_(&izero, &izero, &ictxt);
    blacs_gridinit_(&ictxt, &layout, &ione, &n_proc);

    int n_dist_loc_max = 0; // maximum number of columns in distributed matrix
    if (d_rank < n_proc) {
        // m is the maximum number of columns in the global matrix
        n_dist_loc_max = ((m/nb + 1)/n_proc + 1)*nb;
    }

    std::vector<double> mu_max_array(d_num_procs);
    std::vector<unsigned int> proc_index;
    std::vector<unsigned int> nz_ind(m);
    Vector res_glob(m, false);
    Vector mu(n, true);
    Vector mu2(n, false);
    int n_nz_ind = 0;
    int n_glob = 0;
    int m_update;
    int min_nnz_cap = std::min(static_cast<int>(min_nnz_), std::min(m,n_tot));
    int info;
    Vector l2_res_hist(n_outer_, false);
    std::vector<unsigned int> stalled_indices;
    int stalledFlag = 0;
    int num_stalled = 0;
    int num_stalled_glob = 0;
    int nz_ind_zero = 0;
    int proc_zero = -1;

    Vector soln_nz_glob;
    Vector soln_nz_glob_up;

    // The following matrices are stored in column-major format as Vectors
    Vector mat_0_data(m * n_dist_loc_max, false);
    Vector mat_qr_data(m * n_dist_loc_max, false);

    int mat_qr_desc[9];
    Vector tau(n_dist_loc_max, false);
    Vector vec1;
    int vec1_desc[9];

    if (d_rank == 0)
    {
        const int n_glob_max = m;
        vec1.setSize(n_glob_max);
        soln_nz_glob.setSize(n_glob_max);
        soln_nz_glob_up.setSize(n_glob_max);
        proc_index.resize(n_glob_max);
    }

    // temporary work arrays
    int lwork;
    std::vector<double> work;
    int n_outer_iter = 0;
    int n_total_inner_iter = 0;
    int i_qr_start;
    int i_qr_start_f;
    int n_update;
    unsigned int iiter;
    // 0 = converged; 1 = maximum iterations reached; 2 = NNLS stalled (no change in residual for many iterations)
    int exit_flag = 1;

    // vector descriptor
    if (d_rank < n_proc) {
        descinit_(vec1_desc, &m, &ione, &m, &ione, &izero, &izero, &ictxt, &m, &info);
    }
    res_glob = rhs_avg_glob;
    Vector qt_rhs_glob = rhs_avg_glob;
    Vector qqt_rhs_glob = qt_rhs_glob;
    int qqt_rhs_desc[9];

    // compute threshold tolerance for the Lagrange multiplier mu
    //double mu_tol = (n > 0) ? 1.0e-15*arma::max(mat.t()*rhs_halfgap_glob) : 0.0; //1.0e-15*arma::norm(mat,1);
    double mu_tol = 0.0;

    {
        Vector tmp(n, true);
        //mat.transposeMult(rhs_halfgap_glob, tmp);
        matTrans.mult(rhs_halfgap_glob, tmp);
        double maxv = tmp(0);
        for (int i=1; i<n; ++i)
        {
            maxv = std::max(maxv, tmp(i));
        }

        mu_tol = 1.0e-15 * maxv;
    }

    MPI_Allreduce(MPI_IN_PLACE, &mu_tol, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    double mumax_glob;
    double rmax;

    for (unsigned int oiter = 0; oiter < n_outer_; ++oiter) {
        stalledFlag = 0;

        rmax = fabs(res_glob(0)) - rhs_halfgap_glob(0);
        for (int i=1; i<m; ++i)
            rmax = std::max(rmax, fabs(res_glob(i)) - rhs_halfgap_glob(i));

        // This norm of a non-distributed vector is identical on all ranks.
        l2_res_hist(oiter) = res_glob.norm();

        if (verbosity_ > 1 && d_rank == 0) {
            printf("%d %d %d %d %d %.15e %.15e\n", oiter, n_total_inner_iter,
                   m, n_tot, n_glob, rmax, l2_res_hist(oiter));
            fflush(stdout);
        }
        if (rmax <= const_tol_ && n_glob >= min_nnz_cap) {
            if (d_rank == 0 && verbosity_ > 1) {
                printf("target tolerance met\n");
                fflush(stdout);
            }
            exit_flag = 0;
            break;
        }

        if (n_glob >= max_nnz_) {
            if (d_rank == 0 && verbosity_ > 1) {
                printf("target nnz met\n");
                fflush(stdout);
            }
            exit_flag = 0;
            break;
        }

        if (n_glob >= m) {
            if (d_rank == 0 && verbosity_ > 1) {
                printf("system is square... exiting\n");
                fflush(stdout);
            }
            exit_flag = 3;
            break;
        }

        if (oiter > 101) {//we don't check for stall for the first 100 iterations
            //double mean_res_change = (arma::mean(l2_res_hist.rows(oiter-101,oiter-51))/arma::mean(l2_res_hist.rows(oiter-50,oiter))) - 1;
            double mean0 = 0.0;
            double mean1 = 0.0;
            for (int i=0; i<50; ++i)
            {
                mean0 += l2_res_hist(oiter - i);
                mean1 += l2_res_hist(oiter - 50 - i);
            }

            // Omitting division by 50, since a ratio is taken.
            double mean_res_change = (mean1 / mean0) - 1.0;
            if (std::abs(mean_res_change) < res_change_termination_tol_) {
                if (d_rank == 0 && verbosity_ > 1) {
                    printf("NNLS stall detected... exiting\n");
                    fflush(stdout);
                }
                exit_flag = 2;
                break;
            }
        }

        // find the next index
        //mu = mat.t()*res_glob;
        //mat.transposeMult(res_glob, mu);
        matTrans.mult(res_glob, mu);

        for (int i = 0; i < n_nz_ind; ++i) {
            mu(nz_ind[i]) = 0.0;
        }
        for (unsigned int i = 0; i < stalled_indices.size(); ++i) {
            mu(stalled_indices[i]) = 0.0;
        }

        //double mumax = (mu.n_elem > 0) ? arma::max(mu) : 0.0;
        double mumax = mu(0);
        for (int i=1; i<n; ++i)
            mumax = std::max(mumax, mu(i));

        mumax_glob = mumax;

        MPI_Allreduce(MPI_IN_PLACE, &mumax_glob, 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        if (mumax_glob < mu_tol) {
            num_stalled = stalled_indices.size();
            MPI_Allreduce(&num_stalled, &num_stalled_glob, 1, MPI_INT, MPI_SUM,
                          MPI_COMM_WORLD);
            if (num_stalled_glob > 0) {
                if (verbosity_ > 0 && d_rank == 0) {
                    std::cout << "Lagrange multiplier is below the minimum threshold: mumax = " <<
                              mumax_glob << ", mutol = " << mu_tol << "\n" <<
                              " Resetting stalled indices vector of size " << num_stalled_glob << "\n";
                }
                stalled_indices.resize(0);

                //mat.transposeMult(res_glob, mu);
                matTrans.mult(res_glob, mu);

                for (int i = 0; i < n_nz_ind; ++i) {
                    mu(nz_ind[i]) = 0.0;
                }

                //mumax = (mu.n_elem > 0) ? arma::max(mu) : 0.0;
                mumax = mu(0);
                for (int i=1; i<n; ++i)
                    mumax = std::max(mumax, mu(i));

                mumax_glob = mumax;
                MPI_Allreduce(MPI_IN_PLACE, &mumax_glob, 1, MPI_DOUBLE, MPI_MAX,
                              MPI_COMM_WORLD);
            }
        }

        int imax = 0;
        {
            double tmax = mu(0);
            for (int i=1; i<n; ++i)
            {
                if (mu(i) > tmax)
                {
                    tmax = mu(i);
                    imax = i;
                }
            }
        }

        //int imax = (mu.n_elem > 0) ? arma::index_max(mu) : -1;
        MPI_Gather(&mumax, 1, MPI_DOUBLE, mu_max_array.data(), 1, MPI_DOUBLE, 0,
                   MPI_COMM_WORLD);
        int imax_proc; // processor to which the next index belongs
        if (d_rank == 0) {
            //imax_proc = arma::index_max(mu_max_array);
            imax_proc = 0;
            double tmax = mu_max_array[0];
            for (int i=1; i<d_num_procs; ++i)
            {
                if (tmax < mu_max_array[i])
                {
                    tmax = mu_max_array[i];
                    imax_proc = i;
                }
            }

            proc_index[n_glob] = imax_proc;
        }
        MPI_Bcast(&imax_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (imax_proc == d_rank) {
            // record the local value of the next index
            nz_ind[n_nz_ind] = imax;
            ++n_nz_ind;
            // previous info
        }

        if (d_rank == imax_proc && verbosity_ > 2) {
            printf("found next index: %d %d %.15e\n", imax_proc, imax, mumax);
            fflush(stdout);
        }

        // send and recv the next column
        // recall that mat_0 and mat_qr are distributed using the cyclic format
        int proc_to_recv = (n_glob/nb) % n_proc; // proc to receive the new column
        if (proc_to_recv == imax_proc) {
            if (imax_proc == d_rank) {
                // local copy
                int n_orig = numroc_(&n_glob, &nb, &d_rank, &izero, &n_proc);
                for (int i=0; i<m; ++i)
                {
                    //mat_0_data(i + (n_orig*m)) = mat(i,imax);
                    mat_0_data(i + (n_orig*m)) = matTrans(imax,i);
                    mat_qr_data(i + (n_orig*m)) = mat_0_data(i + (n_orig*m));
                }
            }
        } else {
            // exchange data
            if (proc_to_recv == d_rank) {
                // recieve the matrix entry
                MPI_Status mpi_stat;
                int n_orig = numroc_(&n_glob, &nb, &d_rank, &izero, &n_proc);
                MPI_Recv(mat_0_data.getData() + m*n_orig, m, MPI_DOUBLE, imax_proc, 189,
                         MPI_COMM_WORLD, &mpi_stat);
                // copy the entry to the qr matrix
                //mat_qr.col(n_orig) = mat_0.col(n_orig);
                for (int i=0; i<m; ++i)
                    mat_qr_data(i + (n_orig*m)) = mat_0_data(i + (n_orig*m));
            }
            if (imax_proc == d_rank) {
                // send the partial matrix
                MPI_Send(matTrans.getData() + m*imax, m, MPI_DOUBLE, proc_to_recv, 189,
                         MPI_COMM_WORLD);
            }
        }

        i_qr_start = n_glob;
        ++n_glob; // increment the size of the global matrix

        if (d_rank == 0 && verbosity_ > 2) {
            printf("updated matrix with new index\n");
            fflush(stdout);
        }

        for (iiter = 0; iiter < n_inner_; ++iiter) {
            ++n_total_inner_iter;

            // initialize
            if (d_rank < n_proc) { // TODO: how can this condition not be true?
                descinit_(mat_qr_desc, &m, &n_glob, &m, &nb, &izero, &izero, &ictxt, &m, &info);
                CAROM_VERIFY(info == 0); // mat_qr descriptor initialization failed

                const bool incremental_update = true;
                i_qr_start_f = i_qr_start + 1;
                n_update = n_glob - i_qr_start;
                m_update = m - i_qr_start;
                if (incremental_update) {
                    // apply householder reflectors to compute Q^T new_cols
                    lwork = -1;
                    work.resize(10);

                    pdormqr_(&lside, &trans, &m, &n_update, &i_qr_start,
                             mat_qr_data.getData(), &ione, &ione, mat_qr_desc, tau.getData(),
                             mat_qr_data.getData(), &ione, &i_qr_start_f, mat_qr_desc,
                             work.data(), &lwork, &info);
                    CAROM_VERIFY(info == 0); // Q^T A update work calculation failed
                    lwork = static_cast<int>(work[0]);
                    work.resize(lwork);
                    pdormqr_(&lside, &trans, &m, &n_update, &i_qr_start,
                             mat_qr_data.getData(), &ione, &ione, mat_qr_desc, tau.getData(),
                             mat_qr_data.getData(), &ione, &i_qr_start_f, mat_qr_desc,
                             work.data(), &lwork, &info);
                    CAROM_VERIFY(info == 0); // Q^T A update failed
                    // compute QR factorization of the submatrix
                    lwork = -1;
                    work.resize(10);
                    pdgeqrf_(&m_update, &n_update,
                             mat_qr_data.getData(), &i_qr_start_f, &i_qr_start_f, mat_qr_desc, tau.getData(),
                             work.data(), &lwork, &info);
                    CAROM_VERIFY(info == 0); // QR update factorization work calculation
                    lwork = static_cast<int>(work[0]);
                    work.resize(lwork);
                    pdgeqrf_(&m_update, &n_update,
                             mat_qr_data.getData(), &i_qr_start_f, &i_qr_start_f, mat_qr_desc, tau.getData(),
                             work.data(), &lwork, &info);
                    CAROM_VERIFY(info == 0); // QR update factorization failed

                } else {
                    // copy everything to mat_qr then do full QR
                    int n_loc = numroc_(&n_glob, &nb, &d_rank, &izero, &n_proc);
                    for (int i=0; i<m; ++i)
                        for (int j=0; j<n_loc; ++j)
                            mat_qr_data(i + (j*m)) = mat_0_data(i + (j*m));

                    // compute qr factorization (first find the size of work and then perform qr)
                    lwork = -1;
                    work.resize(10);
                    pdgeqrf_(&m, &n_glob,
                             mat_qr_data.getData(), &ione, &ione, mat_qr_desc, tau.getData(),
                             work.data(), &lwork, &info);
                    CAROM_VERIFY(info == 0); // QR factorization work calculation
                    lwork = static_cast<int>(work[0]);
                    work.resize(lwork);
                    pdgeqrf_(&m, &n_glob,
                             mat_qr_data.getData(), &ione, &ione, mat_qr_desc, tau.getData(),
                             work.data(), &lwork, &info);
                    CAROM_VERIFY(info == 0); // QR factorization failed
                }

                if (d_rank == 0 && verbosity_ > 2) {
                    printf("updated QR %d\n", iiter);
                    fflush(stdout);
                }

                // apply householder reflectors to compute Q^T b
                if (incremental_update && iiter == 0) {
                    lwork = -1;
                    work.resize(10);
                    pdormqr_(&lside, &trans, &m_update, &ione, &ione,
                             mat_qr_data.getData(), &i_qr_start_f, &i_qr_start_f, mat_qr_desc, tau.getData(),
                             qt_rhs_glob.getData(), &i_qr_start_f, &ione, vec1_desc,
                             work.data(), &lwork, &info);
                    CAROM_VERIFY(info == 0); // H_last y work calculation failed
                    lwork = static_cast<int>(work[0]);
                    work.resize(lwork);
                    pdormqr_(&lside, &trans, &m_update, &ione, &ione,
                             mat_qr_data.getData(), &i_qr_start_f, &i_qr_start_f, mat_qr_desc, tau.getData(),
                             qt_rhs_glob.getData(), &i_qr_start_f, &ione, vec1_desc,
                             work.data(), &lwork, &info);
                    CAROM_VERIFY(info == 0); // H_last y failed
                } else {
                    // compute Q^T b from scratch
                    qt_rhs_glob = rhs_avg_glob;
                    lwork = -1;
                    work.resize(10);
                    pdormqr_(&lside, &trans, &m, &ione, &n_glob,
                             mat_qr_data.getData(), &ione, &ione, mat_qr_desc, tau.getData(),
                             qt_rhs_glob.getData(), &ione, &ione, vec1_desc,
                             work.data(), &lwork, &info);
                    CAROM_VERIFY(info == 0); // Q^T b work calculation failed
                    lwork = static_cast<int>(work[0]);
                    work.resize(lwork);
                    pdormqr_(&lside, &trans, &m, &ione, &n_glob,
                             mat_qr_data.getData(), &ione, &ione, mat_qr_desc, tau.getData(),
                             qt_rhs_glob.getData(), &ione, &ione, vec1_desc,
                             work.data(), &lwork, &info);
                    CAROM_VERIFY(info == 0); // Q^T b failed
                }

                if (d_rank == 0 && verbosity_ > 2) {
                    printf("updated rhs %d\n", iiter);
                    fflush(stdout);
                }

                // apply R^{-1}; first n_glob entries of vec1 are overwritten
                char upper = 'U';
                char nounit = 'N';
                vec1 = qt_rhs_glob;
                pdtrsm_(&lside, &upper, &notrans, &nounit,
                        &n_glob, &ione, &fone,
                        mat_qr_data.getData(), &ione, &ione, mat_qr_desc,
                        vec1.getData(), &ione, &ione, vec1_desc);

                if (d_rank == 0 && verbosity_ > 2) {
                    printf("solved triangular system %d\n", iiter);
                    fflush(stdout);
                }
            } // end of d_rank < n_proc

            // check if all entries are positive
            int pos_ibool = 0;
            if (d_rank == 0) {
                for (int i=0; i<n_glob; ++i)
                    soln_nz_glob_up(i) = vec1(i);

                if (soln_nz_glob_up.localMin(n_glob) > zero_tol_)
                {
                    pos_ibool = 1;
                    for (int i=0; i<n_glob; ++i)
                        soln_nz_glob(i) = soln_nz_glob_up(i);
                }
            }

            MPI_Bcast(&pos_ibool, 1, MPI_INT, 0, MPI_COMM_WORLD);

            if (pos_ibool == 1) {
                break;
            }

            if (d_rank == 0 && verbosity_ > 2) {
                printf("start pruning %d\n", iiter);
                for (int i = 0; i < n_glob; ++i) {
                    if (soln_nz_glob_up(i) <= zero_tol_) {
                        printf("%d %d %.6e\n", i, n_glob, soln_nz_glob_up(i));
                    }
                }
                fflush(stdout);
            }

            if (d_rank == 0) {
                if (soln_nz_glob_up(n_glob - 1) <= zero_tol_) {
                    stalledFlag = 1;
                    if (verbosity_ > 2) {
                        if (qr_residual_mode_ == QRresidualMode::hybrid) {
                            printf("Detected stall due to adding and removing same column. Switching to QR residual calculation method.\n");
                        } else {
                            printf("Detected stall due to adding and removing same column. Exiting now.\n");
                        }
                        fflush(stdout);
                    }
                }
            }
            MPI_Bcast(&stalledFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (stalledFlag == 1 && qr_residual_mode_ == QRresidualMode::hybrid) {
                NNLS_qrres_on_ = true;
                break;
            }

            if (d_rank == 0) { // we assume soln_nz_glob is stored on the root process
                double alpha = 1.0e300;
                // find maximum permissible step
                for (int i = 0; i < n_glob; ++i) {
                    if (soln_nz_glob_up(i) <= zero_tol_) {
                        alpha = std::min(alpha, soln_nz_glob(i)/(soln_nz_glob(i) - soln_nz_glob_up(i)));
                    }
                }
                // update solution
                double smin = 0.0;
                for (int i = 0; i < n_glob; ++i)
                {
                    soln_nz_glob(i) += alpha*(soln_nz_glob_up(i) - soln_nz_glob(i));
                    if (i == 0 || soln_nz_glob(i) < smin)
                        smin = soln_nz_glob(i);
                }

                while (smin > zero_tol_) {
                    // This means there was a rounding error, as we should have
                    // a zero element by definition. Recalcualte alpha based on
                    // the index that corresponds to the element that should be
                    // zero.

                    //int index_min = arma::index_min(soln_nz_glob.rows(0,n_glob-1));
                    int index_min = 0;
                    smin = soln_nz_glob(0);
                    for (int i = 1; i < n_glob; ++i)
                    {
                        if (soln_nz_glob(i) < smin)
                        {
                            smin = soln_nz_glob(i);
                            index_min = i;
                        }
                    }

                    // TODO: there are many max/min computations. Refactor by
                    // implementing functions for this in Vector. Also for
                    // min/max index.

                    alpha = soln_nz_glob(index_min)/(soln_nz_glob(index_min)
                                                     - soln_nz_glob_up(index_min));

                    // reupdate solution
                    //soln_nz_glob.rows(0,n_glob-1) += alpha*(soln_nz_glob_up.rows(0,n_glob-1) - soln_nz_glob.rows(0,n_glob-1));
                    for (int i = 0; i < n_glob; ++i)
                        soln_nz_glob(i) += alpha*(soln_nz_glob_up(i) - soln_nz_glob(i));
                }
            }

            // clean up zerod entry
            i_qr_start = n_glob+1;
            while (true) {
                // check if there is a zero entry
                int zero_ibool;
                if (d_rank == 0) {
                    //if (arma::min(soln_nz_glob.rows(0,n_glob-1)) < zero_tol_) {
                    if (soln_nz_glob.localMin(n_glob) < zero_tol_) {
                        zero_ibool = 1;
                    } else {
                        zero_ibool = 0;
                    }
                }
                MPI_Bcast(&zero_ibool, 1, MPI_INT, 0, MPI_COMM_WORLD);
                if (zero_ibool == 0) { // break if there is no more zero entry
                    break;
                }

                int ind_zero = -1; // index where the first zero is encountered
                proc_zero = -1;
                nz_ind_zero = 0;
                if (d_rank == 0) {
                    // identify global index of the zerod element
                    for (int i = 0; i < n_glob; ++i) {
                        if (soln_nz_glob(i) < zero_tol_) {
                            ind_zero = i;
                            break;
                        }
                    }
                    CAROM_VERIFY(ind_zero != -1); // something went wrong
                    // identify the process to which the zerod entry belongs
                    proc_zero = proc_index[ind_zero];
                    // identify the local index for nz_ind to which the zerod entry belongs
                    for (int i = 0; i < ind_zero; ++i) {
                        if (proc_index[i] == proc_zero) {
                            ++nz_ind_zero;
                        }
                    }
                } // end of d_rank == 0
                int zero_array[3] = {ind_zero, proc_zero, nz_ind_zero};
                MPI_Bcast(&zero_array, 3, MPI_INT, 0, MPI_COMM_WORLD);
                ind_zero = zero_array[0];
                proc_zero = zero_array[1];
                nz_ind_zero = zero_array[2];

                if (d_rank < n_proc) {
                    // copy mat_0.cols[ind_zero+1,n_glob) to mat_qr.cols[ind_zero,n_glob-1)
                    int ind_zero_f = ind_zero + 1; // FORTRAN index
                    int n_shift = n_glob - ind_zero_f;
                    int ind_zero_p1_f = ind_zero_f + 1;
                    pdgemr2d_(&m, &n_shift,
                              mat_0_data.getData(), &ione, &ind_zero_p1_f, mat_qr_desc,
                              mat_qr_data.getData(), &ione, &ind_zero_f, mat_qr_desc, &ictxt);
                    // copy mat_qr.cols[ind_zero,n_glob-1) to mat_qr.cols[ind_zero,n_glob-1)
                    char allmat = 'A';
                    pdlacpy_(&allmat, &m, &n_shift,
                             mat_qr_data.getData(), &ione, &ind_zero_f, mat_qr_desc,
                             mat_0_data.getData(), &ione, &ind_zero_f, mat_qr_desc);
                }
                blacs_freebuff_(&ictxt, &ione);

                // remove the zerod entry from the local matrix index
                if (d_rank == proc_zero) {
                    // remove the zerod entry
                    for (int i = nz_ind_zero; i < n_nz_ind-1; ++i) {
                        nz_ind[i] = nz_ind[i+1];
                    }
                    --n_nz_ind;
                }

                if (d_rank == 0) {
                    // shift soln_nz_glob and proc_index
                    for (int i = ind_zero; i < n_glob-1; ++i) {
                        soln_nz_glob(i) = soln_nz_glob(i+1);
                        proc_index[i] = proc_index[i+1];
                    }
                }
                i_qr_start = std::min(i_qr_start, ind_zero);
                --n_glob;
            } // end of pruning loop

            if (d_rank == 0 && verbosity_ > 2) {
                printf("finished pruning %d\n", iiter);
                fflush(stdout);
            }
        } // end of inner loop

        // check if we have stalled
        if (stalledFlag == 1) {
            iiter = 1;
            --n_glob;
            if (d_rank == imax_proc) {
                --n_nz_ind;
                num_stalled = stalled_indices.size();
                stalled_indices.resize(num_stalled + 1);
                stalled_indices[num_stalled] = imax;
                if (verbosity_ > 2) {
                    std::cout << "Adding index " << imax << " from processor " << imax_proc <<
                              " to stalled index list of size " << num_stalled << "\n";
                }
            }
        }

        // compute residual
        if (!NNLS_qrres_on_) {
            if (d_rank == 0) {
                res_glob = rhs_avg_glob;
            }
            if (d_rank < n_proc) {
                double fmone = -1.0;
                pdgemv_( &notrans, &m, &n_glob, &fmone,
                         mat_0_data.getData(), &ione, &ione, mat_qr_desc,
                         soln_nz_glob.getData(), &ione, &ione, vec1_desc, &ione, &fone,
                         res_glob.getData(), &ione, &ione, vec1_desc, &ione);
            }
            MPI_Bcast(res_glob.getData(), m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        } else {
            if (d_rank < n_proc) {
                // compute residual using res = b - Q*Q^T*b, where Q is from economical QR decomposition
                lwork = -1;
                work.resize(10);
                qqt_rhs_glob = 0.0;
                //qqt_rhs_glob(arma::span(0,n_glob - 1)) = qt_rhs_glob(arma::span(0,n_glob - 1));
                for (int i=0; i<n_glob; ++i)
                    qqt_rhs_glob(i) = qt_rhs_glob(i);

                descinit_(qqt_rhs_desc, &m, &ione, &m, &ione, &izero, &izero, &ictxt, &m,
                          &info);
                pdormqr_(&lside, &notrans, &m, &ione, &n_glob, mat_qr_data.getData(), &ione,
                         &ione,
                         mat_qr_desc, tau.getData(), qqt_rhs_glob.getData(), &ione, &ione, qqt_rhs_desc,
                         work.data(), &lwork, &info);
                CAROM_VERIFY(info == 0); // Q Q^T b work calculation failed.
                lwork = static_cast<int>(work[0]);
                work.resize(lwork);
                pdormqr_(&lside, &notrans, &m, &ione, &n_glob, mat_qr_data.getData(), &ione,
                         &ione,
                         mat_qr_desc, tau.getData(), qqt_rhs_glob.getData(), &ione, &ione, qqt_rhs_desc,
                         work.data(), &lwork, &info);
                CAROM_VERIFY(info == 0); // Q Q^T b calculation failed.
                //res_glob = rhs_avg_glob - qqt_rhs_glob;
                res_glob = rhs_avg_glob;
                res_glob -= qqt_rhs_glob;
            }
            MPI_Bcast(res_glob.getData(), m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

        if (d_rank == 0 && verbosity_ > 2) {
            printf("computed residual\n");
            fflush(stdout);
        }

        ++n_outer_iter;
    } // end of outer loop

    /* TODO (skipping this for now, as it's just verbose output)
    if (d_rank == 0 && verbosity_ > 1) {
      printf("finished NNLS loop\n");
      if (n_glob > 0) {
        Vector soln_sorted =  arma::sort(soln_nz_glob.rows(0, n_glob-1), "descend");
        for (unsigned int i = 0; i < soln_sorted.n_elem; ++i) {
          printf("%.6e ", soln_sorted(i));
        }
        printf("\n");
        fflush(stdout);
      }
    }
    */

    // send the solutions back to individual process
    std::vector<int> nnz_array;
    std::vector<int> nnz_cnt;
    Vector soln_nz_glob_psort;
    if (d_rank == 0) {
        nnz_array.assign(d_num_procs, 0);
        for (int i = 0; i < n_glob; ++i) {
            ++nnz_array[proc_index[i]];
        }
        nnz_cnt.assign(d_num_procs+1, 0);
        nnz_cnt[0] = 0;
        //nnz_cnt.rows(1,nnz_cnt.n_rows-1) = arma::cumsum(nnz_array);
        for (int i=1; i <= d_num_procs; ++i)
            nnz_cnt[i] = nnz_cnt[i - 1] + nnz_array[i - 1];

        soln_nz_glob_psort.setSize(n_glob);
        nnz_array.assign(d_num_procs, 0);
        for (int i = 0; i < n_glob; ++i) {
            const int proc = proc_index[i];
            soln_nz_glob_psort(nnz_cnt[proc] + nnz_array[proc]) = soln_nz_glob(i);
            ++nnz_array[proc];
        }
    }

    Vector soln_nz(std::max(n_nz_ind, 1), false);
    MPI_Scatterv(soln_nz_glob_psort.getData(), nnz_array.data(), nnz_cnt.data(),
                 MPI_DOUBLE, soln_nz.getData(), n_nz_ind, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // insert the solutions
    soln = 0.0;
    for (int i = 0; i < n_nz_ind; ++i) {
        soln(nz_ind[i]) = soln_nz(i);
    }

    blacs_freebuff_(&ictxt, &ione);

    if (d_rank < n_proc) {
        blacs_gridexit_(&ictxt);
    }

    if (d_rank == 0 && verbosity_ > 0) {
        printf("NNLS solver: m = %d, n = %d, outer_iter = %d, inner_iter = %d", m,
               n_tot, n_outer_iter, n_total_inner_iter);
        if (exit_flag == 0) {
            printf(": converged\n");
        } else {
            printf("\n");
            printf("Warning: NNLS unconverged. Stalled = %d\n", exit_flag == 2);
            printf("resErr = %.8e vs tol = %.8e; mumax = %.6e vs tol = %.6e\n", rmax,
                   const_tol_, mumax_glob, mu_tol);
        }
        fflush(stdout);
    }
}

} // namespace CAROM
