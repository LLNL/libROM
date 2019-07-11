/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

/*!
 * @file C_interface.h
 * 
 * Forward declaration of my C wrappers around ScaLAPACK routines.
 * 
 * The implementation in `C_interface.cpp` contains declarations or includes of
 * MKL headers as appropriate as well as translating these calls to Fortran
 * calling convention.
 */

#ifndef C_interface_h
#define C_interface_h

#ifdef __cplusplus
extern "C" {
#endif

int sl_init_wrapper(int nprow, int npcol);

void blacs_gridinfo_wrapper(int ctxt, int* nprow, int* npcol, int* pi, int* pj);

int blacs_gridmap_wrapper(const int* usermap, int ldumap, int nprow, int npcol);

void blacs_gridexit_wrapper(int ctxt);

int numroc_wrapper(int n, int nb, int iproc, int isrcproc, int nprocs);

double dlamch_wrapper(char what);

void pdgemr2d_wrapper(int m, int n, const double* src, int srci, int srcj,
                      const int* srcdesc, double* dstdata, int dsti, int dstj,
                      const int* dstdesc, int ctxt);

void pdgeadd_wrapper(char trans, int m, int n, double alpha, const double* a,
                     int ia, int ja, const int* desca, double beta, double* c,
                     int ic, int jc, const int* descc);

int dsyev_wrapper(char jobz, char uplo, int n, double* a, int ia, int ja,
                  const int* desca, double* w, double* z, int iz, int jz,
                  const int* descz, double* work, int lwork);

int pdsygvx_wrapper(int ibtype, char jobz, char range, char uplo, int n,
                    double* a, int ia, int ja, const int* desca, double* b,
                    int ib, int jb, const int* descb, double vl, double vu,
                    int il, int iu, double abstol, int* m, int* nz, double* w,
                    double orfac, double* z, int iz, int jz, const int* descz,
                    double* work, int lwork, int* iwork, int liwork,
                    int* ifail, int* iclustr, double* gap);

int pdgesvd_wrapper(char jobu, char jobvt, int m, int n, double* a, int ia,
                    int ja, const int* desca, double* s, double* u, int iu,
                    int ju, const int* descu, double* vt, int ivt, int jvt,
                    const int* descvt, double* work, int lwork);

void pdgemm_wrapper(char transa, char transb, int m, int n, int k, double alpha,
                    const double* a, int ia, int ja, const int* desca,
                    const double* b, int ib, int jb, const int* descb,
                    double beta, double* c, int ic, int jc, const int* descc);

void dgemm_wrapper(char transa, char transb, int m, int n, int k, double alpha,
                   const double* a, int lda, const double* b, int ldb,
                   double beta, double* c, int ldc);

void axpby_wrapper(int n, double a, const double* x, int incx, double b,
                   double* y, int incy);

double dnrm2_wrapper(int n, const double* x, int incx);

double ddot_wrapper(int n, const double* x, int incx, const double* y,
                    int incy);

void dscal_wrapper(int n, double a, double* x, int incx);

#ifdef __cplusplus
}
#endif

#endif /* C_interface_h */