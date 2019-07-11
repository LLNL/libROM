/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

#include "C_interface.h"
#include "fc_interface.h"

/* SL_INIT is not exported in mkl_scalapack.h for some reason.
 * If that changes this will need to be moved to the block for non-MKL
 * and in the MKL block, #define SL_INIT sl_init.
 */
#define SL_INIT FC_GLOBAL_(sl_init, SL_INIT)

extern "C"
{
    void SL_INIT(int *, const int *, const int *);
}

/*******************************************************************************
 * If the MKL is available it exports ScaLAPACK C declarations in
 * mkl_scalapack.h.
 ******************************************************************************/
#ifdef HAVE_MKL
#include "mkl_scalapack.h"
#include "mkl_blacs.h"
#include "mkl_pblas.h"
#include "mkl_lapack.h"
#include "mkl_blas.h"

#define BLACS_GET blacs_get
#define BLACS_GRIDINFO blacs_gridinfo
#define BLACS_GRIDMAP blacs_gridmap
#define BLACS_GRIDEXIT blacs_gridexit
#define NUMROC numroc
#define DLAMCH dlamch
#define PDGEMR2D pdgemr2d
#define PDGEADD pdgeadd
#define PDSYEV pdsyev
#define PDSYGVX pdsygvx
#define PDGESVD pdgesvd
#define PDGEMM pdgemm
#define DAXPBY daxpby
#define DNRM2 dnrm2
#define DDOT ddot
#define DSCAL dscal

/*******************************************************************************
 * If the MKL is not available we declare the ScaLAPACK API ourselves.
 ******************************************************************************/
#else /* no MKL */

#define BLACS_GET FC_GLOBAL_(blacs_get, BLACS_GET)
extern "C"
{
    void BLACS_GET(const int *, const int *what, int *val);
}

#define BLACS_GRIDINFO FC_GLOBAL_(blacs_gridinfo, BLACS_GRIDINFO)
extern "C"
{
    void BLACS_GRIDINFO(const int *, int *, int *, int *, int *);
}

#define BLACS_GRIDMAP FC_GLOBAL_(blacs_gridmap, BLACS_GRIDMAP)
extern "C"
{
    void BLACS_GRIDMAP(int *, const int *, const int *, const int *, const int *);
}

#define BLACS_GRIDEXIT FC_GLOBAL_(blacs_gridexit, BLACS_GRIDEXIT)
extern "C"
{
    void BLACS_GRIDEXIT(const int *);
}

#define NUMROC FC_GLOBAL(numroc, NUMROC)
extern "C"
{
    int NUMROC(const int *, const int *, const int *, const int *, const int *);
}

#define DLAMCH FC_GLOBAL(dlamch, DLAMCH)
extern "C"
{
    double DLAMCH(const char *cmach);
}

#define PDGEMR2D FC_GLOBAL(pdgemr2d, PDGEMR2D)
extern "C"
{
    void PDGEMR2D(const int *, const int *, const double *, const int *, const int *,
                  const int *, double *, const int *, const int *, const int *,
                  const int *);
}

#define PDGEADD FC_GLOBAL(pdgeadd, PDGEADD)
extern "C"
{
    void PDGEADD(const char *, const int *, const int *, const double *,
                 const double *, const int *, const int *, const int *,
                 const double *, double *, const int *, const int *, const int *);
}

#define PDSYEV FC_GLOBAL(pdsyev, PDSYEV)
extern "C"
{
    void PDSYEV(const char *, const char *, const int *, const double *, const int *,
                const int *, const int *, double *, double *, const int *,
                const int *, const int *, double *, const int *, int *);
}

#define PDSYGVX FC_GLOBAL(pdsygvx, PDSYGVX)
extern "C"
{
    void PDSYGVX(const int *, const char *, const char *, const char *, const int *,
                 double *, const int *, const int *, const int *, double *,
                 const int *, const int *, const int *, const double *,
                 const double *, const int *, const int *, const double *, int *,
                 int *, double *, const double *, double *, const int *, const int *,
                 const int *, double *, const int *, int *, const int *, int *, int *,
                 double *, int *);
}

#define PDGESVD FC_GLOBAL(pdgesvd, PDGESVD)
extern "C"
{
    void PDGESVD(const char *jobu, const char *jobvt, const int *m, const int *n,
                 const double *a, const int *ia, const int *ja, const int *desca,
                 double *s, double *u, const int *iu, const int *ju,
                 const int *descu, double *vt, const int *ivt, const int *jvt,
                 const int *descvt, double *work, const int *lwork, int *info);
}

#define PDGEMM FC_GLOBAL(pdgemm, PDGEMM)
extern "C"
{
    void PDGEMM(const char *, const char *, const int *, const int *, const int *,
                const double *, const double *, const int *, const int *, const int *,
                const double *, const int *, const int *, const int *, const double *,
                double *, const int *, const int *, const int *);
}

#define DAXPBY FC_GLOBAL(daxpby, DAXPBY)
extern "C"
{
    void DAXPBY(const int *, const double *, const double *, const int *,
               const double *, double *, const int *);
}

#define DDOT FC_GLOBAL(ddot, DDOT)
extern "C"
{
    double DDOT(const int*, const double*, const int*, const double*, const int*);
}

#define DSCAL FC_GLOBAL(dscal, DSCAL)
extern "C"
{
    void DSCAL(const int*, const double*, double*, const int*);
}

#define DNRM2 FC_GLOBAL(dnrm2, DNRM2)
extern "C"
{
    double DNRM2(const int*, const double*, const int*);
}

#endif

/*******************************************************************************
 * Implementation of the wrapper functions is here; macros that resolve to the
 * appropriate function name should have been added in the MKL/non-MKL block as
 * appropriate above.
 ******************************************************************************/

int sl_init_wrapper(int nprow, int npcol)
{
    int ctxt;
    SL_INIT(&ctxt, &nprow, &npcol);
    return ctxt;
}

void blacs_gridinfo_wrapper(int ctxt, int *nprow, int *npcol, int *pi, int *pj)
{
    BLACS_GRIDINFO(&ctxt, nprow, npcol, pi, pj);
}

int blacs_gridmap_wrapper(const int *usermap, int ldumap, int nprow, int npcol)
{
    int ctxt, dummy, what = 0;
    BLACS_GET(&dummy, &what, &ctxt);
    BLACS_GRIDMAP(&ctxt, usermap, &ldumap, &nprow, &npcol);
    return ctxt;
}

void blacs_gridexit_wrapper(int ctxt)
{
    BLACS_GRIDEXIT(&ctxt);
}

int numroc_wrapper(int n, int nb, int iproc, int isrcproc, int nprocs)
{
    return NUMROC(&n, &nb, &iproc, &isrcproc, &nprocs);
}

double dlamch_wrapper(char what)
{
    return DLAMCH(&what);
}

void pdgemr2d_wrapper(int m, int n, const double *src, int srci, int srcj,
                      const int *srcdesc, double *dst, int dsti, int dstj,
                      const int *dstdesc, int ctxt)
{
    PDGEMR2D(&m, &n, src, &srci, &srcj, srcdesc, dst, &dsti, &dstj, dstdesc,
             &ctxt);
}

void pdgeadd_wrapper(char trans, int m, int n, double alpha, const double *a,
                     int ia, int ja, const int *desca, double beta, double *c,
                     int ic, int jc, const int *descc)
{
    if (desca[1] != -1)
    {
        PDGEADD(&trans, &m, &n, &alpha, a, &ia, &ja, desca, &beta, c, &ic, &jc,
                descc);
    }
}

int pdsygvx_wrapper(int ibtype, char jobz, char range, char uplo, int n,
                    double *a, int ia, int ja, const int *desca, double *b,
                    int ib, int jb, const int *descb, double vl, double vu,
                    int il, int iu, double abstol, int *m, int *nz, double *w,
                    double orfac, double *z, int iz, int jz, const int *descz,
                    double *work, int lwork, int *iwork, int liwork,
                    int *ifail, int *iclustr, double *gap)
{
    if (desca[1] != -1)
    {
        int info;
        PDSYGVX(&ibtype, &jobz, &range, &uplo, &n, a, &ia, &ja, desca, b, &ib,
                &jb, descb, &vl, &vu, &il, &iu, &abstol, m, nz, w, &orfac, z,
                &iz, &jz, descz, work, &lwork, iwork, &liwork, ifail, iclustr,
                gap, &info);
        return info;
    }
    return 0;
}

int pdgesvd_wrapper(char jobu, char jobvt, int m, int n, double *a, int ia,
                    int ja, const int *desca, double *s, double *u, int iu,
                    int ju, const int *descu, double *vt, int ivt, int jvt,
                    const int *descvt, double *work, int lwork)
{
    if (desca[1] != -1)
    {
        int info;
        PDGESVD(&jobu, &jobvt, &m, &n, a, &ia, &ja, desca, s, u, &iu, &ju,
                descu, vt, &ivt, &jvt, descvt, work, &lwork, &info);
        return info;
    }
    return 0;
}

void pdgemm_wrapper(char transa, char transb, int m, int n, int k, double alpha,
                    const double *a, int ia, int ja, const int *desca,
                    const double *b, int ib, int jb, const int *descb,
                    double beta, double *c, int ic, int jc, const int *descc)
{
    if (desca[1] != -1)
    {
        PDGEMM(&transa, &transb, &m, &n, &k, &alpha, a, &ia, &ja, desca, b, &ib,
               &jb, descb, &beta, c, &ic, &jc, descc);
    }
}

void axpby_wrapper(int n, double alpha, const double* x, int incx, double beta,
                   double* y, int incy)
{
    DAXPBY(&n, &alpha, x, &incx, &beta, y, &incy);
}

double dnrm2_wrapper(int n, const double* x, int incx)
{
    return DNRM2(&n, x, &incx);
}

double ddot_wrapper(int n, const double* x, int incx, const double* y, int incy)
{
    return DDOT(&n, x, &incx, y, &incy);
}

void dscal_wrapper(int n, double a, double* x, int incx)
{
    DSCAL(&n, &a, x, &incx);
}
