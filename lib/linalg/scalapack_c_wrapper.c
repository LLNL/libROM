/*******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT file for
 * details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 ******************************************************************************/

#include "scalapack_wrapper.h"
#include "mpi.h"

#include <stdlib.h>
#include <stdio.h>

void get_local_storage(struct SLPK_Matrix* A)
{
    A->mdata = malloc(sizeof(REAL_TYPE) * A->mm * A->mn);
}

void free_matrix_data(struct SLPK_Matrix* A)
{
    free(A->mdata);
    A->mdata = NULL;
}

void svd_init(struct SVDManager* mgr, struct SLPK_Matrix* A)
{
    mgr->A = A;
    mgr->U = NULL;
    mgr->S = NULL;
    mgr->V = NULL;
    mgr->dou = 1;
    mgr->dov = 0;
    mgr->done = 0;
}

void ls_init(struct LSManager* mgr, struct SLPK_Matrix* A,
             struct SLPK_Matrix* B)
{
    mgr->A = A;
    mgr->B = B;
    mgr->ipiv = NULL;
    mgr->ipivSize = 0;
}

void qr_init(struct QRManager* mgr, struct SLPK_Matrix* A)
{
    mgr->A = A;
    mgr->tau = NULL;
    mgr->ipiv = NULL;
    mgr->tauSize = 0;
    mgr->ipivSize = 0;
}

void factorize_prep(struct SVDManager* mgr)
{
    struct SLPK_Matrix* A = mgr->A;
    // Nomenclature here is from the SLUG.
    int SIZE = (A->m < A->n ? A->m : A->n);
    if (mgr->dou != 0 && mgr->U == NULL) {
        struct SLPK_Matrix* U = malloc(sizeof(struct SLPK_Matrix));
        make_similar_matrix(U, A->m, SIZE, A->ctxt, A->mb, A->nb);
        mgr->U = U;
    }

    if (mgr->dov != 0 && mgr->V == NULL) {
        struct SLPK_Matrix* V = malloc(sizeof(struct SLPK_Matrix));
        make_similar_matrix(V, SIZE, A->n, A->ctxt, A->mb, A->nb);
        mgr->V = V;
    }

    if (mgr->S == NULL) {
        mgr->S = malloc(sizeof(REAL_TYPE) * SIZE);
    }
}

void ls_prep(struct LSManager* mgr)
{
    if (mgr->ipiv == NULL && mgr->ipivSize > 0) {
        mgr->ipiv = malloc(sizeof(int) * mgr->ipivSize);
    }
}

void qrfactorize_prep(struct QRManager* mgr)
{
    if (mgr->ipiv == NULL && mgr->ipivSize > 0) {
        mgr->ipiv = malloc(sizeof(int) * mgr->ipivSize);
    }
    if (mgr->tau == NULL && mgr->tauSize > 0) {
        mgr->tau = malloc(sizeof(REAL_TYPE) * mgr->tauSize);
    }
}

typedef void (* void_func)(void*);

static void ordered_dowork(void_func f, void* arg)
{
    int mrank, commsize;
    int ierr = MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mrank);

    if (mrank > 0) {
        int dst;
        ierr = MPI_Recv((void*) &dst, 1, MPI_INT, mrank-1, 0, MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE);
    }

    f(arg);

    if (mrank+1 < commsize) {
        int src = 1;
        ierr = MPI_Send((const void*) &src, 1, MPI_INT, mrank+1, 0, MPI_COMM_WORLD);
    }
}

static void print_local(const struct SLPK_Matrix* A)
{
    int i, j;
    printf("Local data: (");
    // Remember it's column major.
    for (i = 0; i < A->mm; ++i) {
        for (j = 0; j < A->mn; ++j) {
            printf("%8.4f  ", A->mdata[j*A->mm + i]);
        }
        printf("\n             ");
    }
    printf(")\n");
}

static void show_info(void *ptr)
{
    struct SLPK_Matrix* A = (struct SLPK_Matrix*) ptr;
    int ierr, rank;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("Process rank %d\n", rank);
    printf("=================\n");
    printf("A->ctxt = %d\n", A->ctxt);
    if (A->ctxt == -1)
        return;
    printf("A->(m, n) = (%d, %d)\n", A->m, A->n);
    printf("A->(nprow, npcol) = (%d, %d)\n", A->nprow, A->npcol);
    printf("A->(pi, pj) = (%d, %d)\n", A->pi, A->pj);
    printf("A->(mb, nb) = (%d, %d)\n", A->mb, A->nb);
    printf("A->(mm, mn) = (%d, %d)\n", A->mm, A->mn);
    //print_local(A);
    printf("\n\n");
}

void print_debug_info(struct SLPK_Matrix* A)
{
    print_local(A);
    ordered_dowork(show_info, A);
}

void scatter_block(struct SLPK_Matrix* dst, int dsti, int dstj,
                   const REAL_TYPE* src, int m, int n, int proc)
{
    struct SLPK_Matrix srcmat;
    create_local_matrix(&srcmat, (REAL_TYPE*) src, m, n, proc);
    copy_matrix(dst, dsti, dstj, &srcmat, 1, 1, m, n);
    release_context(&srcmat);
}

void gather_block(REAL_TYPE* dst, const struct SLPK_Matrix* src, int srci,
                  int srcj, int m, int n, int proc)
{
    struct SLPK_Matrix dstmat;
    create_local_matrix(&dstmat, dst, m, n, proc);
    copy_matrix(&dstmat, 1, 1, src, srci, srcj, m, n);
    release_context(&dstmat);
}

void gather_transposed_block(REAL_TYPE* dst, const struct SLPK_Matrix* src,
                             int srci, int srcj, int m, int n, int proc)
{
    struct SLPK_Matrix dstmat;
    create_local_matrix(&dstmat, dst, n, m, proc);
    transpose_submatrix(&dstmat, 1, 1, src, srci, srcj, m, n);
    release_context(&dstmat);
}
