/*******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT file for
 * details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
*******************************************************************************/

/*
 * Interface to the ScaLAPACK wrapper module implemented in Fortran.
 * Declarations in this file must be careful to respect the Fortran calling
 * conventions, and struct declarations must correspond to a Fortran type
 * declared with `bind(C)` and have the same order of fields.
 */

#ifndef included_scalapack_wrapper_h
#define included_scalapack_wrapper_h

#ifdef __cplusplus
extern "C" {
#else
#include <stdbool.h>
#endif

/* Ensure that this typedef matches the parameter definition in
 * 'scalapack_wrapper.f90' or strange things will happen.
 */
typedef double REAL_TYPE;

/**
 * This structure encapsulates all of the information required to describe a
 * distributed matrix in ScaLAPACK's block cyclic layout. Manipulating it
 * manually is probably not a good idea. Some fields that could be kept as
 * internal are provided here for the exceptions, notably: IO. If the matrix is
 * initialized using `initialize_matrix(&A, m, n, comm_size, 1, m / comm_size,
 * n)` you will obtain a matrix in a block row layout with pieces dealt out in
 * order of increasing rank; i.e. in rank 0 you will have A(1:m/comm_size, :),
 * up to rank N which will have the last rows. If m does not evenly divide
 * comm_size, the `mm` field will give the actual number of rows stored.
 * Of course the same thing could be applied for a block column layout.
 *
 * If accessing local data, recall that the storage order is _column major_
 *
 * Fields
 * ======
 * ctxt: The integer handle to the BLACS context on which this matrix lives.
 *       This one really is an opaque object! Exception: if ctxt is -1, this
 *       process does not own any part of the matrix.
 * m, n: *Global* dimensions of the matrix.
 * nprow, npcol: The number of process rows and columns in the BLACS process
 *               grid.
 * pi, pj: This process's coordinate (0-based) in the process grid.
 * mb, nb: The blocking factors in the row and column dimensions, respectively.
 *         32 x 32 or 64 x 64 is probably a good place to start for perf.
 * int mm, mn: The dimension of the local data stored. If distributed in true
 *             block cyclic format, it's probably best to leave manipulations of
 *             this array to ScaLAPACK functions.
 * mdata: Pointer to the local data space used to store this process's local
 *        portion of the matrix.
 */
struct SLPK_Matrix
{
    // Parameters describing global matrix layout.
    int ctxt;
    int m, n;
    int nprow, npcol;
    int pi, pj;
    int mb, nb;

    // Parameters describing the local storage layout
    int mm, mn;
    REAL_TYPE* mdata;
};

/**
 * @brief Initialize a brand new ScaLAPACK matrix.
 *
 * `A` is filled with the data for a new ScaLAPACK matrix. Its local data arrays
 * are allocated internally and a new BLACS context is created. No checking is
 * performed to make sure that `A` is not already initialized. This routine
 * must be called on each process.
 *
 * @pre m >= 0
 * @pre n >= 0
 * @pre nprow >= 1
 * @pre npcol >= 1
 * @pre nprow*npcol <= nprocs
 * @pre mb >= 1
 * @pre nb >= 1
 *
 * @param[in] m The row dimension of the global matrix
 * @param[in] n The column dimension of the global matrix
 * @param[in] nprow The number of process rows in the BLACS process grid
 * @param[in] npcol The number of process columns in the BLACS process grid
 * @param[in] mb The blocking factor in the row direction
 * @param[in] nb The blocking factor in the column direction
 * @param[out] A Filled with the data for a new ScaLAPACK matrix with given
 * parameters.
 */
void initialize_matrix(struct SLPK_Matrix* A, int m, int n, int nprow,
                       int npcol, int mb, int nb);

/**
 * @brief Free `A`'s local workspace.
 * @param[in] A Pointer to a `SLPK_Matrix`; `A`'s local array is freed, but its
 * BLACS context is not touched.
 */
void free_matrix_data(struct SLPK_Matrix* A);

/**
 * @brief Release the BLACS context where `A` is stored. Careful!
 * @param[in] A Pointer to `SLPK_Matrix`; `A->ctxt` is released. If other
 * matrices are described by this context, they will now be broken. If in doubt,
 * don't do this.
 */
void release_context(struct SLPK_Matrix* A);

/**
 * @brief Copy an m x n submatrix from `src` to `dst`.
 *
 * This routine performs the operation
 * `dst[dsti:dsti+m, dstj:dstj+n] = src[srci:srci+m, srcj:srcj+n]` using Matlab
 * notation. All indices are 1-based since this is just a wrapper around the
 * corresponding ScaLAPACK subroutine. `copy_matrix` must be called on every
 * process; if the process doesn't own a part of `dst` (`src`), make sure
 * that `dst->ctxt == -1` (`src->ctxt == -1`) to indicate this to the BLACS.
 *
 * @param[in] src Source matrix, dimension at least (srci+m)x(srcj+n)
 * @param[in] srci Row offset of the submatrix to copy from src. One based!
 * @param[in] srcj Column offset of the submatrix to copy from src. One based!
 * @param[in] m Row dimension of the submatrix to copy.
 * @param[in] n Column dimension of the submatrix to copy.
 * @param[in] dsti Row offset of the submatrix to fill in dst
 * @param[in] dstj Column offset of the submatrix to fill in dst
 * @param[inout] dst The indicated submatrix in `dst` is overwritten by the
 * indicated submatrix from `src`.
 */
void copy_matrix(struct SLPK_Matrix* dst, int dsti, int dstj,
                 const struct SLPK_Matrix* src, int srci, int srcj,
                 int m, int n);

/**
 * @brief Scatter an array from rank `proc` to a distributed matrix.
 *
 * This is a simple wrapper around `copy_matrix` to simplify code. It must be
 * called on all processes, since it calls `copy_matrix` internally. Performs
 * the operation `dst[dsti:dsti+m, dstj:dstj+n] = src` where `src` is
 * interpreted as an m x n _column major_ array.
 */
void scatter_block(struct SLPK_Matrix* dst, int dsti, int dstj,
                   const REAL_TYPE* src, int m, int n, int proc);

/**
 * @brief Gather a block from distributed matrix `src` to a block on rank `proc`
 *
 * This is a simple wrapper around `copy_matrix` to simplify code. It must be
 * called on all processes, since it calls `copy_matrix` internally. Performs
 * the operation `dst = src[srci:srci+m, srcj:srcj+n]` where `dst` is
 * interpreted as an m x n _column major_ array.
 */
void gather_block(REAL_TYPE* dst, const struct SLPK_Matrix* src, int srci,
                  int srcj, int m, int n, int proc);

/**
 * @brief Transpose a submatrix of `src` and store in a submatrix of `dst`.
 *
 * If `dst->ctxt != src->ctxt`, this routine must make a copy of `src` to a
 * matrix in `dst`'s context, then transpose that matrix. This is a limitation
 * of ScaLAPACK. Since the routine *may* copy a matrix, it needs to be called
 * from every process. Performs the operation
 * `dst[dsti:dsti+n, dstj:dstj+m] = src[srci:srci+m, srcj:srcj+n]'` using
 * Matlab notation. Indices are 1-based.
 */
void transpose_submatrix(struct SLPK_Matrix* dst, int dsti, int dstj,
                         const struct SLPK_Matrix* src, int srci, int srcj,
                         int m, int n);

/**
 * @brief Gather a block of distributed matrix `src` to a local array on rank
 * `proc`, and transpose it in between.
 *
 * This is a simple wrapper around `transpose_submatrix` to simplify code. Like
 * that routine, it needs to be called on all ranks. Performs the operation
 * `dst = src[srci:srci+m, srcj:srcj+n]' where `dst` is interpreted as an
 * n x m _column major_ array. This operation is also useful to change a
 * column major Fortran array to a row major C array.
 */
void gather_transposed_block(REAL_TYPE* dst, const struct SLPK_Matrix* src,
                             int srci, int srcj, int m, int n, int proc);
/**
 * @brief Construct a matrix from data stored on a single process.
 *
 * This function is provided as a convenience to use ScaLAPACK operations on
 * an array located on a single process. It needs to be called on all processes,
 * however. Used internally by the `copy` and `transpose` routines that take a
 * raw pointer to REAL_TYPE as an argument.
 *
 * @param[in] x Pointer to the data that will be contained in the new matrix.
 * Recall that this will be interpreted as a column major array.
 * @param[in] m Row dimension of x
 * @param[in] n Column dimension of x
 * @param[in] process MPI rank of the process that is calling.
 * Processes should all be calling these interface functions. If the calling
 * process is not the indicated rank, its A will have ctxt == -1 and not
 * participate in ScaLAPACK operations.
 * @param[out] A A new SLPK_Matrix that lives in a BLACS context encompassing
 * only the indicated process.
 */
void create_local_matrix(struct SLPK_Matrix* A, double* x, int m, int n,
                         int process);

/**
 * @brief Initialize dst with a new (uninitialized) m x n matrix, reusing the
 * indicated BLACS context.
 *
 * `dst` is initialized as an m x n block cyclic matrix in the BLACS
 * context given. If the matrix from which `ctxt` was obtained is still around,
 * don't call `release_context` until you're done with all matrices!
 *
 * @param[in] m Row dimension of dst
 * @param[in] n Column dimension of dst
 * @param[in] ctxt Handle to an already-initialized BLACS context to re-use;
 * obtain from an SLPK_Matrix initialized using another function from this
 * interface.
 * @param[in] mb Blocking factor for rows.
 * @param[in] nb Blocking factor for columns.
 * @param[out] dst Overwritten with the new matrix
 */
void make_similar_matrix(struct SLPK_Matrix* dst, int m, int n, int ctxt,
                         int mb, int nb);

/**
 * @brief Structure managing the call to the ScaLAPACK SVD.
 *
 * The SVDManager layout corresponds to the structure of a Fortran derived type;
 * don't rearrange it. The only matrix provided by the user is `A`; the others
 * will be allocated automatically. They will be reusable for multiple
 * factorizations of an `A` with the same layout if the structure is not
 * destroyed. Use `svd_init` to initialize the structure the first time. By
 * default, `dov` is set to 0 so the right singular vectors are not computed;
 * set to 1 to compute V.
 *
 * If reusing for a multiple factorization, you must set the `done` field to 0
 * to indicate that you really know what you're doing.
 */
struct SVDManager
{
    struct SLPK_Matrix* A, *U;
    double* S;
    struct SLPK_Matrix* V;
    // These are int not bool because of differences in C and Fortran handling
    // of logical variables.
    int dou, dov, done;
};

/**
 * @brief Initialize an SVDManager with the matrix to be factorized.
 */
void svd_init(struct SVDManager*, struct SLPK_Matrix* A);

/**
 * @brief Do the factorization. Allocates the U, V, S fields if they aren't
 * already initialized.
 */
void factorize(struct SVDManager*);

/**
 * @brief Release the wrapper's internal BLACS global context and finalize the
 * BLACS.
 */
void wrapper_finalize();

/**
 * @brief Prints the structure of the matrix A on each process to stdout
 */
void print_debug_info(struct SLPK_Matrix* A);

/**
 * @brief Structure managing the call to the ScaLAPACK linear solve (pdgesv).
 *
 * SLPK_Matrix A is the lhs.
 * SLPK_Matrix B is the rhs. The results (x in Ax=b) will be stored in B.
 */
struct LSManager
{
    struct SLPK_Matrix* A, *B;
    int* ipiv;
    int ipivSize;
};

/**
 * @brief Initialize the LSManager struct with SLPK_Matrices A and B.
 */
void ls_init(struct LSManager* mgr, struct SLPK_Matrix* A,
             struct SLPK_Matrix* B);

/**
 * @brief Perform the linear solve (pdgesv) with the data contained in the LSManager.
 */
void linear_solve(struct LSManager*);

/**
 * @brief Structure managing the call to the ScaLAPACK QR factorization.
 *
 * SLPK_Matrix A is the matrix to be factorized.
 */
struct QRManager
{
    struct SLPK_Matrix* A;
    double* tau;
    int* ipiv;
    int tauSize;
    int ipivSize;
};

/**
 * @brief Initialize the QRManager struct with SLPK_Matrix A.
 */
void qr_init(struct QRManager* mgr, struct SLPK_Matrix* A);

/**
 * @brief Perform the QR factorization (pdgelqf) with the data contained in the QRManager.
 *        The true Q is not returned, but Q is returned as a product of elementary
 *        reflectors. qrcompute(struct QRManager*) must be called afterwards to
 *        obtain Q.
 */
void qrfactorize(struct QRManager*);

/**
 * @brief Given the elementary reflectors stored in the QRManager, as a result of
 *        calling qrfactorize(struct QRManager*), obtain the action of q
 *        on SLPK_Matrix A.
 * @param[in] A The matrix to compute q's action on.
 * @param[in] S Which side to compute q's action on A (left or right side).
 * @param[in] T Whether to transpose Q before computing the action.
 */
void qaction(struct QRManager*, struct SLPK_Matrix* A, int S, int T);

/**
 * @brief Given the elementary reflectors stored in the QRManager, as a result of
 *        calling qrfactorize(struct QRManager*), compute the factorized Q.
 */
void qrcompute(struct QRManager*);

/**
 * @brief Perform the LQ factorization (pdgelqf) with the data contained in the QRManager.
 *        The true Q is not returned, but Q is returned as a product of elementary
 *        reflectors. lqcompute(struct QRManager*) must be called afterwards to
 *        obtain Q.
 */
void lqfactorize(struct QRManager*);

/**
 * @brief Given the elementary reflectors stored in the QRManager, as a result of
 *        calling lqfactorize(struct QRManager*), compute the factorized Q.
 */
void lqcompute(struct QRManager*);

#ifdef __cplusplus
}
#endif

#endif /* included_scalapack_wrapper_h */
