/*******************************************************************************
 * 
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
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

    // Parameters desribing the local storage layout
    int mm, mn;
    REAL_TYPE* mdata;
};

/**
 * @brief Initialize a brand new ScaLAPACK matrix.
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
 *
 * @return In `A`, a matrix stored on the first `nprow * npcol` available
 * processes (in order of MPI rank). A new BLACS context is initialized to hold
 * the matrix.
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
 * @param[in] src Source matrix, dimension at least (srci+m)x(srcj+n)
 * @param[in] srci Row offset of the submatrix to copy from src. One based!
 * @param[in] srcj Column offset of the submatrix to copy from src. One based!
 * @param[in] m Row dimension of the submatrix to copy.
 * @param[in] n Column dimension of the submatrix to copy.
 * @param[in] dsti Row offset of the submatrix to fill in dst
 * @param[in] dstj Column offset of the submatrix to fill in dst
 *
 * @return `dst` with the m x n submatrix beginning at (dsti, dstj) filled with
 * the submatrix from `src`.
 */
// dst[dsti:dsti+m, dstj:dstj+n] = src[srci:srci+m, srcj:srcj+n]
void copy_matrix(struct SLPK_Matrix* dst, int dsti, int dstj,
                 const struct SLPK_Matrix* src, int srci, int srcj,
                 int m, int n);

/*
 * Scatter an array from rank `proc` to a distributed matrix. Must be called
 * on all processors.
 */
void scatter_block(struct SLPK_Matrix* dst, int dsti, int dstj,
                   const REAL_TYPE* src, int m, int n, int proc);

void gather_block(REAL_TYPE* dst, const struct SLPK_Matrix* src, int srci,
                  int srcj, int m, int n, int proc);

void transpose_submatrix(struct SLPK_Matrix* dst, int dsti, int dstj, 
                         const struct SLPK_Matrix* src, int srci, int srcj,
                         int m, int n);

void gather_transposed_block(REAL_TYPE* dst, const struct SLPK_Matrix* src,
                             int srci, int srcj, int m, int n, int proc);
/**
 * @brief Construct a matrix from data stored on a single process.
 *
 * @param[in] x Pointer to the data that will be contained in the new matrix.
 * Recall that this will be interpreted as a column major array.
 * 
 * @param[in] m Row dimension of x
 * @param[in] n Column dimension of x
 * @param[in] process MPI rank of the process that is calling.
 * Processes should all be calling these interface functions. If the calling
 * process is not the indicated rank, its A will have ctxt == -1 and not
 * participate in ScaLAPACK operations.
 *
 * @return `A` is initialized with its own BLACS context containing only the
 * indicated process, and associated to the data x. Treat `x` as owned by `A`
 * after calling this function; it is not copied.
 */
void create_local_matrix(struct SLPK_Matrix* A, double* x, int m, int n,
                         int process);

/**
 * @brief Initialize dst with a new (uninitialized) m x n matrix, reusing the
 * indicated BLACS context.
 *
 * @param[in] m Row dimension of dst
 * @param[in] n Column dimension of dst
 * @param[in] ctxt Handle to an already-initialized BLACS context to re-use;
 * obtain from an SLPK_Matrix initialized using another function from this
 * interface.
 * @param[in] mb Blocking factor for rows.
 * @param[in] nb Blocking factor for columns.
 *
 * @return `dst` is initialized as an m x n block cyclic matrix in the BLACS
 * context given. If the matrix from which `ctxt` was obtained is still around,
 * don't call `release_context` until you're done with all matrices!
 */
void make_similar_matrix(struct SLPK_Matrix* dst, int m, int n, int ctxt,
                         int mb, int nb);

struct SVDManager
{
    struct SLPK_Matrix* A, *U;
    double* S;
    struct SLPK_Matrix* V;
    // These are int not bool because of differences in C and Fortran handling
    // of logical variables.
    int dou, dov, done;
};

void svd_init(struct SVDManager*, struct SLPK_Matrix* A);
void factorize(struct SVDManager*);
void wrapper_finalize();

void print_debug_info(struct SLPK_Matrix* A);

#ifdef __cplusplus
}
#endif

#endif /* included_scalapack_wrapper_h */
