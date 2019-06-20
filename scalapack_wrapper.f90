!*******************************************************************************
!
! Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
! and other libROM project developers. See the top-level COPYRIGHT file for
! details.
!
! SPDX-License-Identifier: (Apache-2.0 OR MIT)
!
!*******************************************************************************

! This module is a C-callable wrapper for ScaLAPACK implementing an API to
! perform a distributed static SVD. It is interoperable with C code without
! name-mangling detection; the associated header file 'scalapack_wrapper.h'
! declares the C interface to this functionality with documentation.

module scalapack_wrapper
    use, intrinsic :: ISO_C_BINDING
    implicit none

!*******************************************************************************
! Begin declarations
!*******************************************************************************

    integer, parameter :: REAL_KIND = C_DOUBLE

    ! This integer variable will hold a handle to a BLACS context with all
    ! processes.
    integer :: GLOBAL_CTXT = -1

    ! This interface block contains the declarations for all ScaLAPACK routines
    ! used. I have added intent specifications on the arguments as appropriate,
    ! since this doesn't change the calling convention and conveys extra
    ! information both to the developer reading this code and to the compiler.
    interface
        subroutine sl_init( ictxt, nprow, npcol )
            integer, intent(out) :: ictxt
            integer, intent(in) :: nprow, npcol
        end subroutine

        subroutine blacs_gridinfo(ictxt, nprow, npcol, myprow, mypcol)
            integer, intent(in) :: ictxt
            integer, intent(out) :: nprow, npcol, myprow, mypcol
        end subroutine

        subroutine blacs_gridmap(ictxt, usermap, ldumap, nprow, npcol)
            integer, intent(inout) :: ictxt
            integer, intent(in) :: ldumap, nprow, npcol
            integer, intent(in) :: usermap(ldumap, npcol)
        end subroutine

        subroutine blacs_gridexit(ictxt)
            integer, intent(in) :: ictxt
        end subroutine

        subroutine blacs_exit(cont)
            integer, intent(in) :: cont
        end subroutine

        subroutine descinit(desc, m, n, mb, nb, irsrc, icsrc, ictxt, lld, info)
            integer, intent(in) :: m, n, mb, nb, irsrc, icsrc, ictxt, lld
            integer, intent(out) :: desc(*), info
        end subroutine

        subroutine blacs_get(icontxt, what, val)
            integer, intent(in) :: icontxt
            integer, intent(in) :: what
            integer, intent(out) :: val
        end subroutine

        subroutine pdgemr2d( m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt)
            integer, intent(in) :: ia, ib, ictxt, ja, jb, m, n
            integer, intent(in) :: desca(*), descb(*)
            double precision, intent(in) :: a(*)
            double precision, intent(out) :: b(*)
        end subroutine

        subroutine pdgesvd(jobu, jobvt, m, n, a, ia, ja, desca, s, u, &
                         & iu, ju, descu, vt, ivt, jvt, descvt, work, lwork, &
                         & info)
            character jobu, jobvt
            integer ia, info, iu, ivt, ja, ju, jvt, lwork, m, n
            integer desca(*), descu(*), descvt(*)
            double precision a(*), u(*), vt(*), work(*)
            double precision s(*)
        end subroutine

        subroutine pdgeadd(trans, m, n, alpha, A, ia, ja, desca, beta, C, &
                         & ic, jc, descc)
            character trans
            integer m, n, ia, ja, desca(*), ic, jc, descc(*)
            double precision alpha, A(*), beta, C(*)
        end subroutine

        function numroc(n, nb, iproc, isrcproc, nprocs)
            integer :: numroc
            integer :: n, nb, iproc, isrcproc, nprocs
        end function
    end interface
    
    ! An interoperable struct encapsulating a ScaLAPACK matrix - contains all of
    ! the information needed to use ScaLAPACK operations on the matrix.
    type, bind(C) :: SLPK_Matrix
        integer(C_INT) :: ctxt ! The integer handle to the ScaLAPACK context
                               ! that handles the matrix.
        integer(C_INT) :: m, n ! The dimensions of the matrix
        integer(C_INT) :: nprow, npcol ! The number of rows and columns in the
                                       ! ScaLAPACK processor grid
        integer(C_INT) :: pi, pj ! This process's coordinates (0-based) in the
                                 ! processor grid.
        integer(C_INT) :: mb, nb ! The blocking factor in the row- and column
                                 ! directions. 32 x 32 is probably a good choice
                                 ! for performance.
        integer(C_INT) :: mm, mn ! The size of the local data stored on this
                                 ! processor
        type(C_PTR) :: mdata ! Pointer obtained from C to an array of data of
                             ! type real(REAL_KIND) and size mm x mn
    end type SLPK_Matrix

    ! An interoperable struct that implements the computation of a distributed
    ! SVD. The user simply calls `svd_init` with the matrix A to factorize, and
    ! may set the values of the flags `dou` and `dov` as desired, then calls
    ! `factorize` to compute the SVD.
    type, bind(C) :: SVDManager
        type(C_PTR) :: A, U, S, V ! These are pointers to SLPK_Matrix in C.
        integer(C_INT) :: dou, dov, done ! Integer to avoid logical interop
    end type SVDManager

    ! Declarations of functionality implemented in the associated C code.
    interface
        subroutine get_local_storage(A) bind(C)
            import SLPK_Matrix
            type(SLPK_Matrix), intent(inout) :: A
        end subroutine

        subroutine free_matrix_data(A) bind(C)
            import SLPK_Matrix
            type(SLPK_Matrix), intent(inout) :: A
        end subroutine
        
        subroutine svd_init(mgr, A) bind(C)
            import SVDManager, SLPK_Matrix
            type(SVDManager), intent(out) :: mgr
            type(SLPK_Matrix), intent(inout) :: A
        end subroutine
    end interface

!*******************************************************************************
! End declarations; begin implementation.
!*******************************************************************************
    
contains

subroutine check_init()
    use MPI
    integer :: sz, ierr
    logical :: initialized

    call MPI_Initialized(initialized, ierr)
    if (.not. initialized) then
        call MPI_Init(ierr)
    endif

    call MPI_Comm_size(MPI_COMM_WORLD, sz, ierr)
    if (GLOBAL_CTXT .eq. -1) then
        call sl_init(GLOBAL_CTXT, sz, 1)
    endif
end subroutine

!*******************************************************************************
! Subroutine: initialize_matrix(A, m, n, nprow, npcol, mb, nb) bind(C)
!
! Synopsis: Initialize a SLPK_Matrix structure with a processor grid and
! allocated storage.
!
! Detail: `A` is treated as a new, uninitialized object and contents are
! replaced blindly. This initializer performs a few tasks: it initializes MPI if
! a call to MPI_Init has not been made yet; it checks that the number of
! processors available is enough to initialize the requested grid. Then it fills
! appropriate fields in `A` and calls `sl_init` to initialize a BLACS process
! grid, then queries the grid to get this process's coordinates. With that
! information it computes how much local storage is needed and acquires that
! storage, and precomputes the coordinates, global and local, of the blocks that
! this processor will store.
!
! Arguments
! =========
!   - A: A (pointer to when calling from C) SLPK_Matrix. On exit, `A` is
!     initialized with all of the information required to proceed with a
!     computation.
!   - m, n: Number of rows and columns of the global matrix. For the ROM
!     application, the data dimension and the number of shapshot vectors
!     respectively.
!   - nprow, npcol: Number of rows and columns in the processor grid;
!     MPI_Comm_size must give a number _greater than or equal to_ nprow*npcol.
!     Excess processes are allowed.
!   - mb, nb: Blocking factor for rows and columns, respectively. If nb == n and
!     m == mb * nprow, recover the block row storage scheme.
!*******************************************************************************
subroutine initialize_matrix( A, m, n, nprow, npcol, mb, nb ) bind(C)
    use mpi, only: MPI_Initialized, MPI_Init, MPI_Comm_size, MPI_COMM_WORLD
    use ISO_FORTRAN_ENV, only: error_unit

    type(SLPK_Matrix), intent(out) :: A
    integer(C_INT), value :: m, n, nprow, npcol, mb, nb
    logical :: initialized
    integer :: ierr, csize

    call check_init()

    call MPI_Comm_size(MPI_COMM_WORLD, csize, ierr)
    if (csize .lt. nprow*npcol) then
        write(error_unit, *) "Error: MPI_COMM_WORLD size < nprows * npcols"
        call exit(1)
    end if

    ! Create the process grid, and get my coordinates in it.
    A%m = m
    A%n = n
    A%mb = mb
    A%nb = nb
    call sl_init( A%ctxt, nprow, npcol )
    call blacs_gridinfo(A%ctxt, A%nprow, A%npcol, A%pi, A%pj)

    ! Compute the size of the local storage that I need, and allocate it (call to
    ! the C function to provision storage).
    A%mm = numroc(m, mb, A%pi, 0, nprow)
    A%mn = numroc(n, nb, A%pj, 0, npcol)
    call get_local_storage(A)
end subroutine

!*******************************************************************************
! End subroutine initialize_matrix
!*******************************************************************************




!*******************************************************************************
! Subroutine: make_similar_matrix(dst, m, n, ctxt, mb, nb) bind(C)
!
! Synopsis: Creates an m x n matrix in `dst` with the given BLACS parameters.
!
! Detail: This is provided as a tool to easily re-use the same BLACS context as
! has already been obtained. If you have initialized multiple matrices in the
! same context in this way, remember not to `release_context` on any of them
! unless you're done with all of them.
!
! Arguments
! =========
!   - dst: An SLPK_Matrix (pointer to, from C)
!   - m, n: Dimensions of the matrix to create.
!   - ctxt: BLACS context to re-use.
!   - mb, nb: Blocking factors for rows and columns
!*******************************************************************************
subroutine make_similar_matrix(dst, m, n, ctxt, mb, nb) bind(C)
    type(SLPK_Matrix), intent(out) :: dst
    integer(C_INT), value :: m, n, ctxt, mb, nb

    dst%m = m; dst%n = n; dst%mb = mb; dst%nb = nb; dst%ctxt = ctxt
    call blacs_gridinfo(ctxt, dst%nprow, dst%npcol, dst%pi, dst%pj)
    dst%mm = numroc(m, mb, dst%pi, 0, dst%nprow)
    dst%mn = numroc(n, nb, dst%pj, 0, dst%npcol)
    call get_local_storage(dst)
end subroutine

!*******************************************************************************
! End subroutine make_similar_matrix
!*******************************************************************************




!*******************************************************************************
! Subroutine: wrapper_finalize() bind(C)
!
! Synopsis: Release the global BLACS context initialized by this wrapper and
! finalize the BLACS, freeing its resources. You will not be able to use
! ScaLAPACK again after this function is called.
!*******************************************************************************
subroutine wrapper_finalize() bind(C)
    integer :: i
    if (GLOBAL_CTXT .ne. -1) then
        call blacs_gridexit(GLOBAL_CTXT)
    endif
    
    call blacs_exit(1)
end subroutine

!*******************************************************************************
! End subroutine wrapper_finalize
!*******************************************************************************




!*******************************************************************************
! Subroutine: release_context(A) bind(C)
! 
! Synopsis: Free the resources used by A's associated BLACS context.
!
! Detail: *WARNING* Make sure the context associated with A is not still in use
! for computations on other matrices, or that this subroutine is not called on
! multiple matrices sharing a context! I don't know what will happen if you do.
! Use `free_matrix_data` to release the local storage of a matrix but not touch
! its context.
!*******************************************************************************
subroutine release_context(A) bind(C)
    use MPI
    type(SLPK_Matrix) :: A
    integer :: rank, ierr

    if (A%ctxt .ne. -1) then
        call blacs_gridexit(A%ctxt)
    endif
end subroutine

!*******************************************************************************
! end subroutine release_context
!*******************************************************************************



!*******************************************************************************
! These two subroutines are internal, just present to reduce the amount of
! boilerplate code.
!*******************************************************************************

! Initialize the descriptor for the SLPK_Matrix A. This just avoids getting
! error messages from descinit by not calling it if this process is not part of
! A's context. If A has junk data then you will probably get an invalid argument
! error from descinit that is a warning flag.
subroutine maybe_descinit(desc, A)
    type(SLPK_Matrix), intent(in) :: A
    integer, intent(out) :: desc(9)
    integer :: info

    if (A%ctxt .eq. -1) then
        desc(1) = 1
        desc(2) = -1
    else
        call descinit(desc, A%m, A%n, A%mb, A%nb, 0, 0, A%ctxt, A%mm, info)
    endif
end subroutine

! Get a Fortran (fat) pointer to the data segment of A.
subroutine get_ptr(dst, A)
    real(REAL_KIND), pointer :: dst(:, :)
    type(SLPK_Matrix), intent(in) :: A

    dst => null()
    if (A%ctxt .ne. -1) then
        call c_f_pointer(A%mdata, dst, [A%mm, A%mn])
    endif
end subroutine

function make_local_context(src)
    integer :: make_local_context
    integer, intent(in) :: src
    integer :: usermap(1, 1)

    call blacs_get(0, 0, make_local_context)
    usermap(1, 1) = src
    call blacs_gridmap(make_local_context, usermap, 1, 1, 1)
    return
end function
!*******************************************************************************
!*******************************************************************************



!*******************************************************************************
subroutine copy_matrix(dst, dsti, dstj, src, srci, srcj, m, n) bind(C)
    type(SLPK_Matrix), intent(inout) :: dst
    type(SLPK_Matrix), intent(in) :: src
    integer(C_INT), value :: dsti, dstj, srci, srcj, m, n

    integer :: dst_desc(9), src_desc(9), info
    real(REAL_KIND), pointer :: dstdata(:, :), srcdata(:, :)

    call check_init()
    call maybe_descinit(dst_desc, dst)
    call maybe_descinit(src_desc, src)
    call get_ptr(dstdata, dst)
    call get_ptr(srcdata, src)

    call pdgemr2d(m, n, srcdata, srci, srcj, src_desc, &
                & dstdata, dsti, dstj, dst_desc, GLOBAL_CTXT)
end subroutine

recursive subroutine transpose_submatrix(dst, dsti, dstj, src, srci, srcj, m, &
                                       & n) bind(C)
    use ISO_FORTRAN_ENV, only: error_unit
    use MPI
    type(SLPK_Matrix), intent(in) :: src
    type(SLPK_Matrix), intent(inout) :: dst
    integer(C_INT), value :: dsti, dstj, srci, srcj, m, n
    integer :: rank, ierr, dst_desc(9), src_desc(9)
    real(REAL_KIND), pointer :: srcdata(:, :), dstdata(:, :)
    type(SLPK_Matrix) :: tmp

    if (dst%ctxt .ne. src%ctxt) then
        call make_similar_matrix(tmp, n, m, src%ctxt, src%nb, src%mb)
        write(*, *) "Recursive call here"
        call transpose_submatrix(tmp, 1, 1, src, srci, srcj, m, n)
        call copy_matrix(dst, dsti, dstj, tmp, 1, 1, n, m)
        call free_matrix_data(tmp)
        return
    endif

    call maybe_descinit(dst_desc, dst)
    call maybe_descinit(src_desc, src)
    call get_ptr(dstdata, dst)
    call get_ptr(srcdata, src)

    write(*, *) "Calling pdgeadd"
    call pdgeadd('T', m, n, 1.0_REAL_KIND, srcdata, srci, srcj, src_desc, &
               & 0.0_REAL_KIND, dstdata, dsti, dstj, dst_desc)
end subroutine

subroutine create_local_matrix(A, x, m, n, src) bind(C)
    use MPI, only: MPI_Comm_rank, MPI_COMM_WORLD
     
    type(SLPK_Matrix), intent(out) :: A
    type(C_PTR), value :: x
    integer(C_INT), value :: m, n, src
    integer :: rank, ierr

    call check_init()
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    A%ctxt = make_local_context(src)
    A%m = m; A%n = n
    A%mdata = C_NULL_PTR
    if (rank .eq. src) then
        A%mb = m; A%mm = m
        A%nb = n; A%mn = n
        A%nprow = 1; A%npcol = 1;
        A%pi = 0; A%pj = 0
        A%mdata = x
    endif
end subroutine

subroutine factorize(mgr) bind(C)
    use mpi
    use ISO_FORTRAN_ENV, only: error_unit

    interface
        subroutine factorize_prep(mgr) bind(C)
            import SVDManager
            type(SVDManager), intent(inout) :: mgr
        end subroutine
    end interface

    type(SVDManager) :: mgr
    integer :: mrank, ierr

    real(REAL_KIND), pointer :: Udata(:, :), Vdata(:, :), Sdata(:), Adata(:, :)
    type(SLPK_Matrix), pointer :: A, U, V, S
    integer :: desca(9), descu(9), descv(9), lwork
    character :: dou, dov
    real(REAL_KIND), allocatable :: work(:)
    real(REAL_KIND) :: bestwork(1)

    call MPI_Comm_rank(MPI_COMM_WORLD, mrank, ierr)
    if (mgr%done .ne. 0) then
        if (mrank .eq. 0) then
            write(*, *) "Multiple call to factorize(mgr), doing nothing"
        endif
        return
    endif

    call c_f_pointer(mgr%A, A)
    if (A%ctxt .eq. -1) then
        write(error_unit, *) "Rank ", mrank, ": This process not associated &
                             &with A, but we assume all processes hold part of A"
    endif
    call descinit(desca, A%m, A%n, A%mb, A%nb, 0, 0, A%ctxt, A%mm, ierr)
    call c_f_pointer(A%mdata, Adata, [A%mm, A%mn])

    call factorize_prep(mgr)
    call c_f_pointer(mgr%S, Sdata, [min(A%m, A%n)])
    dou = 'N'
    if (mgr%dou .ne. 0) then
        dou = 'V'
        call c_f_pointer(mgr%U, U)
        call maybe_descinit(descu, U)
        call get_ptr(Udata, U)
    endif

    dov = 'N'
    if (mgr%dov .ne. 0) then
        dov = 'V'
        call c_f_pointer(mgr%V, V)
        call maybe_descinit(descv, V)
        call get_ptr(Vdata, V)
    endif

    ! The first call to pdgesvd is a workspace query.
    call pdgesvd(dou, dov, A%m, A%n, Adata, 1, 1, desca, Sdata, Udata, 1, 1, &
               & descu, Vdata, 1, 1, descv, bestwork, -1, ierr)
    lwork = bestwork(1)
    allocate(work(lwork))
    call pdgesvd(dou, dov, A%m, A%n, Adata, 1, 1, desca, Sdata, Udata, 1, 1, &
               & descu, Vdata, 1, 1, descv, work, lwork, ierr)
    
    if (mrank .eq. 0) then
        if (ierr .lt. 0) then
            write(error_unit, *) "SVD: Warning: parameter ", -ierr, " had an illegal value"
        elseif (ierr .gt. 0) then
            write(error_unit, *) "SVD: Warning: DBDSQR did not converge"
            if (ierr .eq. min(A%m, A%n)+1) then
                write(error_unit, *) "Singular values not identical across processes"
            endif
        endif
    endif
    deallocate(work)
    mgr%done = 1
end subroutine

end module
