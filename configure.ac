dnl
dnl Copyright (c) 2013-2017, Lawrence Livermore National Security, LLC.
dnl Produced at the Lawrence Livermore National Laboratory
dnl Written by William Arrighi wjarrighi@llnl.gov
dnl CODE-686965
dnl All rights reserved.
dnl
dnl This file is part of libROM.
dnl For details, see https://computation.llnl.gov/librom
dnl Please also read README_BSD_NOTICE.
dnl
dnl Redistribution and use in source and binary forms, with or without
dnl modifications, are permitted provided that the following conditions are
dnl met:
dnl
dnl    o Redistributions of source code must retain the above copyright notice,
dnl      this list of conditions and the disclaimer below.
dnl    o Redistribution in binary form must reproduce the above copyright
dnl      notice, this list of conditions and the disclaimer (as noted below) in
dnl      the documentation and/or other materials provided with the
dnl      distribution.
dnl    o Neither the name of the LLNS/LLNL nor the names of its contributors
dnl      may be used to endorse or promote products derived from this software
dnl      without specific prior written permission.
dnl
dnl THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
dnl AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
dnl IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
dnl ARE DISCLAIMED.  IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
dnl LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
dnl DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
dnl (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
dnl SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
dnl CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
dnl LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
dnl OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OR
dnl SUCH DAMAGE.
dnl

dnl Description: autoconf configuration file


AC_INIT(libROM, 0.1)

# Initialize Automake
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects])

AC_CONFIG_HEADER(CAROM_config.h:config/CAROM_config.h.in)
dnl
dnl Autoheader stuff
AH_TOP([

#ifndef INCLUDED_CAROM_CONFIG_H
#define INCLUDED_CAROM_CONFIG_H

#define CAROM_VERSION_MAJOR 0
#define CAROM_VERSION_MINOR 1
#define CAROM_VERSION_PATCHLEVEL 0
])

AH_BOTTOM([

#endif
])

AC_ARG_ENABLE([check_assertions],
   [AS_HELP_STRING([--enable-check-assertions],
      [build with assertion checking on (default is disabled; is enabled for when debugging is enabled)])],
   [user_specified_check_assertions="yes"],
   [enable_check_assertions="no"]
   )

dnl
dnl Initialize common macros to empty strings
dnl

INC_PATH=
LIB_PATH=
LIB_NAME=

AC_PROG_CXX
AC_PROG_RANLIB
AC_FC_WRAPPERS
AM_PROG_AR

dnl
dnl Check for C++11 support (required for Elemental)
dnl

AX_CXX_COMPILE_STDCXX_11([], [mandatory])

dnl Now pick up the required dependencies
dnl hdf5
dnl mpi

dnl Begin HDF5 section
m4_include([m4/ax_lib_hdf5.m4])

AX_LIB_HDF5([serial])	 

if test "$with_hdf5" = "no"; then
	AC_MSG_ERROR([Unable to find HDF5, we need serial HDF5.])
else

dnl Manually add rpaths for each linker search path
LIBS="$LIBS $HDF5_LIBS"
HDF5_RPATH=`echo $HDF5_LDFLAGS | sed 's/-L/-Wl,-rpath,/g'`

LDFLAGS="$LDFLAGS $HDF5_LDFLAGS $HDF5_RPATH"
CPPFLAGS="$CPPFLAGS $HDF5_CPPFLAGS"
fi
dnl End HDF5 section

dnl Begin mpi compiler wrappers section
m4_include([m4/ax_mpi.m4])
AC_ARG_WITH(mpi, [  --with-mpi=<path>              prefix of MPI installation],
                 [  MPICXX="$withval/bin/mpic++"
                    CPPFLAGS="$CPPFLAGS -I $withval/include"
                    LDFLAGS="$LDFLAGS -L $withval/lib"])

if test "x$with_mpi" != "xno"; then
   AX_MPI([], AC_MSG_ERROR([Could not find required MPI library.
                         Consider using --with-mpi=<my-path-to-mpi>.]))
   CXX="$MPICXX"
   LIBS="$MPILIBS $LIBS"
fi
dnl End mpi compiler wrappers section

dnl Start lapack section
m4_include([m4/ax_blas.m4])
m4_include([m4/ax_lapack.m4])

AX_BLAS
if test "$ax_blas_ok" = "no"; then
   AC_MSG_ERROR([Cannot find BLAS libraries])
else
   SPLIT_LIBS_STRING(BLAS_LIBS,path1,name1)
   SPLIT_LIBS_STRING(LIBS,path2,name2)
   SPLIT_LIBS_STRING(FLIBS,path3,name3)
   LIB_PATH="$LIB_PATH $path1 $path2 $path3"
   LIB_NAME="$LIB_NAME $name1 $name2 $name3"
fi

AX_LAPACK
if test "$ax_lapack_ok" = "no"; then
   AC_MSG_ERROR([Cannot find LAPACK libraries])
fi
End lapack section


dnl START


dnl dnl
dnl dnl Check for Fortran compiler (required to detect Fortran name-mangling
dnl dnl convention); prefer a Fortran 90 compiler for LAPACK & BLAS.
dnl dnl

dnl if test -n "$with_FC"; then
dnl    FC=$with_FC
dnl else
dnl    case "$host_os" in
dnl      linux*) FCC=" gfortran ";;
dnl      solaris*) FCC=" gfortran ";;
dnl      osf*) FCC=" gfortran ";;
dnl      *) FCC=
dnl    esac
dnl fi
dnl CASC_AC_LOG_VAR(with_FC FC FCC, before running compiler search)
dnl dnl The intention behind this search order for Fortran compilers is to
dnl dnl put the user-specified compiler (if available) first. If the OS is
dnl dnl detected properly, then common OSes will put gfortran first. If
dnl dnl the OS is exotic, then the OS detection above will set FCC to the
dnl dnl empty string, in which case, the search order proceeds from
dnl dnl (probably) rarer compilers to more common compilers.  Fewer
dnl dnl systems have IBM's xl compilers than PGI. If a system has a PGI
dnl dnl compiler, it usually also has a Intel compiler, but the reverse
dnl dnl isn't necessarily true. Intel compilers are still pretty common,
dnl dnl but flang is not. However, flang is free, and the xl, PGI, and
dnl dnl Intel compilers generally require a license, so we put flang after
dnl dnl the compilers that require a license, and before gfortran, which
dnl dnl is very common.

dnl AC_PROG_FC([$FCC xlflang xlf90 pg90 ifort flang gfortran])
dnl CASC_AC_LOG_VAR(with_FC FC FCC, after running compiler search)

dnl dnl Determine Fortran libraries for linking; on some systems, these
dnl dnl are needed for linking LAPACK, which is Fortran-based. This logic
dnl dnl doesn't work on the CORAL EA systems, which use a beta xl compiler

dnl case "$FC" in
dnl     *xlf) ;;
dnl     *xlf90) ;;
dnl     *xlf-gpu) ;;
dnl     *) dnl This branch should include xlflang, for which Fortran linkage
dnl        dnl detection works.
dnl         AC_FC_LIBRARY_LDFLAGS
dnl         SPLIT_LIBS_STRING(FCLIBS,fc_LIB_PATH,fc_LIB_NAME)
dnl     ;;
dnl esac


dnl dnl
dnl dnl Generate macro wrappers for Fortran functions to make the Fortran
dnl dnl name-mangling convention portable.
dnl dnl

dnl AC_FC_WRAPPERS

dnl dnl
dnl dnl A way to manually add flags to the loader (ld) command.
dnl dnl
dnl AC_ARG_WITH([extra-ld-flags],
dnl [AS_HELP_STRING([--with-extra-ld-flags=STRING],
dnl    [Specify additional flags for the loader.])],
dnl ,with_extra_ld_flags='no')
dnl case "$with_extra_ld_flags" in
dnl   no)
dnl     : No op
dnl   ;;
dnl   yes)
dnl     : No op
dnl   ;;
dnl   *)
dnl     LDFLAGS="$LDFLAGS $with_extra_ld_flags"
dnl   ;;
dnl esac

dnl dnl
dnl dnl Find ranlib
dnl dnl
dnl AC_PROG_RANLIB

dnl dnl
dnl dnl Find ar for archiving
dnl dnl

dnl AC_CHECK_PROGS(AR, $with_AR ar,)
dnl if test -z "$AR"; then
dnl    AC_MSG_ERROR([ar archiver not found - specify with --with-AR=ARG])
dnl fi

dnl dnl
dnl dnl Determine what compilers are being used.
dnl dnl

dnl CASC_INFO_CC_CXX_ID
dnl # Warn user about unsupported compilers:
dnl case $CXX_ID in
dnl gnu|xlc|intel) ;;
dnl *)
dnl   AC_MSG_WARN([
dnl You are using a C++ compiler unsupported by CAROM
dnl You are on your own!
dnl ])
dnl   ;;
dnl esac

dnl dnl
dnl dnl Check various C++ compiler features.
dnl dnl

dnl case $CXX_ID in
dnl    gnu)
dnl       CXXFLAGS="-Wall -Wextra -Wconversion $CXXFLAGS"
dnl       case $ARCH in
dnl          IRIX64)
dnl             CXXFLAGS="-mabi=64 $CXXFLAGS"
dnl             ;;
dnl          solaris)
dnl             ;;
dnl          CYGWIN)
dnl             CXXFLAGS="$CXXFLAGS"
dnl             ;;
dnl       esac
dnl       ;;
dnl    xlc)
dnl       CXXFLAGS="-qrtti -g $CXXFLAGS"
dnl       ;;
dnl    intel)
dnl       # Ignore some Intel warnings that seem harmless and very hard to get rid of
dnl       CXXFLAGS="$CXXFLAGS -Wall -w1 -wd981 -wd383 -wd1418 -wd522"
dnl       ;;
dnl     *)
dnl       ;;
dnl esac

dnl AC_ARG_ENABLE([opt],
dnl    [AS_HELP_STRING([--enable-opt=STRING],
dnl       [Set compiler optimization switches.])],
dnl [case "$enable_opt" in
dnl    no) opt_switches= ;;
dnl    yes)
dnl       case $CXX_ID in
dnl          gnu)   opt_switches="-O2" ;;
dnl          intel) opt_switches="-O3" ;;
dnl          xlc)   opt_switches="-O3" ;;
dnl          *)     opt_switches="-O" ;;
dnl       esac ;;
dnl    *) opt_switches="$enable_opt" ;;
dnl esac],
dnl opt_switches=)

dnl if test -n "$opt_switches" ; then
dnl    CXXFLAGS="$opt_switches $CXXFLAGS"
dnl fi

dnl AC_ARG_ENABLE([debug],
dnl [AS_HELP_STRING([--enable-debug=STRING],[Set compiler debug switches.])],
dnl [case "$enable_debug" in
dnl    no) debug_switches= ;;
dnl    yes)
dnl       case $CXX_ID in
dnl          *)   debug_switches="-g" ;;
dnl       esac ;;
dnl    *) debug_switches="$enable_debug" ;;
dnl esac
dnl ],
dnl debug_switches=-g)

dnl case "$enable_debug" in
dnl    no) ;;
dnl    *)
dnl    # Enable assertion checking and disable inlining when debugging is enabled
dnl    # if user does not override.
dnl    if test "$user_specified_check_assertions" != "yes"; then
dnl       enable_check_assertions="yes"
dnl    fi
dnl esac

dnl dnl
dnl dnl Check if compiling with assertions on
dnl dnl
dnl if test "$enable_check_assertions" = "yes"; then
dnl    AC_DEFINE([DEBUG_CHECK_ASSERTIONS],[1],[Enable assertion checking])
dnl fi

dnl if test -n "$debug_switches" ; then
dnl    CXXFLAGS="$debug_switches $CXXFLAGS"
dnl fi

dnl CASC_AC_LOG_VAR(debug_switches opt_switches)

dnl dnl Add linkage to Fortran libraries, which is required because
dnl dnl LAPACK is written in Fortran, and it is being linked by a
dnl dnl C++ compiler, which does not include Fortran libraries by default
dnl LIB_PATH="$fc_LIB_PATH $LIB_PATH"
dnl LIB_NAME="$fc_LIB_NAME $LIB_PATH"

dnl dnl
dnl dnl Set variables for lapack support
dnl dnl
dnl CASC_SUPPORT_LAPACK
dnl if test "${lapack_PREFIX+set}" = set ||
dnl    test "${lapack_INCLUDES+set}" = set ||
dnl    test "${lapack_LIBS+set}" = set ; then
dnl   SPLIT_LIBS_STRING(lapack_LIBS,lapack_LIB_PATH,lapack_LIB_NAME)
dnl   LIB_PATH="$lapack_LIB_PATH $LIB_PATH"
dnl   LIB_NAME="$lapack_LIB_NAME $LIB_NAME"
dnl   INC_PATH="$lapack_INCLUDES $INC_PATH"
dnl else
dnl   AC_MSG_ERROR([LAPACK is required.])
dnl fi
dnl CASC_AC_LOG_VAR(lapack_LIB_PATH lapack_LIB_NAME LIB_PATH LIB_NAME)

dnl dnl
dnl dnl Determine the proper MPI include and lib paths.
dnl dnl
dnl if test "${with_mpi-include+set}" = set || test "${with_mpi-lib-dirs+set}" = set || test "${with_mpi_libs+set}" = set; then
dnl    CASC_FIND_MPI

dnl    if test -n "$MPIINCLUDE" -o -n "$MPILIBS" -o -n "$MPILIBDIRS"; then
dnl      INC_PATH="$MPIINCLUDE $INC_PATH"
dnl      LIB_NAME="$MPILIBS $LIB_NAME"
dnl      LIB_PATH="$MPILIBDIRS $LIB_PATH"
dnl      CASC_AC_LOG_VAR(MPIINCLUDE MPILIBS MPILIBDIRS, After successfully finding MPI using casc_find_mpi)
dnl    fi
dnl fi

dnl dnl
dnl dnl Support linking with Z lib if HDF requires it
dnl dnl
dnl CHECK_ZLIB()
dnl SPLIT_LIBS_STRING(zlib_LIBS,zlib_LIB_PATH,zlib_LIB_NAME)
dnl LIB_PATH="$zlib_LIB_PATH $LIB_PATH"
dnl LIB_NAME="$zlib_LIB_NAME $LIB_NAME"
dnl INC_PATH="$zlib_INCLUDES $INC_PATH"

dnl dnl
dnl dnl Support linking with SZ lib if HDF requires it
dnl dnl
dnl CHECK_SZLIB()
dnl SPLIT_LIBS_STRING(szlib_LIBS,szlib_LIB_PATH,szlib_LIB_NAME)
dnl LIB_PATH="$szlib_LIB_PATH $LIB_PATH"
dnl LIB_NAME="$szlib_LIB_NAME $LIB_NAME"
dnl INC_PATH="$szlib_INCLUDES $INC_PATH"

# Defines elemental_PREFIX elemental_INCLUDES and elemental_LIBS if
# with-elemental is specified.
AC_ARG_WITH(elemental,
[  --with-elemental[=PATH]  Use Elemental and optionally specify where Elemental is installed.],
[],
[with_elemental=no])

case "$with_elemental" in
  no)
    AC_MSG_NOTICE([configuring without Elemental support])
    : Do nothing
  ;;
  yes)
    # Elemental install path was not specified.
    # Look in a couple of standard locations to probe if
    # Elemental header files are there.
    AC_MSG_CHECKING([for Elemental installation])
    for dir in /usr /usr/local; do
      if test -f ${dir}/include/El.hpp; then
        elemental_PREFIX=${dir}
	elemental_INCLUDES="-I${elemental_PREFIX}/include"
	if test -e ${dir}/lib/libEl.dylib -o -e ${dir}/lib/libEl.so -o -e ${dir}/lib/libEl.a; then
           elemental_LIBS="-L${elemental_PREFIX}/lib -lEl"
	elif test -e ${dir}/lib64/libEl.dylib -o -e ${dir}/lib64/libEl.so -o -e ${dir}/lib64/libEl.a; then
	   elemental_LIBS="-L${elemental_PREFIX}/lib64 -lEl"
	else
	   AC_MSG_ERROR([Elemental's libEl not found in $elemental_PREFIX])
	fi
        break
      fi
    done
    AC_MSG_RESULT([$elemental_PREFIX])
  ;;
  *)
    # Elemental install path was specified.
    AC_MSG_CHECKING([for Elemental installation])

    if test -f ${with_elemental}/include/El.hpp; then
        elemental_PREFIX=$with_elemental
        elemental_INCLUDES="-I${elemental_PREFIX}/include"
	if test -e ${with_elemental}/lib/libEl.dylib -o -e ${with_elemental}/lib/libEl.so -o -e ${with_elemental}/lib/libEl.a; then
           elemental_LIBS="-L${elemental_PREFIX}/lib -lEl"
	elif test -e ${with_elemental}/lib64/libEl.dylib -o -e ${with_elemental}/lib64/libEl.so -o -e ${with_elemental}/lib64/libEl.a; then
	   elemental_LIBS="-L${elemental_PREFIX}/lib64 -lEl"
	else
	   AC_MSG_ERROR([Elemental's libEl not found in $elemental_PREFIX])
	fi
        AC_MSG_RESULT([$elemental_PREFIX])
    else
        AC_MSG_RESULT([$elemental_PREFIX])
        AC_MSG_ERROR([Elemental not found in $with_elemental])
    fi
  ;;
esac

if test "x${with_elemental}" != xno; then
    if test "${elemental_PREFIX+set}" = set; then
        SPLIT_LIBS_STRING(elemental_LIBS,elemental_LIB_PATH,elemental_LIB_NAME)
        LIB_PATH="$elemental_LIB_PATH $LIB_PATH"
        LIB_NAME="$elemental_LIB_NAME $LIB_NAME"
        INC_PATH="$elemental_INCLUDES $INC_PATH"
        AC_DEFINE(CAROM_HAS_ELEMENTAL,1,[Have Elemental library.])
else
        AC_MSG_NOTICE([Elemental not found - disabling distributed QR decomposition code])
    fi
    
    CASC_AC_LOG_VAR(
    [elemental_INCLUDES INC_PATH elemental_LIB_PATH elemental_LIB_NAME LIB_PATH LIB_NAME])
fi

dnl dnl
dnl dnl Detect Google Test Framework support
dnl dnl

dnl AC_ARG_WITH(gtest,
dnl [  --with-gtest[=PATH]  Use Google Test Framework and optionally specify where Google Test is installed.],
dnl , with_gtest=no)

dnl case "$with_gtest" in
dnl   no)
dnl     AC_MSG_NOTICE([configuring without Google Test Framework support])
dnl     : Do nothing
dnl   ;;
dnl   yes)
dnl     # Google Test Framework install path was not specified.
dnl     # Look in a couple of standard locations to probe if
dnl     # Google Test Framework header files are there.
dnl     AC_MSG_CHECKING([for Google Test Framework installation])
dnl     for dir in /usr /usr/local; do
dnl       if test -f ${dir}/include/gtest/gtest.h; then
dnl         gtest_PREFIX=${dir}
dnl         break
dnl       fi
dnl     done
dnl     AC_MSG_RESULT([$gtest_PREFIX])
dnl   ;;
dnl   *)
dnl     # Google Test Framework install path was specified.
dnl     AC_MSG_CHECKING([for Google Test Framework installation])

dnl     if test -f ${with_gtest}/include/gtest/gtest.h; then
dnl         gtest_PREFIX=$with_gtest
dnl         gtest_INCLUDES="-I${gtest_PREFIX}/include"
dnl         gtest_LIBS="-L${gtest_PREFIX}/lib -lgtest -lpthread"
dnl         AC_MSG_RESULT([$gtest_PREFIX])
dnl     else
dnl         AC_MSG_RESULT([$gtest_PREFIX])
dnl         AC_MSG_ERROR([Google Test Framework not found in $with_gtest])
dnl     fi
dnl   ;;
dnl esac

dnl if test "x${with_gtest}" != xno; then
dnl     if test "${gtest_PREFIX+set}" = set; then
dnl         SPLIT_LIBS_STRING(gtest_LIBS,gtest_LIB_PATH,gtest_LIB_NAME)
dnl         LIB_PATH="$gtest_LIB_PATH $LIB_PATH"
dnl         LIB_NAME="$gtest_LIB_NAME $LIB_NAME"
dnl         INC_PATH="$gtest_INCLUDES $INC_PATH"
dnl         CXXFLAGS="$CXXFLAGS -DCAROM_HAS_GTEST"
dnl     else
dnl         AC_MSG_NOTICE([Google Test Framework not found - disabling unit tests.])
dnl     fi
dnl     CASC_AC_LOG_VAR(
dnl     [gtest_INCLUDES INC_PATH gtest_LIB_PATH gtest_LIB_NAME LIB_PATH LIB_NAME])
dnl fi
dnl dnl
dnl dnl Support rpath syntax by adding it to the LIB_PATH variable.
dnl dnl
dnl CASC_AC_LOG_VAR(LIB_PATH, before adding rpath syntax)
dnl CASC_LIBS_ADD_RPATH(LIB_PATH,LIB_PATH,1)
dnl CASC_AC_LOG_VAR(LIB_PATH, after adding rpath syntax)

dnl dnl
dnl dnl Set DOXYGEN_BIN to the doxygen binary.
dnl dnl

dnl CASC_PATH_PROG(DOXYGEN_BIN,doxygen)
dnl CASC_AC_LOG_VAR(DOXYGEN_BIN)

dnl dnl
dnl dnl Set DOT_PATH to the directory of the dot binary, for use with doxygen.
dnl dnl

dnl AC_ARG_WITH([dot-dir],
dnl [AS_HELP_STRING([--with-dot-dir=DIR],
dnl    [Specify bin directory where dot binaries live.])],,
dnl unset with_dot_dir)
dnl if   test "${with_dot_dir}" = yes || test ! "${with_dot_dir+set}" = set; then
dnl   # Find dot path, and use if found
dnl   DOT_PATH=`which dot >& /dev/null`
dnl   if test -z "${DOT_PATH}"; then
dnl     unset DOT_PATH
dnl     HAVE_DOT=NO
dnl   else
dnl     [DOT_PATH=`echo $DOT_PATH | sed 's!/[^/]\{1,\}$!!']`
dnl     HAVE_DOT=YES
dnl   fi
dnl elif test "${with_dot_dir}" = no; then
dnl   # Disable dot.
dnl   unset DOT_PATH
dnl   HAVE_DOT=NO
dnl else
dnl   # Make sure valid path is specified.  Enable dot.
dnl   if test ! -d ${with_dot_dir}; then
dnl      AC_MSG_WARN([
dnl Invalid path given for dot executable: ${with_dot_dir}.
dnl Path does not exist.
dnl You can expect doxygen documentation generation to not work properly.
dnl ]);
dnl   fi
dnl   DOT_PATH=${with_dot_dir}
dnl   HAVE_DOT=YES
dnl fi
dnl CASC_AC_LOG_VAR(with_dot_dir HAVE_DOT DOT_PATH)
dnl AC_SUBST(HAVE_DOT)
dnl AC_SUBST(DOT_PATH)

dnl dnl
dnl dnl Set the directory for doxygen documentation.
dnl dnl
dnl AC_ARG_WITH([dox-output-dir],
dnl [AS_HELP_STRING([--with-dox-output-dir=DIR],
dnl    [Specify full where Doxygen documentation is generated.  Default is the docs
dnl    directory in the compile tree.])])
dnl if test -n "$with_dox_output_dir"; then
dnl   if echo "$with_dox_output_dir" | grep -v '^/' ; then
dnl     AC_MSG_ERROR(You must give a full path with --with-dox-output-dir=DIR)
dnl   fi
dnl   dox_output_dir="$with_dox_output_dir"
dnl   # Make sure dox_output_dir ends with / so because the doxygen
dnl   # configuration files simply append a subdirectory name.
dnl   [dox_output_dir="`echo $dox_output_dir | sed 's:\([^/]\)$:\1/:'`"]
dnl fi
dnl AC_SUBST(dox_output_dir)

dnl dnl
dnl dnl Set ROM_PERL to the perl binary for use in doxygen.
dnl dnl

dnl CASC_PATH_PROG(ROM_PERL,perl)
dnl CASC_AC_LOG_VAR(ROM_PERL)

dnl END

AC_SUBST(AR)

AC_SUBST(CXX)
AC_SUBST(CXXFLAGS)

AC_SUBST(LDFLAGS)

AC_SUBST(INC_PATH)
AC_SUBST(LIB_PATH)
AC_SUBST(LIB_NAME)

AC_CONFIG_FILES([
   Makefile
   docs/Makefile
   docs/Doxyfile
])
AC_OUTPUT
