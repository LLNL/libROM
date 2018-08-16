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



AC_INIT(README_BSD_NOTICE)
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
AC_ARG_WITH([CXX],
   [AS_HELP_STRING([--with-CXX=ARG],[manually set C++ compiler to ARG])]
   )
AC_ARG_WITH([AR],
   [AS_HELP_STRING([--with-AR=ARG],[manually set AR to ARG])]
   )

dnl
dnl Initialize common macros to empty strings
dnl

CXX_FLAGS="$CXXFLAGS"
LD_FLAGS="$LDFLAGS"

INC_PATH=
LIB_PATH=
LIB_NAME=

dnl
dnl Guess the machine architecture and set the variable ARCH
dnl

CASC_GUESS_ARCH
CASC_AC_LOG_VAR(ARCH)

dnl
dnl Set up the C++ compiler.
dnl

if test -n "$with_CXX"; then
   CXX=$with_CXX
else
   case "$host_os" in
     linux*) CCC=" g++ ";;
     solaris*) CCC=" g++ ";;
     osf*) CCC=" g++ ";;
     *) CCC=
   esac
fi
CASC_AC_LOG_VAR(with_CXX CXX CCC, before running compiler search)
AC_PROG_CXX($CCC)
CASC_AC_LOG_VAR(with_CXX CXX CCC, after running compiler search)

dnl
dnl Check for C++11 support (required for Elemental)
dnl

AX_CXX_COMPILE_STDCXX_11([], [mandatory])

dnl
dnl Check for Fortran compiler (required to detect Fortran name-mangling
dnl convention); prefer a Fortran 90 compiler for LAPACK & BLAS.
dnl

if test -n "$with_FC"; then
   FC=$with_FC
else
   case "$host_os" in
     linux*) FCC=" gfortran ";;
     solaris*) FCC=" gfortran ";;
     osf*) FCC=" gfortran ";;
     *) FCC=
   esac
fi
CASC_AC_LOG_VAR(with_FC FC FCC, before running compiler search)
AC_PROG_FC($FCC)
CASC_AC_LOG_VAR(with_FC FC FCC, after running compiler search)

dnl
dnl Generate macro wrappers for Fortran functions to make the Fortran
dnl name-mangling convention portable.
dnl

AC_FC_WRAPPERS

dnl
dnl A way to manually add flags to the loader (ld) command.
dnl
AC_ARG_WITH([extra-ld-flags],
[AS_HELP_STRING([--with-extra-ld-flags=STRING],
   [Specify additional flags for the loader.])],
,with_extra_ld_flags='no')
case "$with_extra_ld_flags" in
  no)
    : No op
  ;;
  yes)
    : No op
  ;;
  *)
    LD_FLAGS="$LD_FLAGS $with_extra_ld_flags"
  ;;
esac


dnl
dnl Find ar for archiving
dnl

AC_CHECK_PROGS(AR, $with_AR ar,)
if test -z "$AR"; then
   AC_MSG_ERROR([ar archiver not found - specify with --with-AR=ARG])
fi

dnl
dnl Determine what compilers are being used.
dnl

CASC_INFO_CC_CXX_ID
# Warn user about unsupported compilers:
case $CXX_ID in
gnu|xlc|intel) ;;
*)
  AC_MSG_WARN([
You are using a C++ compiler unsupported by CAROM
You are on your own!
])
  ;;
esac

dnl
dnl Check various C++ compiler features.
dnl

case $CXX_ID in
   gnu)
      CXX_FLAGS="-Wall -Wextra -Wconversion $CXX_FLAGS"
      case $ARCH in
         IRIX64)
            CXX_FLAGS="-mabi=64 $CXX_FLAGS"
            ;;
         solaris)
            ;;
         CYGWIN)
            CXX_FLAGS="$CXX_FLAGS"
            ;;
      esac
      ;;
   xlc)
      CXX_FLAGS="-qrtti -g $CXX_FLAGS"
      ;;
   intel)
      # Ignore some Intel warnings that seem harmless and very hard to get rid of
      CXX_FLAGS="$CXX_FLAGS -Wall -w1 -wd981 -wd383 -wd1418 -wd522"
      ;;
    *)
      ;;
esac

AC_ARG_ENABLE([opt],
   [AS_HELP_STRING([--enable-opt=STRING],
      [Set compiler optimization switches.])],
[case "$enable_opt" in
   no) opt_switches= ;;
   yes)
      case $CXX_ID in
         gnu)   opt_switches="-O2" ;;
         intel) opt_switches="-O3" ;;
         xlc)   opt_switches="-O3" ;;
         *)     opt_switches="-O" ;;
      esac ;;
   *) opt_switches="$enable_opt" ;;
esac],
opt_switches=)

if test -n "$opt_switches" ; then
   CXX_FLAGS="$opt_switches $CXX_FLAGS"
fi

AC_ARG_ENABLE([debug],
[AS_HELP_STRING([--enable-debug=STRING],[Set compiler debug switches.])],
[case "$enable_debug" in
   no) debug_switches= ;;
   yes)
      case $CXX_ID in
         *)   debug_switches="-g" ;;
      esac ;;
   *) debug_switches="$enable_debug" ;;
esac
],
debug_switches=-g)

case "$enable_debug" in
   no) ;;
   *)
   # Enable assertion checking and disable inlining when debugging is enabled
   # if user does not override.
   if test "$user_specified_check_assertions" != "yes"; then
      enable_check_assertions="yes"
   fi
esac

dnl
dnl Check if compiling with assertions on
dnl
if test "$enable_check_assertions" = "yes"; then
   AC_DEFINE([DEBUG_CHECK_ASSERTIONS],[1],[Enable assertion checking])
fi

if test -n "$debug_switches" ; then
   CXX_FLAGS="$debug_switches $CXX_FLAGS"
fi

CASC_AC_LOG_VAR(debug_switches opt_switches)

dnl
dnl Set variables for lapack support
dnl
CASC_SUPPORT_LAPACK
if test "${lapack_PREFIX+set}" = set ||
   test "${lapack_INCLUDES+set}" = set ||
   test "${lapack_LIBS+set}" = set ; then
  SPLIT_LIBS_STRING(lapack_LIBS,lapack_LIB_PATH,lapack_LIB_NAME)
  LIB_PATH="$lapack_LIB_PATH $LIB_PATH"
  LIB_NAME="$lapack_LIB_NAME $LIB_NAME"
  INC_PATH="$lapack_INCLUDES $INC_PATH"
else
  AC_MSG_ERROR([LAPACK is required.])
fi
CASC_AC_LOG_VAR(lapack_LIB_PATH lapack_LIB_NAME LIB_PATH LIB_NAME)

dnl
dnl Determine the proper MPI include and lib paths.
dnl
if test "${with_mpi-include+set}" = set || test "${with_mpi-lib-dirs+set}" = set || test "${with_mpi_libs+set}" = set; then
   CASC_FIND_MPI

   if test -n "$MPIINCLUDE" -o -n "$MPILIBS" -o -n "$MPILIBDIRS"; then
     INC_PATH="$MPIINCLUDE $INC_PATH"
     LIB_NAME="$MPILIBS $LIB_NAME"
     LIB_PATH="$MPILIBDIRS $LIB_PATH"
     CASC_AC_LOG_VAR(MPIINCLUDE MPILIBS MPILIBDIRS, After successfully finding MPI using casc_find_mpi)
   fi
fi

dnl
dnl Support linking with Z lib if HDF requires it
dnl
CHECK_ZLIB()
SPLIT_LIBS_STRING(zlib_LIBS,zlib_LIB_PATH,zlib_LIB_NAME)
LIB_PATH="$zlib_LIB_PATH $LIB_PATH"
LIB_NAME="$zlib_LIB_NAME $LIB_NAME"
INC_PATH="$zlib_INCLUDES $INC_PATH"

dnl
dnl Support linking with SZ lib if HDF requires it
dnl
CHECK_SZLIB()
SPLIT_LIBS_STRING(szlib_LIBS,szlib_LIB_PATH,szlib_LIB_NAME)
LIB_PATH="$szlib_LIB_PATH $LIB_PATH"
LIB_NAME="$szlib_LIB_NAME $LIB_NAME"
INC_PATH="$szlib_INCLUDES $INC_PATH"

dnl
dnl Set variables for hdf5 support
dnl
CASC_SUPPORT_HDF5
if test "${hdf5_PREFIX+set}" = set; then
  SPLIT_LIBS_STRING(hdf5_LIBS,hdf5_LIB_PATH,hdf5_LIB_NAME)
  LIB_PATH="$hdf5_LIB_PATH $LIB_PATH"
  LIB_NAME="$hdf5_LIB_NAME $LIB_NAME"
  INC_PATH="$hdf5_INCLUDES $INC_PATH"
else
  AC_MSG_ERROR([HDF5 not found - you must specify a valid location via --with-hdf5=PATH.])
fi
CASC_AC_LOG_VAR(
[hdf5_INCLUDES INC_PATH hdf5_LIB_PATH hdf5_LIB_NAME LIB_PATH LIB_NAME])

# Defines elemental_PREFIX elemental_INCLUDES and elemental_LIBS if
# with-elemental is specified.
AC_ARG_WITH(elemental,
[  --with-elemental[=PATH]  Use Elemental and optionally specify where Elemental is installed.],
, with_elemental=no)

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

if test "${elemental_PREFIX+set}" = set; then
  SPLIT_LIBS_STRING(elemental_LIBS,elemental_LIB_PATH,elemental_LIB_NAME)
  LIB_PATH="$elemental_LIB_PATH $LIB_PATH"
  LIB_NAME="$elemental_LIB_NAME $LIB_NAME"
  INC_PATH="$elemental_INCLUDES $INC_PATH"
  CXX_FLAGS="$CXX_FLAGS -DCAROM_HAS_ELEMENTAL"
else
  AC_MSG_NOTICE([Elemental not found - disabling distributed QR decomposition code])
fi
CASC_AC_LOG_VAR(
[elemental_INCLUDES INC_PATH elemental_LIB_PATH elemental_LIB_NAME LIB_PATH LIB_NAME])


dnl
dnl Detect Google Test Framework support
dnl

AC_ARG_WITH(gtest,
[  --with-gtest[=PATH]  Use Google Test Framework and optionally specify where Google Test is installed.],
, with_gtest=no)

case "$with_gtest" in
  no)
    AC_MSG_NOTICE([configuring without Google Test Framework support])
    : Do nothing
  ;;
  yes)
    # Google Test Framework install path was not specified.
    # Look in a couple of standard locations to probe if
    # Google Test Framework header files are there.
    AC_MSG_CHECKING([for Google Test Framework installation])
    for dir in /usr /usr/local; do
      if test -f ${dir}/include/gtest/gtest.h; then
        gtest_PREFIX=${dir}
        break
      fi
    done
    AC_MSG_RESULT([$gtest_PREFIX])
  ;;
  *)
    # Google Test Framework install path was specified.
    AC_MSG_CHECKING([for Google Test Framework installation])

    if test -f ${with_gtest}/include/gtest/gtest.h; then
        gtest_PREFIX=$with_gtest
        gtest_INCLUDES="-I${gtest_PREFIX}/include"
        gtest_LIBS="-L${gtest_PREFIX}/lib -lgtest -lpthread"
        AC_MSG_RESULT([$gtest_PREFIX])
    else
        AC_MSG_RESULT([$gtest_PREFIX])
        AC_MSG_ERROR([Google Test Framework not found in $with_gtest])
    fi
  ;;
esac

if test "${gtest_PREFIX+set}" = set; then
  SPLIT_LIBS_STRING(gtest_LIBS,gtest_LIB_PATH,gtest_LIB_NAME)
  LIB_PATH="$gtest_LIB_PATH $LIB_PATH"
  LIB_NAME="$gtest_LIB_NAME $LIB_NAME"
  INC_PATH="$gtest_INCLUDES $INC_PATH"
  CXX_FLAGS="$CXX_FLAGS -DCAROM_HAS_ELEMENTAL"
else
  AC_MSG_NOTICE([Google Test Framework not found - disabling unit tests.])
fi
CASC_AC_LOG_VAR(
[gtest_INCLUDES INC_PATH gtest_LIB_PATH gtest_LIB_NAME LIB_PATH LIB_NAME])

dnl
dnl Support rpath syntax by adding it to the LIB_PATH variable.
dnl
CASC_AC_LOG_VAR(LIB_PATH, before adding rpath syntax)
CASC_LIBS_ADD_RPATH(LIB_PATH,LIB_PATH,1)
CASC_AC_LOG_VAR(LIB_PATH, after adding rpath syntax)

dnl
dnl Set DOXYGEN_BIN to the doxygen binary.
dnl

CASC_PATH_PROG(DOXYGEN_BIN,doxygen)
CASC_AC_LOG_VAR(DOXYGEN_BIN)

dnl
dnl Set DOT_PATH to the directory of the dot binary, for use with doxygen.
dnl

AC_ARG_WITH([dot-dir],
[AS_HELP_STRING([--with-dot-dir=DIR],
   [Specify bin directory where dot binaries live.])],,
unset with_dot_dir)
if   test "${with_dot_dir}" = yes || test ! "${with_dot_dir+set}" = set; then
  # Find dot path, and use if found
  DOT_PATH=`which dot >& /dev/null`
  if test -z "${DOT_PATH}"; then
    unset DOT_PATH
    HAVE_DOT=NO
  else
    [DOT_PATH=`echo $DOT_PATH | sed 's!/[^/]\{1,\}$!!']`
    HAVE_DOT=YES
  fi
elif test "${with_dot_dir}" = no; then
  # Disable dot.
  unset DOT_PATH
  HAVE_DOT=NO
else
  # Make sure valid path is specified.  Enable dot.
  if test ! -d ${with_dot_dir}; then
     AC_MSG_WARN([
Invalid path given for dot executable: ${with_dot_dir}.
Path does not exist.
You can expect doxygen documentation generation to not work properly.
]);
  fi
  DOT_PATH=${with_dot_dir}
  HAVE_DOT=YES
fi
CASC_AC_LOG_VAR(with_dot_dir HAVE_DOT DOT_PATH)
AC_SUBST(HAVE_DOT)
AC_SUBST(DOT_PATH)

dnl
dnl Set the directory for doxygen documentation.
dnl
AC_ARG_WITH([dox-output-dir],
[AS_HELP_STRING([--with-dox-output-dir=DIR],
   [Specify full where Doxygen documentation is generated.  Default is the docs
   directory in the compile tree.])])
if test -n "$with_dox_output_dir"; then
  if echo "$with_dox_output_dir" | grep -v '^/' ; then
    AC_MSG_ERROR(You must give a full path with --with-dox-output-dir=DIR)
  fi
  dox_output_dir="$with_dox_output_dir"
  # Make sure dox_output_dir ends with / so because the doxygen
  # configuration files simply append a subdirectory name.
  [dox_output_dir="`echo $dox_output_dir | sed 's:\([^/]\)$:\1/:'`"]
fi
AC_SUBST(dox_output_dir)

dnl
dnl Set ROM_PERL to the perl binary for use in doxygen.
dnl

CASC_PATH_PROG(ROM_PERL,perl)
CASC_AC_LOG_VAR(ROM_PERL)

AC_SUBST(AR)

AC_SUBST(CXX)
AC_SUBST(CXX_FLAGS)

AC_SUBST(LD_FLAGS)

AC_SUBST(INC_PATH)
AC_SUBST(LIB_PATH)
AC_SUBST(LIB_NAME)

AC_OUTPUT([Makefile docs/Makefile docs/Doxyfile],[],[])
