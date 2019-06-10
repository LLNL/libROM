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

dnl Import some utility macros
m4_include([m4/configure-log.m4])
m4_include([m4/split-libs-string.m4])

AC_CONFIG_HEADER(CAROM_config.h)
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
dnl These will be built up with the flags from each package as they
dnl are discovered.  Final values are used in Makefile.am.
dnl

INC_PATH=
LIB_PATH=
LIB_NAME=

dnl Check for some basic functionality
AC_PROG_CXX
AC_PROG_RANLIB
AC_FC_WRAPPERS
AM_PROG_AR

dnl
dnl Check for C++11 support (required for Elemental)
dnl

m4_include([m4/ax_require_defined.m4])
m4_include([m4/ax_cxx_compile_stdcxx.m4])
m4_include([m4/ax_cxx_compile_stdcxx_11.m4])
AX_CXX_COMPILE_STDCXX_11([], [mandatory])
dnl Above macro sets cxxflags for some reason. Unset it.
CXXFLAGS=

dnl
dnl Check for the required dependencies: hdf5, mpi, lapack/blas
dnl

dnl Begin HDF5 section
m4_include([m4/ax_lib_hdf5.m4])

AX_LIB_HDF5([serial])	 

if test "$with_hdf5" = "no"; then
	AC_MSG_ERROR([Unable to find HDF5, we need serial HDF5.])
else

LIBS="$LIBS $HDF5_LIBS"
LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
CPPFLAGS="$CPPFLAGS $HDF5_CPPFLAGS"
fi
dnl End HDF5 section

dnl Begin mpi compiler wrappers section
m4_include([m4/ax_mpi.m4])
AC_ARG_WITH(mpi, [  --with-mpi=<path>              prefix of MPI installation],
                 [  MPICXX="$withval/bin/mpic++"
                    CPPFLAGS="$CPPFLAGS -I $withval/include"
                    LDFLAGS="$LDFLAGS -L$withval/lib"])

if test "x$with_mpi" != "xno"; then
   AX_MPI([], AC_MSG_ERROR([Could not find required MPI library.
                         Consider using --with-mpi=<my-path-to-mpi>.]))
   CXX="$MPICXX"
   LIBS="$MPILIBS $LIBS"
fi
dnl End mpi compiler wrappers section

dnl Start lapack/blas section
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
dnl End lapack/blas section

dnl
dnl Check for optional dependencies: elemental, google-test
dnl

dnl Begin Elemental check

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
dnl End Elemental check

dnl Begin Google Test Framework support
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

if test "x${with_gtest}" != xno; then
    if test "${gtest_PREFIX+set}" = set; then
        SPLIT_LIBS_STRING(gtest_LIBS,gtest_LIB_PATH,gtest_LIB_NAME)
        LIB_PATH="$gtest_LIB_PATH $LIB_PATH"
        LIB_NAME="$gtest_LIB_NAME $LIB_NAME"
        INC_PATH="$gtest_INCLUDES $INC_PATH"
        CXXFLAGS="$CXXFLAGS -DCAROM_HAS_GTEST"
    else
        AC_MSG_NOTICE([Google Test Framework not found - disabling unit tests.])
    fi
    CASC_AC_LOG_VAR(
    [gtest_INCLUDES INC_PATH gtest_LIB_PATH gtest_LIB_NAME LIB_PATH LIB_NAME])
fi
dnl End Google Test Framework support


AC_ARG_ENABLE([opt],
   [AS_HELP_STRING([--enable-opt=STRING],
      [Set compiler optimization switches.])],
[case "$enable_opt" in
   no) opt_switches= ;;
   yes) opt_switches="-O3";;
   *) opt_switches="$enable_opt" ;;
esac],
opt_switches=)

if test -n "$opt_switches" ; then
   CXXFLAGS="$opt_switches $CXXFLAGS"
fi

AC_ARG_ENABLE([debug],
[AS_HELP_STRING([--enable-debug=STRING],[Set compiler debug switches.])],
[case "$enable_debug" in
   no) debug_switches= ;;
   yes) debug_switches="-g" ;;
   *) debug_switches="$enable_debug" ;;
esac
],
debug_switches=)

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
   CXXFLAGS="$debug_switches $CXXFLAGS"
fi

CASC_AC_LOG_VAR(debug_switches opt_switches)

dnl
dnl Support rpath syntax by adding it to the LIB_PATH variable.
dnl
m4_include([m4/compiler-id.m4])
m4_include([m4/support-rpath.m4])
CASC_AC_LOG_VAR(LIB_PATH, before adding rpath syntax)
CASC_LIBS_ADD_RPATH(LIB_PATH,LIB_PATH,1)
CASC_AC_LOG_VAR(LIB_PATH, after adding rpath syntax)

dnl Begin Doxygen section

AC_CHECK_PROGS([DOXYGEN_BIN], [doxygen])
if test -z "$DOXYGEN_BIN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi
CASC_AC_LOG_VAR(DOXYGEN_BIN)

dnl
dnl Set DOT_PATH to the directory of the dot binary, for use with doxygen.
dnl dot is a graph layout program, used for making class hierarchies.
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

AC_CHECK_PROGS([ROM_PERL], [perl])
if test -z "$ROM_PERL";
   then AC_MSG_WARN([perl not found - continuing without perl support])
fi
CASC_AC_LOG_VAR(ROM_PERL)

dnl End doxygen section

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
