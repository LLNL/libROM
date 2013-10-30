dnl $Id$

AC_DEFUN([CASC_VAR_SET_LAPACK],[
dnl Provides support for the lapack library.
dnl
dnl Arguments are:
dnl 1. Name of variable to set to path where lapack are installed.
dnl    Nothig is done if this variable is unset.
dnl 2. Name of the INCLUDES variable similar to the automake INCLUDES variable.
dnl    This variable is modified ONLY if it is NOT set.
dnl 3. Name of the LIBS variable similar to the automake LIBS variable.
dnl    This variable is modified ONLY if it is NOT set.
dnl
dnl If arg1 is defined, assume that the user wants lapack
dnl support.  Do so by assigning arg2 and arg3 if they are not defined.
dnl
if test "${$1+set}" = set ; then
  # Modify the output INCLUDES variable, if it is not set.
  if test ! "${$2+set}" = set ; then
    test -n "${$1}" && $2="-I${$1}/include"
  fi
  # Modify the output LIBS variable, if it is not set.
  if test ! "${$3+set}" = set ; then
    # Save LIBS for later recovery.
    save_LIBS="$LIBS";
    # Extra libraries, if any, required by this check.
    extra_libs="$libz_LIBS -lm"
    # If path is given, add path to extra flag for library search.
    test -n "${$1}" && extra_libs="-L${$1}/lib $extra_libs"
    # Look for library.
    AC_SEARCH_LIBS([xerbla_],lapack,[
      CASC_AC_LOG_VAR(LIBS,After finding lapack flag)
      # Action if found ...
      # Extract modifications to LIB into library-specific LIBS variable.
      $3=`echo " $LIBS" | sed "s! $save_LIBS!!"`;
      test -n "${$1}" && $3="-L${$1}/lib ${$3}"
      CASC_AC_LOG_VAR($3, Found lapack library flag)
      ],[
      # Action if NOT found ...
      CASC_AC_LOG_VAR($3, Did not find lapack library flag)
      AC_MSG_WARN(
[I could not systematically find the name of
the lapack library so I am using -llapack instead.])
      $3="-llapack"
      test -n "${$1}" &&	\
	$3="-L${$1}/lib ${$3}"	# Add path flag to output variable.
      ],[$extra_libs])
    LIBS="$save_LIBS";	# Restore global-use variable.
    unset extra_libs
    unset save_LIBS
  else
    CASC_AC_LOG(Not looking for lapack because $3 is already set)
  fi
fi
])dnl




AC_DEFUN([CASC_SUPPORT_LAPACK],[
dnl Support lapack library by setting the variables
dnl lapack_PREFIX, lapack_INCLUDES, and lapack_LIBS.
dnl Arg1: non-empty if you want the default to be on.
dnl
# Begin macro CASC_SUPPORT_LAPACK

CASC_ARG_WITH_ENV_WRAPPER(lapack, lapack_PREFIX,
ifelse($1,,
[  --with-lapack[=PATH]
			Use lapack and optionally specify where
			they are installed.],
[  --without-lapack	Do not use the lapack library.]),
if test "${with_lapack+set}" = set; then
  lapack_PREFIX=
else
ifelse($1,,unset lapack_PREFIX,lapack_PREFIX=)
fi
)

CASC_ARG_WITH_PREFIX(lapack-includes,lapack_INCLUDES,
[  --with-lapack-includes=STRING
			Specify the INCLUDES flags for lapack.
			If not specified, and --with-lapack=PATH is,
			this defaults to "-IPATH/include".])dnl

CASC_ARG_WITH_PREFIX(lapack-libs,lapack_LIBS,
[  --with-lapack-libs=STRING
			Specify LIBS flags for lapack.
			If not specified, and --with-lapack=PATH is,
			this defaults to "-LPATH/lib -llapack".])dnl

CASC_VAR_SET_LAPACK(lapack_PREFIX,lapack_INCLUDES,lapack_LIBS)

CASC_AC_LOG_VAR(lapack_PREFIX lapack_INCLUDES lapack_LIBS)
# End macro CASC_SUPPORT_LAPACK
])
