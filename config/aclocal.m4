# generated automatically by aclocal 1.11.1 -*- Autoconf -*-

# Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004,
# 2005, 2006, 2007, 2008, 2009  Free Software Foundation, Inc.
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

dnl File arg-with-environment.m4

dnl IMPORTANT NOTE: This macro was originally written to
dnl let configure macros check environments so that developers
dnl can set up "make distcheck" to activate or deactivate
dnl certain packages.  This is largely not needed anymore
dnl because recent automake versions (1.5+, maybe?) allows
dnl you to specify configure options for the distcheck rule.
dnl
dnl Therefore, you are encouraged to use the plain autoconf
dnl macros (AC_ARC_WITH and AC_ARG_ENABLE).

AC_DEFUN([CASC_ARG_WITH_ENV_WRAPPER],[
dnl This is a high-level macro similar to AC_ARG_WITH but it does
dnl   a few extra things.
dnl
dnl It is meant for setting a shell variable using either the
dnl   --with-feature configure flag or by setting a shell variable
dnl   in the environment.  But its primary goal it to set or unset
dnl   the shell variable (arg2).
dnl
dnl One of several things can happen to the shell variable
dnl   when you use this macro, depending first on the configure
dnl   option issued:
dnl   |
dnl   `-- no option given
dnl   |   `-- leave shell variable alone, regardless of whether
dnl   |       it is set (This is how you avoid
dnl   |       having to use the configure option, such as in
dnl   |       running the check rule by automake.)
dnl   `-- with-feature=no or without-feature
dnl   |   `-- unset shell variable, regardless of whether it is set
dnl   `-- with-feature=string
dnl   |   `-- set shell variable to the string
dnl   `-- with-feature or with-feature=yes
dnl       `-- if shell variable already set
dnl       |   `-- leave it a lone
dnl       `-- else if developer gave optional arg4
dnl       |   `-- execute optional arg4 to set shell variable
dnl       `-- else
dnl           `-- set shell variable to blank
dnl
dnl One of two things can happen to the with_feature variable,
dnl   assuming the developer does not change it using arg4.
dnl   `-- no option given
dnl       `-- with_feature is unset
dnl   `-- one of the options referring to "feature" is given
dnl       `-- with_feature is set
dnl
dnl In addition to running AC_ARG_WITH and caching the result, it:
dnl   Allows the variable to be set by the environment.  This is
dnl     for avoiding having to manually issue configure options
dnl     or when manual configure options are not permissible, as
dnl     in running "make check".  The environment variable is
dnl     checked if the --with-something=something_else option
dnl     is not given or given without the equal sign.
dnl   Lets you specify command to run if --with-blah is issued
dnl     without the equal sign or not issued at all.  In this
dnl     case, the environment variable is consulted.  An unset
dnl     environment variable is the same as --without-bla.  A set
dnl     variable is the same as --with-blah=$value.  If $value is an
dnl     empty string, runs optional command (arg5) to set it.
dnl   Lets you specify command (arg4) to check the value chosen
dnl     to make it is good, before caching it.
dnl The arguments to this macro are:
dnl   1: Name of what is being sought (the first argument in
dnl     AC_ARG_WITH).
dnl   2: Name of variable to set (also name of environment variable
dnl      to look for if configure option is not issued).
dnl   3(optional): Help message.  A generic message is provided if
dnl     this argument is empty.
dnl   4(optional): Commands to run if configure flag is not specific
dnl     and environment variable is not set.  These commands are
dnl     run if with_blah is "yes" or "".  They should set or unset
dnl     the variable named in arg2, depending on what you want
dnl     the default behavior to be in these cases.
dnl   5(optional): Quality checking commands, to check if arg2 is good.
dnl     This is run before caching result.  Generally, this would issue
dnl     a warning or error as appropriate.  For example, if this macro
dnl     is used to set the path to a program, you may want to check
dnl     if that program exist and is executable.
# Start macro CASC_ARG_WITH_ENV_WRAPPER with args $1 and $2
AC_CACHE_CHECK(for $1,_cv_prog_[]translit($1,-,_),[
AC_ARG_WITH($1,
ifelse($3,,[  --with-$1	Specify $1 (same as setting $2 in environment)],[$3]))
# Set $2, using environment setting if available
#   and if command line is ambiguous.
case "$with_[]translit($1,-,_)" in
  no[)]
    # User explictly turned off $1.
    # Ignore value of $2, even if set in the environment.
    unset $2
    ;;
  yes|''[)]
    # Flag unissued or ambiguously issued using --with-$1.
    # Because the user did not explicitly turn if off,
    # try to set $2.
    # If environment variable $2 is available, use it.
    # If not, try the user-supplied commands to set it.
    if test -n "${$2}" ;  then
      : Nothing to do here actually, because $2 is already in the environment.
    else
      ifelse($4,,:,$4)
    fi
    ;;
  *)
    # User issued a specific string using --with-$1=non-null-string.
    # so that is used to set $2.
    $2=$with_[]translit($1,-,_)
    ;;
esac
dnl if test ! "${$2+set}" = set ; then
dnl   # $2 is still unset, after processing --with-$1 flag,
dnl   # and possibly using optional command to set it.
dnl   # At this point, check to make sure it is not required.
dnl   # if it is, then we have an error.
dnl   case "$with_[]translit($1,-,_)" in
dnl     no|'')
dnl       : $2 is not set but it is ok because user did not
dnl       : explicitly ask for it by issuing --with-$1=something.
dnl       ;;
dnl     *)
dnl       # If user explicitly asked for $1 and we cannot find it[,]
dnl       # that is an error
dnl       AC_MSG_ERROR(cannot find appropriate value for $2)
dnl       ;;
dnl   esac
dnl fi
if test "${$2+set}" = set ; then
  # This block executes the quality check commands, if any, for $2.
  ifelse($5,,:,$5)
fi
# Cache the value if it was found.
if test "${$2+set}" = set ; then
  _cv_prog_[]translit($1,-,_)=${$2}
fi
])
dnl This part runs if $2 should be set from cache.
# Set $2 from cache.
# $2 is not yet set if we grabbed it from cache.
if test "${_cv_prog_[]translit($1,-,_)+set}" = set ; then
  $2=$_cv_prog_[]translit($1,-,_)
else
  unset $2
fi
# End macro CASC_ARG_WITH_ENV_WRAPPER with args $2 and $1
])




AC_DEFUN([CASC_PATH_PROG],[
dnl This is a high-level macro to find paths to certain programs.
dnl In addition to (possibly) running AC_ARG_WITH and AC_PATH_PROG it:
dnl   Allows the path to be set in an environment variable ($1),
dnl     useful for setting configuration during "make check" and
dnl     for avoiding manual configure options setting.
dnl   Makes sure that the program is executable, if the user explicitly
dnl     specified it.
dnl The arguments are (similar to AC_PATH_PROG):
dnl   1: Variable name to set to the path (used in CASC_PATH_PROG).
dnl   2: Name of program being sought (used in CASC_PATH_PROG).
dnl   3(optional): Commands to set $1 in case neither environment
dnl      nor command line options are given.  Defaults to a call to
dnl      AC_PATH_PROG($1,$2).
dnl   4(optional): Quality check commands to make sure that a
dnl      sufficiently good program is found.  Defaults to a simple
dnl      check that the program is executable.
CASC_ARG_WITH_ENV_WRAPPER($2,$1,
[[  --with-$2=PATH	Specify path to $2 program
			(equivalent to setting $1 in environment)]]dnl
,
[
dnl Commands to run if user was not specific.
# Just set the variable to blank and check later.
$1=
],
dnl Quality check commands.
ifelse($4,,[
  # if $1 is an absolute path, make sure it is executable.
  if echo "${$1}" | grep '^/' > /dev/null && test ! -x "${$1}"; then
    AC_MSG_WARN($2 program ${$1} is not executable.)
  fi],$4)
)dnl
if test "${$1+set}" = set; then
  ifelse($3,,[AC_PATH_PROG($1,$2)],$3)
fi
])




AC_DEFUN([CASC_ARG_WITH_PREFIX],[
dnl This is a high-level macro to set the prefix for a
dnl previously installed package.
dnl The macro arguments are:
dnl 1. package name
dnl 2. variable to contain installation prefix.
dnl 3(optional): Help message.  A generic message is provided if
dnl   this argument is empty.
dnl 4(optional): Commands to run if configure flag is not specific
dnl   and environment variable is not set.  These commands are
dnl   run if with_blah is "yes" or "".  They should set or unset
dnl   the variable named in arg2, depending on what you want
dnl   the defaul behavior to be in these cases.  The default is
dnl   to exit with an error.
# Start macro CASC_ARG_WITH_PREFIX
CASC_ARG_WITH_ENV_WRAPPER($1,$2,
ifelse([$3],,
[[  --with-$1=PATH	Specify prefix where $1 is installed
			(equivalent to setting $2 in the environment)]]
,[[[$3]]]),
ifelse([$4],,
[[if test "${with_[]translit($1,-,_)}" = yes ; then
  AC_MSG_ERROR([[If you specify --with-$1, you must give it the path as in --with-$1=/installation/path]])
fi
CASC_AC_LOG(environment $2 not defined)
]],[[[$4]]])
)dnl
# End macro CASC_ARG_WITH_PREFIX
])




dnl Define a macro for supporting HDF5

AC_DEFUN([CASC_SUPPORT_HDF5],[

# Begin CASC_SUPPORT_HDF5
# Defines hdf5_PREFIX hdf5_INCLUDES and hdf5_LIBS if with-hdf5 is specified.
AC_ARG_WITH(hdf5,
[  --with-hdf5[=PATH]  Use HDF5 and optionally specify where HDF5 is installed.],
, with_hdf5=no)

case "$with_hdf5" in
  no)
    AC_MSG_NOTICE([configuring without HDF5 support])
    : Do nothing
  ;;
  yes)
    # HDF5 install path was not specified.
    # Look in a couple of standard locations to probe if 
    # HDF5 header files are there.
    AC_MSG_CHECKING([for HDF5 installation])
    for dir in /usr /usr/local; do
      if test -f ${dir}/include/hdf5.h; then
        hdf5_PREFIX=${dir}
        break
      fi
    done
    AC_MSG_RESULT([$hdf5_PREFIX])
  ;;
  *)
    # HDF5 install path was specified.
    AC_MSG_CHECKING([for HDF5 installation])

    if test -f ${with_hdf5}/include/hdf5.h; then
        hdf5_PREFIX=$with_hdf5
        hdf5_INCLUDES="-I${hdf5_PREFIX}/include"
        hdf5_LIBS="-L${hdf5_PREFIX}/lib -lhdf5"
        AC_MSG_RESULT([$hdf5_PREFIX])
    else
        AC_MSG_RESULT([$hdf5_PREFIX])
        AC_MSG_ERROR([HDF5 not found in $with_hdf5])
    fi
  ;;
esac



# Test compiling an HDF application

# NOTE that AC_SEARCH_LIBS didn't work completely so
# use a more complicated example program to see
# if that will catch when HDF is not working.
if test "${hdf5_PREFIX+set}" = set; then

   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether HDF5 link works)
   AC_LANG_PUSH(C++)
   CASC_PUSH_COMPILER_STATE
   LIBS="${LIBS} ${hdf5_LIBS} $zlib_LIBS $szlib_LIBS -lm "
   CXXFLAGS="${CXXFLAGS} ${hdf5_INCLUDES}"
   AC_LINK_IFELSE([
      #include "hdf5.h"
      #define FILE "file.h5"

      int main() {

         hid_t       file_id;   /* file identifier */
         herr_t      status;

         /* Create a new file using default properties. */
         file_id = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

         /* Terminate access to the file. */
         status = H5Fclose(file_id); 
     }
      ], 
      casc_hdf5_compile=yes,
      casc_hdf5_compile=no)
   CASC_POP_COMPILER_STATE
   AC_LANG_POP
   AC_MSG_RESULT($casc_hdf5_compile)

   if test "$casc_hdf5_compile" = no; then
      AC_MSG_ERROR([HDF5 compile/link test failed])
   fi
fi

# END CASC_SUPPORT_HDF5

])dnl End definition of CASC_SUPPORT_HDF5




dnl *********************************************************************
dnl * CASC_ADD_LIB(LIBRARY, FUNCTION, DIRECTORY-LIST[, PREFIX[, 
dnl *              ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl * checks first if LIBRARY is available on the linking search path and
dnl * if FUNCTION can be linked with LIBRARY.  If so, -lLIBRARY is added
dnl * to the variable [PREFIX]LIBS. (i.e., if prefix is LD, -llibrary is
dnl * added to LDLIBS.)  If not, checks whitespace-separated
dnl * DIRECTORY-LIST to see if LIBRARY exists in a specified directory and
dnl * can be linked with FUNCTION.  If so, the first directory where
dnl * linking is successful is added to the front of [PREFIX]LIBDIRS, and
dnl * -lLIBRARY is added to the end of [PREFIX]LIBS.  If no prefix is
dnl * specified, the directories and libraries are added to LIBS and
dnl * LIBDIRS, respectively.  If the order of -l flags on the linking
dnl * lines is important, CASC_ADD_LIB should be called for each library
dnl * in the order they should appear on linking lines.  Mere existence of
dnl * LIBRARY in the search path or in a specified directory can usually
dnl * be determined by entering 'main' for FUNCTION.  Optional argument
dnl * ACTION-IF-FOUND contains additional instructions to execute as soon
dnl * as LIBRARY is found in any directory.  Optional argument
dnl * ACTION-IF-NOT-FOUND contains instructions to execute if LIBRARY is
dnl * not found anywhere.
dnl **********************************************************************

AC_DEFUN([CASC_ADD_LIB],
[
   # define some macros to hopefully improve readability
   define([m_THESE_LIBS],[$4LIBS])
   define([m_THESE_LIBDIRS],[$4LIBDIRS])

   # check for the library from first argument.  If linking is successful
   # the first time, the job is done, otherwise loop through DIRECTORY-LIST
   AC_CHECK_LIB($1, $2, m_THESE_LIBS="$m_THESE_LIBS -l$1"
                          casc_lib_found=yes 
                          ifelse([$5], , , [$5]),

      dnl * If library not found
      for casc_lib_dir in $3; do

         AC_CHECK_LIB($1, $2, 
            m_THESE_LIBDIRS="-L$casc_lib_dir $m_THESE_LIBDIRS"
            m_THESE_LIBS="$m_THESE_LIBS -l$1"
            casc_lib_found=yes
            ifelse([$5], , , [$5])
            break
            , ,
            -L$casc_lib_dir $m_THESE_LIBDIRS $m_THESE_LIBS -l$1, no)
      done
      , $m_THESE_LIBDIRS $m_THESE_LIBS, no)  dnl * last two arguments for
                                             dnl * first check

   # ACTION-IF-NOT_FOUND for when the library is found nowhere
   ifelse([$6], , ,
      if test "$casc_lib_found" != "yes"; then
         [$6]
      fi
   )

   unset casc_lib_found

   undefine([m_THESE_LIBS])
   undefine([m_THESE_LIBDIRS])

])dnl


dnl ***********************************************************************
dnl CASC_CHECK_LIB(LIBRARY, FUNCTION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND
dnl              [, OTHER-LIBRARIES [, CACHE-CHOICE]]]])
dnl * This is the same as AC_CHECK_LIB, except when it tests for LIBRARY
dnl * it puts the flag -lLIBRARY after $LIBS and OTHER-LIBRARIES.  The Sun
dnl * cc compiler does not search for LIBRARY in any directories specified
dnl * by -L in OTHER-LIBRARIES when -lLIBRARY is listed first.  The
dnl * functionality of this macro is the same as that of AC_CHECK_LIB in
dnl * the Autoconf documentation.  
dnl * CACHE-CHOICE [$6]added by N. Elliott, 6-24-98.  If CACHE-CHOICE is 'no',
dnl * the results of this test are not cached.  CACHE-CHOICE should be
dnl * used only when this test is called recursively.
dnl *
dnl * CASC_CHECK_LIB_OLD is an older version of this macro which doesn't
dnl * seem to work with newer versions of autoconf
dnl **********************************************************************

AC_DEFUN([CASC_CHECK_LIB],
[
dnl AC_MSG_CHECKING([for $2 in -l$1])
dnl Use a cache variable name containing both the library and function name,
dnl because the test really is for library $1 defining function $2, not
dnl just for library $1.  Separate tests with the same $1 and different $2s
dnl may have different results.
ac_lib_var=`echo $1['_']$2 | sed 'y%./+-%__p_%'`
AC_CACHE_VAL(ac_cv_lib_$ac_lib_var,
[ac_save_LIBS="$LIBS"
LIBS="$5 $LIBS -l$1"
AC_TRY_LINK(dnl
ifelse(AC_LANG, [FORTRAN77], ,
ifelse([$2], [main], , dnl Avoid conflicting decl of main.
[/* Override any gcc2 internal prototype to avoid an error.  */
]ifelse(AC_LANG, CPLUSPLUS, [#ifdef __cplusplus
extern "C"
#endif
])dnl
[/* We use char because int might match the return type of a gcc2
    builtin and then its argument prototype would still apply.  */
char $2();
])),
            [$2()],
            eval "ac_cv_lib_$ac_lib_var=yes",
            eval "ac_cv_lib_$ac_lib_var=no")   
LIBS="$ac_save_LIBS"
])dnl 
if eval "test \"`echo '$ac_cv_lib_'$ac_lib_var`\" = yes"; then
  :
  dnl AC_MSG_RESULT(yes)
  ifelse([$3], ,
[changequote(, )dnl
  ac_tr_lib=HAVE_LIB`echo $1 | sed -e 's/[^a-zA-Z0-9_]/_/g' \
    -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
changequote([, ])dnl
  AC_DEFINE_UNQUOTED($ac_tr_lib)
  LIBS="-l$1 $LIBS" 
], [$3])
else
 : 
  dnl AC_MSG_RESULT(no)
ifelse([$4], , , [
ifelse([$6], no, unset ac_cv_lib_$ac_lib_var)  
$4
])dnl
fi
ifelse([$6], no, unset ac_cv_lib_$ac_lib_var)
])dnl



AC_DEFUN([CASC_CHECK_LIB_OLD],
[AC_MSG_CHECKING([for -l$1])
dnl Use a cache variable name containing both the library and function name,
dnl because the test really is for library $1 defining function $2, not
dnl just for library $1.  Separate tests with the same $1 and different $2s
dnl may have different results.
ac_lib_var=`echo $1['_']$2 | tr './+\055' '__p_'`
AC_CACHE_VAL(ac_cv_lib_$ac_lib_var,
[ac_save_LIBS="$LIBS"
LIBS="$5 $LIBS -l$1"
AC_TRY_LINK(dnl
ifelse([$2], [main], , dnl Avoid conflicting decl of main.
[/* Override any gcc2 internal prototype to avoid an error.  */
]ifelse(AC_LANG, CPLUSPLUS, [#ifdef __cplusplus 
extern "C"
#endif
])dnl
[/* We use char because int might match the return type of a gcc2
    builtin and then its argument prototype would still apply.  */
char $2();
]),
            [$2()],
            eval "ac_cv_lib_$ac_lib_var=yes",
            eval "ac_cv_lib_$ac_lib_var=no")dnl
LIBS="$ac_save_LIBS"
])dnl
if eval "test \"`echo '$ac_cv_lib_'$ac_lib_var`\" = yes"; then
  AC_MSG_RESULT(yes)  
  ifelse([$3], ,
[changequote(, )dnl
  ac_tr_lib=HAVE_LIB`echo $1 | tr 'abcdefghijklmnopqrstuvwxyz' 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'`
changequote([, ])dnl
  AC_DEFINE_UNQUOTED($ac_tr_lib)
  LIBS="-l$1 $LIBS"
], [
ifelse([$6], no, unset ac_cv_lib_$ac_lib_var)
$3])
else
  AC_MSG_RESULT(no) 
ifelse([$4], , , [
ifelse([$6], no, unset ac_cv_lib_$ac_lib_var)
$4
])dnl
fi
ifelse([$6], no, unset ac_cv_lib_$ac_lib_var)
])



dnl *********************************************************************
dnl * CASC_CHECK_HEADER(HEADER-FILE, DIRECTORY-LIST[, ACTION-IF-FOUND[,
dnl *                   ACTION-IF-NOT-FOUND]])
dnl * This macro is an alternative to AC_CHECK_HEADER.  It does
dnl * essentially the same thing, but it allows the user to specify
dnl * a directory list if HEADER-FILE can not be found in the current path
dnl * for #includes, and it adds to the variable INCLUDES the first
dnl * directory in DIRECTORY-LIST from where HEADER-FILE can be included.
dnl *********************************************************************

AC_DEFUN([CASC_CHECK_HEADER],
[
   dnl * loop through the directory list.  The first iteration leaves the
   dnl * casc_dir variable empty to check if the header can be #included
   dnl * without specifying a directory.
   for casc_dir in '' $2 ; do
      if test -n "$casc_dir"; then
         casc_header=$casc_dir/$1
      else
         casc_header=$1
      fi

      dnl * Check for the header.  Add the necessary -I flag to INCLUDES
      AC_CHECK_HEADER( $casc_header,
         if test -n "$casc_dir"; then
            INCLUDES="$INCLUDES -I$casc_dir"
         fi
         casc_header_found=yes
         ifelse([$3], , , [$3])
         break )

   done

   dnl * This takes care of the action if not found
   ifelse([$4], , ,
      if test "$casc_header_found" != "yes"; then
         [$4]
      fi
   )

   unset casc_header_found
])dnl


dnl **********************************************************************
dnl * CASC_CREATE_PACKAGE_OPTION(PACKAGE-NAME[, DIR-LIST[, FILE]])
dnl * This is a general macro that creates a configure command-line option
dnl * called `--with-PACKAGE-NAME-dir' which will allow the user to
dnl * specify the location of the installation of an outside software
dnl * package, such as PETSc or ISIS++.  After a check to make sure the
dnl * given directory is valid (see below for discussion of validity), the
dnl * directory's path is stored in the shell variable PACKAGE-NAME_DIR.
dnl * For example, to allow the user to specify the location of PETSc,
dnl * place `CASC_CREATE_PACKAGE_OPTION(PETSC)' in configure.in.  Then the
dnl * user, if configuring on the CASC Sun cluster, would type `configure
dnl * --with-PETSC-dir=/home/casc/petsc', and the directory's path would
dnl * be stored in PETSC_DIR.  With this macro, the user is also permitted
dnl * to set the variable PACKAGE-NAME_DIR in the environment before
dnl * running configure, but any choice made on the command line would
dnl * override any preset values.  
dnl *
dnl * This macro takes an optional second argument, DIR-LIST, which is a
dnl * whitespace-separated list of directories where the developer thinks
dnl * PACKAGE-NAME might be installed.  If DIR-LIST is given, and the user
dnl * does not use the `--with' option to give the location of
dnl * PACKAGE-NAME (or if the directory given by the user does not exist),
dnl * then configure will assign to PACKAGE-NAME_DIR the path of the first
dnl * directory in DIR-LIST that is valid.
dnl *
dnl * Validity:  The optional third argument to this macro is FILE, which
dnl * should be either the name of a file in the top directory of the
dnl * package in question or the relative path of a file in a subdirectory
dnl * of the package.  If the argument FILE is given, then configure will
dnl * consider a user specified directory or a directory from DIR-LIST 
dnl * valid only if FILE exists in the directory.  If this argument is not
dnl * given, then configure will consider a directory valid simply if it
dnl * is indeed a directory.  FILE should be a file with a unique name
dnl * that can be expected to exist in the same location in any 
dnl * installation of the package in question.  If you know of no such
dnl * file, do not include a third argument when invoking this macro.
dnl * 
dnl * This macro also gives the user the command-line option
dnl * `--without-PACKAGE-NAME-dir', which, when invoked, will leave the
dnl * variable PACKAGE-NAME_DIR empty.  This option should be invoked when
dnl * the user wants to exclude a package from the configuration.
dnl * 
dnl * NOTE:  Since PACKAGE-NAME is used as part of both a command-line
dnl * option and a variable name, it MUST consist of only alphanumeric
dnl * characters.  PACKAGE-NAME is only a label, so it need not conform to
dnl * any existing directory or file name.  I would recommend that it be
dnl * all caps, as it becomes part of the name of a variable that is
dnl * substituted into the Makefile.
dnl **********************************************************************

AC_DEFUN([CASC_CREATE_PACKAGE_OPTION],
[
   AC_MSG_CHECKING([for $1 directory])

   dnl * $1 stands for the PACKAGE-NAME.  If [$1]_DIR has been set in the
   dnl * environment, give its value to casc_env_[$1]_dir, and clear
   dnl * [$1]_DIR.  The environmental value will ultimately be reassigned
   dnl * to [$1]_DIR if it is valid and no command-line options are able
   dnl * to change [$1]_DIR to a valid directory.  The environmental value
   dnl * will also be used even if it is invalid, if the command-line
   dnl * options and the DIRECTORY-LIST are both unable to generate a
   dnl * valid value.
   casc_result=
   casc_env_[$1]_dir=$[$1]_DIR
   [$1]_DIR=

   AC_ARG_WITH($1-dir, 
[  --with-$1-dir=DIR    $1 is installed in directory DIR
  --without-$1-dir     do not look for $1],

               if test "$withval" = "no"; then
                  casc_result="configuring without [$1]"
                  [$1]_DIR=
               fi
               , )

   dnl * If "--without-$1-dir" was given, then [$1]_DIR is left blank.
   dnl * Otherwise there is the following procedure to try to give
   dnl * [$1]_DIR a valid value:
   dnl *
   dnl * if "--with-$1-dir" was given
   dnl *    if the argument to "--with-$1-dir" is valid
   dnl *       assign the argument to [$1]_DIR
   dnl *    endif
   dnl * endif
   dnl *
   dnl * if a value for [$1]_DIR has not yet been found
   dnl *    if [$1]_DIR from the environment exists and is valid
   dnl *       assign the environmental value to [$1]_DIR
   dnl *    endif
   dnl * endif
   dnl *
   dnl * if [$1]_DIR still has no value
   dnl *    if the macro was given a DIRECTORY-LIST argument
   dnl *       for each directory in the list
   dnl *          if the directory is valid
   dnl *             assign the directory to [$1]_DIR
   dnl *             break loop
   dnl *          else
   dnl *             continue loop
   dnl *          endif
   dnl *       end loop
   dnl *       if [$1]_DIR still doesn't have a value
   dnl *          casc_result="none"
   dnl *       else
   dnl *          casc_result=$[$1]_DIR
   dnl *       endif
   dnl *    else
   dnl *       casc_result="none"
   dnl *    endif
   dnl * endif

   if test "$with_[$1]_dir" != "no"; then

      if test -d "$with_[$1]_dir"; then

         ifelse([$3], , ,
            if test -f $with_[$1]_dir/[$3]; then)

               casc_result="$with_[$1]_dir"
               [$1]_DIR="$casc_result"

         ifelse([$3], , ,
            fi)
      fi

      if test -z "$casc_result"; then

         if test -d "$casc_env_[$1]_dir"; then

            ifelse([$3], , ,
               if test -f $casc_env_[$1]_dir/[$3]; then)

                  casc_result="$casc_env_[$1]_dir"
                  [$1]_DIR="$casc_result"

            ifelse([$3], , ,
               fi)
         fi
      fi



      if test -z "$casc_result"; then
         [$1]_DIR=
   
         ifelse([$2], ,
            casc_result="none" ,

            for casc_dir in $2; do

               if test -d "$casc_dir"; then

                  ifelse([$3], , ,
                     if test -f $casc_dir/[$3]; then)

                        $1_DIR=$casc_dir

                        break

                  ifelse([$3], , ,
                     fi)
               fi
            done

            if test -z "$[$1]_DIR"; then
               casc_result="none"

            else
               casc_result="$[$1]_DIR"
            fi
         )
      fi
   fi

   dnl * $casc_result either is a valid value for [$1]_DIR or "none".
   dnl * if none, then assign the original environmental value of
   dnl * [$1]_DIR, whatever it may be, to casc_result and [$1]_DIR.  If
   dnl * there was no environmental value, then $casc_result remains
   dnl * "none" and [$1]_DIR is left empty.

   if test "$casc_result" = "none"; then

      if test -n "$casc_env_[$1]_dir"; then

         casc_result="$casc_env_[$1]_dir"
         [$1]_DIR="$casc_result"
      fi
   fi

   AC_MSG_RESULT($casc_result)
   AC_SUBST([$1]_DIR)
])


dnl smr_ARG_WITHLIB from FVWM by S. Robbins 
dnl Allow argument for optional libraries; wraps AC_ARG_WITH, to
dnl provide a "--with-foo-lib" option in the configure script, where foo
dnl is presumed to be a library name.  The argument given by the user
dnl (i.e. "bar" in ./configure --with-foo-lib=bar) may be one of four 
dnl things:
dnl     * boolean (no, yes or blank): whether to use library or not
dnl     * file: assumed to be the name of the library
dnl     * directory: assumed to *contain* the library
dnl     * a quoted, space-separated list of linker flags needed to link
dnl       with this library.  (To be used if this library requires
dnl       linker flags other than the normal `-L' and `-l' flags.)
dnl 
dnl The argument is sanity-checked.  If all is well, two variables are
dnl set: "with_foo" (value is yes, no, or maybe), and "foo_LIBFLAGS" (value
dnl is either blank, a file, -lfoo, '-L/some/dir -lfoo', or whatever 
dnl linker flags the user gives). The idea is: the first tells you whether
dnl the library is to be used or not (or the user didn't specify one way
dnl or the other) and the second to put on the command line for linking
dnl with the library.
dnl
dnl Usage:
dnl smr_ARG_WITHLIB(name, libname, description)
dnl 
dnl name                name for --with argument ("foo" for libfoo)
dnl libname             (optional) actual name of library,
dnl                     if different from name
dnl description         (optional) used to construct help string
dnl 
dnl Changes:  Changed some identifier names.
dnl           --with-foo-library is now --with-foo-lib
dnl           foo_LIBS is now foo_LIBFLAGS
dnl           Fourth posibility for argument to --with-foo-lib added
dnl           Documentation above changed to reflect these changes
dnl           Noah Elliott, October 1998


AC_DEFUN([CASC_SMR_ARG_WITHLIB],
[
   smr_ARG_WITHLIB([$1],[$2],[$3])
])dnl

AC_DEFUN([smr_ARG_WITHLIB], [

ifelse($2, , smr_lib=[$1], smr_lib=[$2]) 
    
AC_ARG_WITH([$1]-lib,
ifelse($3, ,
[  --with-$1-lib[=PATH]       use $1 library], 
[  --with-$1-lib[=PATH]       use $1 library ($3)]),
[
    if test "$withval" = yes; then
        with_[$1]=yes
        [$1]_LIBFLAGS="-l${smr_lib}"
    elif test "$withval" = no; then
        with_[$1]=no
        [$1]_LIBFLAGS=
    else
        with_[$1]=yes
        if test -f "$withval"; then
            [$1]_LIBFLAGS=$withval
        elif test -d "$withval"; then
            [$1]_LIBFLAGS="-L$withval -l${smr_lib}"
        else
            case $withval in
            -*)
               [$1]_LIBFLAGS="$withval"
            ;;
            *)
               AC_MSG_ERROR(
                  [argument must be boolean, file, directory, or compiler flags]
                           )
            ;;
            esac
        fi
    fi
], [
    with_[$1]=maybe
    [$1]_LIBFLAGS="-l${smr_lib}"
])])

    
dnl smr_ARG_WITHINCLUDES from FVWM by S. Robbins
dnl Check if the include files for a library are accessible, and
dnl define the variable "name_INCLUDE" with the proper "-I" flag for
dnl the compiler.  The user has a chance to specify the includes
dnl location, using "--with-foo-include".
dnl 
dnl This should be used *after* smr_ARG_WITHLIB *and* AC_CHECK_LIB are
dnl successful.
dnl 
dnl Usage:
dnl smr_ARG_WITHINCLUDES(name, header, extra-flags)
dnl 
dnl name                library name, MUST same as used with smr_ARG_WITHLIB
dnl header              a header file required for using the lib
dnl extra-flags         (optional) flags required when compiling the
dnl                     header, typically more includes; for ex. X_CFLAGS
dnl
dnl Changes:  Changed some identifier names.
dnl           --with-foo-includes is now --with-foo-include
dnl           name_CFLAGS is now name_INCLUDE
dnl           Documentation above changed to reflect these changes
dnl           Noah Elliott, October 1998

AC_DEFUN([CASC_SMR_ARG_WITHINCLUDES],
[
   smr_ARG_WITHINCLUDES([$1], [$2], [$3])
])dnl

AC_DEFUN([smr_ARG_WITHINCLUDES], [

AC_ARG_WITH([$1]-include,
[  --with-$1-include=DIR  set directory for $1 headers],
[
    if test -d "$withval"; then
        [$1]_INCLUDE="-I${withval}"
    else
        AC_MSG_ERROR(argument must be a directory)
    fi])

dnl This bit of logic comes from autoconf's AC_PROG_CC macro.  We need
dnl to put the given include directory into CPPFLAGS temporarily, but
dnl then restore CPPFLAGS to its old value.
dnl 
smr_test_CPPFLAGS="${CPPFLAGS+set}"
smr_save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS ${[$1]_CFLAGS}"

    ifelse($3, , , CPPFLAGS="$CPPFLAGS [$3]")
    AC_CHECK_HEADERS($2)
   
if test "$smr_test_CPPFLAGS" = set; then
    CPPFLAGS=$smr_save_CPPFLAGS
else
    unset CPPFLAGS
fi
])
    
        
dnl smr_CHECK_LIB from FVWM by S. Robbins
dnl Probe for an optional library.  This macro creates both
dnl --with-foo-lib and --with-foo-include options for the configure
dnl script.  If --with-foo-lib is *not* specified, the default is to
dnl probe for the library, and use it if found.
dnl
dnl Usage:
dnl smr_CHECK_LIB(name, libname, desc, func, header, x-libs, x-flags)
dnl 
dnl name        name for --with options
dnl libname     (optional) real name of library, if different from
dnl             above
dnl desc        (optional) short descr. of library, for help string
dnl func        function of library, to probe for
dnl header      (optional) header required for using library
dnl x-libs      (optional) extra libraries, if needed to link with lib
dnl x-flags     (optional) extra flags, if needed to include header files
dnl
dnl Changes:  identifier names and documentation modified to reflect
dnl           changes to smr_ARG_WITHLIB and smr_ARG_WITHINCLUDES
dnl           Noah Elliott, October 1998

AC_DEFUN([CASC_SMR_CHECK_LIB],
[
   smr_CHECK_LIB([$1], [$2], [$3], [$4], [$5], [$6], [$7])
])dnl

AC_DEFUN([smr_CHECK_LIB],
[   
ifelse($2, , smr_lib=[$1], smr_lib=[$2])
ifelse($5, , , smr_header=[$5])
smr_ARG_WITHLIB($1,$2,$3)
if test "$with_$1" != no; then
    AC_CHECK_LIB($smr_lib, $4,
        smr_havelib=yes, smr_havelib=no,
        ifelse($6, , ${$1_LIBFLAGS}, [${$1_LIBFLAGS} $6]))
    if test "$smr_havelib" = yes -a "$smr_header" != ""; then
        smr_ARG_WITHINCLUDES($1, $smr_header, $7)
        smr_safe=`echo "$smr_header" | sed 'y%./+-%__p_%'`
        if eval "test \"`echo '$ac_cv_header_'$smr_safe`\" != yes"; then
            smr_havelib=no
        fi
    fi
    if test "$smr_havelib" = yes; then
        AC_MSG_RESULT(Using $1 library)
    else
        $1_LIBFLAGS=
        $1_INCLUDE=
        if test "$with_$1" = maybe; then
            AC_MSG_RESULT(Not using $1 library)
        else
            AC_MSG_WARN(Requested $1 library not found!)
        fi
    fi
fi])


dnl **********************************************************************
dnl * CASC_CONFIG_OUTPUT_LIST(DIR-LIST[, OUTPUT-FILE])
dnl *
dnl * The intent of this macro is to make configure handle the possibility
dnl * that a portion of the directory tree of a project may not be
dnl * present.  This will modify the argument list of AC_OUTPUT to contain
dnl * only output file names for which corresponding input files exist.
dnl * If you are not concerned about the possible absence of the necessary
dnl * input (.in) files, it is better to not use this macro and to
dnl * explicitly list all of the output files in a call to AC_OUTPUT.
dnl * Also, If you wish to create a file Foo from a file with a name
dnl * other than Foo.in, this macro will not work, and you must use
dnl * AC_OUTPUT.
dnl *
dnl * This macro checks for the existence of the file OUTPUT-FILE.in in
dnl * each directory specified in the whitespace-separated DIR-LIST.  
dnl * (Directories should be specified by relative path from the directory 
dnl * containing configure.in.) If OUTPUT-FILE is not specified, the
dnl * default is 'Makefile'.  For each directory that contains 
dnl * OUTPUT-FILE.in, the relative path of OUTPUT-FILE is added to the 
dnl * shell variable OUTPUT-FILE_list.  When AC_OUTPUT is called,
dnl * '$OUTPUT-FILE_list' should be included in the argument list.  So if
dnl * you have a directory tree and each subdirectory contains a 
dnl * Makefile.in, DIR-LIST should be a list of every subdirectory and
dnl * OUTPUT-FILE can be omitted, because 'Makefile' is the default.  When
dnl * configure runs, it will check for the existence of a Makefile.in in
dnl * each directory in DIR-LIST, and if so, the relative path of each
dnl * intended Makefile will be added to the variable Makefile_list.
dnl *
dnl * This macro can be called multiple times, if there are files other
dnl * than Makefile.in with a .in suffix other that are intended to be 
dnl * processed by configure. 
dnl *
dnl * Example
dnl *     If directories dir1 and dir2 both contain a file named Foo.in, 
dnl *     and you wish to use configure to create a file named Foo in each
dnl *     directory, then call 
dnl *     CASC_CONFIG_OUTPUT_LIST(dir1 dir2, Foo)
dnl *     If you also called this macro for Makefile as described above,
dnl *     you should call
dnl *     AC_OUTPUT($Makefile_list $Foo_list)
dnl *     at the end of configure.in .
dnl *********************************************************************


AC_DEFUN([CASC_CONFIG_OUTPUT_LIST],
[
   dnl * m_OUTPUT_LIST is a macro to store the name of the variable
   dnl * which will contain the list of output files
   define([m_OUTPUT_LIST], ifelse([$2], , Makefile_list, [$2_list]))

   if test -z "$srcdir"; then
      srcdir=.
   fi

   dnl * use "Makefile" if second argument not given
   if test -n "$2"; then
      casc_output_file=$2
   else   
      casc_output_file=Makefile
   fi   
      
   dnl * Add a file to the output list if its ".in" file exists.
   for casc_dir in $1; do
      if test -f $srcdir/$casc_dir/$casc_output_file.in; then
         m_OUTPUT_LIST="$m_OUTPUT_LIST $casc_dir/$casc_output_file"
      fi
   done
])dnl


dnl **********************************************************************
dnl * CASC_GUESS_ARCH
dnl * Guesses a one-word name for the current architecture, unless ARCH
dnl * has been preset.  This is an alternative to the built-in macro
dnl * AC_CANONICAL_HOST, which gives a three-word name.  Uses the utility
dnl * 'tarch', which is a Bourne shell script that should be in the same  
dnl * directory as the configure script.  If tarch is not present or if it
dnl * fails, ARCH is set to the value, if any, of shell variable HOSTTYPE,
dnl * otherwise ARCH is set to "unknown".
dnl **********************************************************************

AC_DEFUN([CASC_GUESS_ARCH],
[
   AC_MSG_CHECKING(the architecture)

   dnl * $ARCH could already be set in the environment or earlier in configure
   dnl * Use the preset value if it exists, otherwise go throug the procedure
   if test -z "$ARCH"; then

      dnl * configure searches for the tool "tarch".  It should be in the
      dnl * same directory as configure.in, but a couple of other places
      dnl * will be checked.  casc_tarch stores a relative path for "tarch".
      casc_tarch_dir=
      for casc_dir in $srcdir $srcdir/.. $srcdir/../.. $srcdir/config; do
         if test -f $casc_dir/tarch; then
            casc_tarch_dir=$casc_dir
            casc_tarch=$casc_tarch_dir/tarch
            break
         fi
      done

      dnl * if tarch was not found or doesn't work, try using env variable
      dnl * $HOSTTYPE
      if test -z "$casc_tarch_dir"; then
         AC_MSG_WARN(cannot find tarch, using \$HOSTTYPE as the architecture)
         ARCH=$HOSTTYPE
      else
         ARCH="`$casc_tarch`"

         if test -z "$ARCH" || test "$ARCH" = "unknown"; then
            ARCH=$HOSTTYPE
         fi
      fi

      dnl * if $ARCH is still empty, give it the value "unknown".
      if test -z "$ARCH"; then
         ARCH=unknown
         AC_MSG_WARN(architecture is unknown)
      else
         AC_MSG_RESULT($ARCH)
      fi    
   else
      AC_MSG_RESULT($ARCH)
   fi

   AC_SUBST(ARCH)

])dnl


dnl **********************************************************************
dnl * CASC_SET_SUFFIX_RULES is not like the other macros in aclocal.m4
dnl * because it does not run any kind of test on the system on which it
dnl * is running.  All it does is create several variables which contain
dnl * the text of some simple implicit suffix rules that can be
dnl * substituted into Makefile.in.  The suffix rules that come from the
dnl * macro all deal with compiling a source file into an object file.  If
dnl * this macro is called in configure.in, then if `@CRULE@' is placed in
dnl * Makefile.in, the following will appear in the generated Makefile:
dnl *
dnl * .c.o:
dnl *         @echo "Making (c) " $@ 
dnl *         @${CC} -o $@ -c ${CFLAGS} $<	
dnl *
dnl * The following is a list of the variables created by this macro and
dnl * the corresponding suffixes of the files that each implicit rule 
dnl * deals with.
dnl *
dnl * CRULE       --   .c
dnl * CXXRULE     --   .cxx
dnl * CPPRULE     --   .cpp
dnl * CCRULE      --   .cc
dnl * CAPCRULE    --   .C
dnl * F77RULE     --   .f
dnl *
dnl * There are four suffix rules for C++ files because of the different
dnl * suffixes that can be used for C++.  Only use the one which
dnl * corresponds to the suffix you use for your C++ files.
dnl *
dnl * The rules created by this macro require you to use the following
dnl * conventions for Makefile variables:
dnl *
dnl * CC        = C compiler
dnl * CXX       = C++ compiler
dnl * F77       = Fortran 77 compiler
dnl * CFLAGS    = C compiler flags
dnl * CXXFLAGS  = C++ compiler flags
dnl * FFLAGS    = Fortran 77 compiler flags
dnl **********************************************************************

AC_DEFUN([CASC_SET_SUFFIX_RULES],
[
   dnl * Things weren't working whenever "$@" showed up in the script, so
   dnl * I made the symbol $at_sign to signify '@'
   at_sign=@

   dnl * All of the backslashes are used to handle the $'s and the
   dnl * newlines which get passed through echo and sed.

   CRULE=`echo ".c.o:\\\\
\t@echo \"Making (c) \" \\$$at_sign \\\\
\t@\\${CC} -o \\$$at_sign -c \\${CFLAGS} \$<"`

   AC_SUBST(CRULE)

   CXXRULE=`echo ".cxx.o:\\\\
\t@echo \"Making (c++) \" \\$$at_sign \\\\
\t@\\${CXX} -o \\$$at_sign -c \\${CXXFLAGS} \$<"`

   AC_SUBST(CXXRULE)

   CPPRULE=`echo ".cpp.o:\\\\
\t@echo \"Making (c++) \" \\$$at_sign \\\\
\t@\\${CXX} -o \\$$at_sign -c \\${CXXFLAGS} \$<"`

   AC_SUBST(CPPRULE)

   CCRULE=`echo ".cc.o:\\\\
\t@echo \"Making (c++) \" \\$$at_sign \\\\
\t@\\${CXX} -o \\$$at_sign -c \\${CXXFLAGS} \$<"`

   AC_SUBST(CCRULE)

   CAPCRULE=`echo ".C.o:\\\\
\t@echo \"Making (c++) \" \\$$at_sign \\\\
\t@\\${CXX} -o \\$$at_sign -c \\${CXXFLAGS} \$<"`

   AC_SUBST(CAPCRULE)

   F77RULE=`echo ".f.o:\\\\
\t@echo \"Making (f) \" \\$$at_sign \\\\
\t@\\${F77} -o \\$$at_sign -c \\${FFLAGS} \$<"`

   AC_SUBST(F77RULE)

])

dnl Macro to save compiler state flags for invoking dnl compiler tests
dnl NOTE that this is NOT currently a stack so can dnl only be called
dnl in push/pop order.  push push pop pop dnl will fail
AC_DEFUN([CASC_PUSH_COMPILER_STATE],[
   casc_save_LIBS=$LIBS
   casc_save_CXXFLAGS=$CXXFLAGS
])

dnl Macro to restore compiler state flags for invoking
dnl compiler tests
AC_DEFUN([CASC_POP_COMPILER_STATE],[
   LIBS=$casc_save_LIBS
   unset casc_save_LIBS
   CXXFLAGS=$casc_save_CXXFLAGS
   unset casc_save_CXXFLAGS
])

dnl ********************************************************************
dnl * CASC_PROG_MPICC searches the PATH for an available MPI C compiler
dnl * wraparound.  It assigns the name to MPICC.
dnl ********************************************************************

AC_DEFUN([CASC_PROG_MPICC],
[
   AC_CHECK_PROGS(MPICC, mpcc mpicc tmcc hcc)
   test -z "$MPICC" && AC_MSG_ERROR([no acceptable mpicc found in \$PATH])
])dnl


dnl ********************************************************************
dnl * CASC_PROG_MPICXX searches the PATH for an available MPI C++
dnl * compiler wraparound.  It assigns the name to MPICXX.
dnl ********************************************************************

AC_DEFUN([CASC_PROG_MPICXX],
[
   AC_CHECK_PROGS(MPICXX, mpKCC mpCC mpig++ mpiCC hcp)
   test -z "$MPICXX" && AC_MSG_ERROR([no acceptable mpic++ found in \$PATH])
])dnl


dnl **********************************************************************
dnl * CASC_PROG_MPIF77 searches the PATH for an available MPI Fortran 77
dnl * compiler wraparound.  It assigns the name to MPIF77.
dnl **********************************************************************

AC_DEFUN([CASC_PROG_MPIF77],
[
   AC_CHECK_PROGS(MPIF77, mpf77 mpxlf mpif77 mpixlf tmf77 hf77)
   test -z "$MPIF77" && AC_MSG_ERROR([no acceptable mpif77 found in \$PATH])
])dnl


dnl ***********************************************************************
dnl * CASC_CHECK_MPIF77_PP checks whether the preprocessor needs to
dnl * be called before calling the compiler for Fortran files with
dnl * preprocessor directives and MPI function calls.  If the preprocessor
dnl * is necessary, MPIF77NEEDSPP is set to "yes", otherwise it is set to
dnl * "no"
dnl ***********************************************************************

AC_DEFUN([CASC_CHECK_MPIF77_PP],
[
   AC_REQUIRE([CASC_PROG_MPIF77])

   rm -f testppmp.o

   AC_MSG_CHECKING(whether $FPP needs to be called before $MPIF77)

   # This follows the same procedur as CASC_CHECK_F77_PP, except it tests
   # $MPIF77 using a test program that includes MPI functions.

   cat > testppmp.F <<EOF
#define FOO 3
	program testppmp
	include 'mpif.h'
	integer rank,size,mpierr,sum
	call MPI_INIT(mpierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,size,mpierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,rank,mpierr)
#ifdef FORTRAN_NO_UNDERSCORE
        sum = rank + size
#else
        sum = rank + rank
#endif
        call MPI_FINALIZE(mpierr)
        end
EOF

   $MPIF77 -DBAR -c testppmp.F
   if test -f testppmp.o; then
      MPIF77NEEDSPP=no
   else
      MPIF77NEEDSPP=yes
   fi

   echo $MPIF77NEEDSPP
   rm -f testppmp.o testppmp.F
   AC_SUBST(MPIF77NEEDSPP)
])dnl


dnl *********************************************************************
dnl * CASC_SET_MPI sets up the needed MPI library and directory flags.
dnl * The location of the file mpi.h is put into the variable MPIINCLUDE
dnl * as a -I flag.  The -l flags that specify the needed libraries and
dnl * the -L flags that specify the paths of those libraries are placed in
dnl * the variables MPILIBS and MPILIBDIRS, respectively.  To set the MPI
dnl * libraries and directories manually, use the --with-mpi-include,
dnl * --with-mpi-libs, and --with-mpi-lib-dirs command-line options when
dnl * invoking configure.  Only one directory should be specified with
dnl * --with-mpi-include, while any number of directories can be specified
dnl * by --with-mpi-lib-dirs.  Any number of libraries can be specified
dnl * with --with-mpi-libs, and the libraries must be referred to by their
dnl * base names, so libmpi.a is just mpi.  It is adviseable to use all
dnl * three --with flags whenever one is used, because it is likely that
dnl * when one is chosen it will mess up the automatic choices for the
dnl * other two.  If the architecture is unknown, or if the needed MPI
dnl * settings for the current architecture are not known, then the naive
dnl * settings of MPILIBS="-lmpi" and MPILIBDIRS="-L/usr/local/mpi/lib"
dnl * are tested, and if they exist they are used, otherwise the MPILIB*
dnl * variables are left blank.  In the case of rs6000, the variable
dnl * MPIFLAGS is also set.
dnl **********************************************************************

AC_DEFUN([CASC_SET_MPI],
        [

   dnl * If called from within CASC_FIND_MPI, then the configure-line
   dnl * options will already exist.  This ifdef creates them otherwise.
   ifdef([AC_PROVIDE_CASC_FIND_MPI], ,
      [AC_ARG_WITH(mpi-include, [  --with-mpi-include=DIR  mpi.h is in DIR],
                  casc_mpi_include_dir=$withval)

      AC_ARG_WITH(mpi-libs,
[  --with-mpi-libs=LIBS    LIBS is space-separated list of library names
                          needed for MPI, e.g. \"nsl socket mpi\"],
                  casc_mpi_libs=$withval)

      AC_ARG_WITH(mpi-lib-dirs,
[  --with-mpi-lib-dirs=DIRS
                          DIRS is space-separated list of directories
                          containing the libraries specified by
                          \`--with-mpi-libs', e.g \"/usr/lib /usr/local/mpi/lib\"],
                  casc_mpi_lib_dirs=$withval)]
   )

   if test -z "$casc_mpi_libs"; then
      AC_REQUIRE([CASC_GUESS_ARCH])

      dnl * Set everything to known values
      case $ARCH in

         sun4 | solaris)
            case $F77 in
               *g77)
                   if test -z "$casc_mpi_include_dir"; then
                      casc_mpi_include_dir=/usr/local/mpi/lam/h
                   fi

                   if test -z "$casc_mpi_lib_dirs"; then
                      casc_mpi_lib_dirs="/usr/local/mpi/lam/lib"
                   fi

                   casc_mpi_libs="socket mpi trillium args tstdio t";;

               *)

                  if test -z "$casc_mpi_include_dir"; then
                     MPIINCLUDE="-I/usr/local/mpi/mpich/include \
                                 -I/usr/local/mpi/mpich/lib/solaris/ch_p4"
                  fi

                  if test -z "$casc_mpi_lib_dirs"; then
                     casc_mpi_lib_dirs="/usr/local/mpi/mpich/lib/solaris/ch_p4 \
                                       /usr/lib"
                  fi

               casc_mpi_libs="nsl socket mpi";;
               esac

            if test -z "$MPIINCLUDE"; then
               AC_CHECK_HEADER($casc_mpi_include_dir/mpi.h,
                               MPIINCLUDE="-I$casc_mpi_include_dir")
            fi
         ;;

         alpha)
            if test -z "$casc_mpi_include_dir"; then
               casc_mpi_include_dir=/usr/local/mpi/include
            fi
            AC_CHECK_HEADER($casc_mpi_include_dir/mpi.h,
                               MPIINCLUDE="-I$casc_mpi_include_dir")

            if test -z "$casc_mpi_lib_dirs"; then
               casc_mpi_lib_dirs="/usr/local/mpi/lib/alpha/ch_shmem \
                                  /usr/local/lib"
            fi

            casc_mpi_libs="mpich gs";;

         rs6000)

dnl            if test -z "$casc_mpi_include_dir"; then
dnl               casc_mpi_include_dir=/usr/lpp/ppe.poe/include
dnl            fi
dnl            AC_CHECK_HEADER($casc_mpi_include_dir/mpi.h,
dnl                               MPIINCLUDE="-I$casc_mpi_include_dir")

dnl            if test -z "$casc_mpi_lib_dirs"; then
dnl               casc_mpi_lib_dirs=/usr/lpp/ppe.poe/lib
dnl            fi

            casc_mpi_libs=mpi

            MPIFLAGS="-binitfini:poe_remote_main";;

         IRIX64 | iris4d)
            if test -z "$casc_mpi_include_dir"; then
               casc_mpi_include_dir=/usr/local/mpi/include
            fi
            AC_CHECK_HEADER($casc_mpi_include_dir/mpi.h,
                               MPIINCLUDE="-I$casc_mpi_include_dir")

            if test -z "$casc_mpi_lib_dirs"; then
               casc_mpi_lib_dirs=/usr/local/mpi/lib/IRIX64/ch_p4
            fi

            casc_mpi_libs=mpi;;

         *)
AC_MSG_WARN([trying naive MPI settings - can use --with flags to change])
            if test -z "$casc_mpi_include_dir"; then
               casc_mpi_include_dir=/usr/local/mpi/include
            fi
            AC_CHECK_HEADER($casc_mpi_include_dir/mpi.h,
                               MPIINCLUDE="-I$casc_mpi_include_dir")

            if test -z "$casc_mpi_lib_dirs"; then
               casc_mpi_lib_dirs=/usr/local/mpi/lib
            fi
            casc_mpi_libs=mpi ;;
      esac

      for casc_lib in $casc_mpi_libs; do
         CASC_ADD_LIB($casc_lib, main, $casc_mpi_lib_dirs, MPI)
      done

   else
      if test -n "$casc_mpi_include_dir"; then
         MPIINCLUDE="-I$casc_mpi_include_dir"
      else
         MPIINCLUDE=
      fi

      if test -n "$casc_mpi_lib_dirs"; then
         for casc_lib_dir in $casc_mpi_lib_dirs; do
            MPILIBDIRS="-L$casc_lib_dir $MPILIBDIRS"
         done
      else
         MPILIBDIRS=
      fi

      for casc_lib in $casc_mpi_libs; do
         MPILIBS="$MPILIBS -l$casc_lib"
      done
   fi
])dnl


dnl ********************************************************************
dnl * CASC_FIND_MPI will determine the libraries, directories, and other
dnl * flags needed to compile and link programs with MPI function calls.
dnl * This macro runs tests on the script found by the CASC_PROG_MPICC
dnl * macro.  If there is no such mpicc-type script in the PATH and
dnl * MPICC is not set manually, then this macro will not work.
dnl *
dnl * One may question why these settings would need to be determined if
dnl * there already is mpicc available, and that is a valid question.  I
dnl * can think of a couple of reasons one may want to use these settings
dnl * rather than using mpicc directly.  First, these settings allow you
dnl * to choose the C compiler you wish to use rather than using whatever
dnl * compiler is written into mpicc.  Also, the settings determined by
dnl * this macro should also work with C++ and Fortran compilers, so you
dnl * won't need to have mpiCC and mpif77 alongside mpicc.  This is
dnl * especially helpful on systems that don't have mpiCC.  The advantage
dnl * of this macro over CASC_SET_MPI is that this one doesn't require
dnl * a test of the machine type and thus will hopefully work on unknown
dnl * architectures.  The main disadvantage is that it relies on mpicc.
dnl *
dnl * --with-mpi-include, --with-mpi-libs, and --with-mpi-lib-dirs can be
dnl * used to manually override the automatic test, just as with
dnl * CASC_SET_MPI.  If any one of these three options are used, the
dnl * automatic test will not be run, so it is best to call all three
dnl * whenever one is called.  In addition, the option --with-mpi-flags is
dnl * available here to set any other flags that may be needed, but it
dnl * does not override the automatic test.  Flags set by --with-mpi-flags
dnl * will be added to the variable MPIFLAGS.  This way, if the macro, for
dnl * whatever reason, leaves off a necessary flag, the flag can be added
dnl * to MPIFLAGS without eliminating anything else.  The other variables
dnl * set are MPIINCLUDE, MPILIBS, and MPILIBDIRS, just as in
dnl * CASC_SET_MPI.  This macro also incorporates CASC_SET_MPI as a backup
dnl * plan, where if there is no mpicc, it will use the settings
dnl * determined by architecture name in CASC_SET_MPI
dnl ********************************************************************

AC_DEFUN([CASC_FIND_MPI],
[

   casc_find_mpi_cache_used=yes

   AC_CACHE_VAL(casc_cv_mpi_include, casc_find_mpi_cache_used=no)
   AC_CACHE_VAL(casc_cv_mpi_libs, casc_find_mpi_cache_used=no)
   AC_CACHE_VAL(casc_cv_mpi_lib_dirs, casc_find_mpi_cache_used=no)
   AC_CACHE_VAL(casc_cv_mpi_flags, casc_find_mpi_cache_used=no)

   if test "$casc_find_mpi_cache_used" = "yes"; then
      AC_MSG_CHECKING(for location of mpi.h)
      MPIINCLUDE=$casc_cv_mpi_include
      AC_MSG_RESULT("\(cached\) $MPIINCLUDE")

      AC_MSG_CHECKING(for MPI library directories)
      MPILIBDIRS=$casc_cv_mpi_lib_dirs
      AC_MSG_RESULT("\(cached\) $MPILIBDIRS")

      AC_MSG_CHECKING(for MPI libraries)
      MPILIBS=$casc_cv_mpi_libs
      AC_MSG_RESULT("\(cached\) $MPILIBS")

      AC_MSG_CHECKING(for other MPI-related flags)
      MPIFLAGS=$casc_cv_mpi_flags
      AC_MSG_RESULT("\(cached\) $MPIFLAGS")
   else


      dnl * Set up user options.  If user uses any of the fist three options,
      dnl * then automatic tests are not run.

      casc_user_chose_mpi=no
      AC_ARG_WITH(mpi-include, [  --with-mpi-include=DIR  mpi.h is in DIR],
                  for mpi_dir in $withval; do
                     MPIINCLUDE="$MPIINCLUDE -I$withval"
                  done; casc_user_chose_mpi=yes)

      AC_ARG_WITH(mpi-libs,
[  --with-mpi-libs=LIBS    LIBS is space-separated list of library names
                          needed for MPI, e.g. \"nsl socket mpi\"],
                  for mpi_lib in $withval; do
                     MPILIBS="$MPILIBS -l$mpi_lib"
                  done; casc_user_chose_mpi=yes)


      AC_ARG_WITH(mpi-lib-dirs,
[  --with-mpi-lib-dirs=DIRS
                          DIRS is space-separated list of directories
                          containing the libraries specified by
                          \`--with-mpi-libs', e.g \"/usr/lib /usr/local/mpi/lib\"],
                  for mpi_lib_dir in $withval; do
                     MPILIBDIRS="-L$mpi_lib_dir $MPILIBDIRS"
                  done; casc_user_chose_mpi=yes)

      dnl * --with-mpi-flags only adds to automatic selections,
      dnl * does not override

      AC_ARG_WITH(mpi-flags,
[  --with-mpi-flags=FLAGS  FLAGS is space-separated list of whatever flags other
                          than -l and -L are needed to link with mpi libraries],
                          MPIFLAGS=$withval)


      if test "$casc_user_chose_mpi" = "no"; then

      dnl * Find an MPICC.  If there is none, call CASC_SET_MPI to choose MPI
      dnl * settings based on architecture name.  If CASC_SET_MPI fails,
      dnl * print warning message.  Manual MPI settings must be used.

         AC_ARG_WITH(MPICC,
[  --with-MPICC=ARG        ARG is mpicc or similar MPI C compiling tool],
            MPICC=$withval,
            [AC_CHECK_PROGS(MPICC, mpcc mpicc tmcc hcc)])

         if test -z "$MPICC"; then
            AC_MSG_WARN([no acceptable mpicc found in \$PATH])
            CASC_SET_MPI
            if test -z "$MPILIBS"; then
             AC_MSG_WARN([MPI not found - must set manually using --with flags])
            fi

         dnl * When $MPICC is there, run the automatic test
         dnl * here begins the hairy stuff

         else

            dnl changequote(, )dnl

            dnl * Create a minimal MPI program.  It will be compiled using
            dnl * $MPICC with verbose output.
            cat > mpconftest.c << EOF
#include "mpi.h"

main(int argc, char **argv)
{
   int rank, size;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Finalize();
   return 0;
}
EOF

            casc_mplibs=
            casc_mplibdirs=
            casc_flags=
            casc_lmpi_exists=no

            dnl * These are various ways to produce verbose output from $MPICC
            dnl * All of their outputs are stuffed into variable
            dnl * $casc_mpoutput

            for casc_command in "$MPICC -show"\
                                "$MPICC -v"\
                                "$MPICC -#"\
                                "$MPICC"; do

               casc_this_output=`$casc_command mpconftest.c -o mpconftest 2>&1`

               dnl * If $MPICC uses xlc, then commas must be removed from output
               xlc_p=`echo $casc_this_output | grep xlcentry`
               if test -n "$xlc_p"; then
                  casc_this_output=`echo $casc_this_output | sed 's/,/ /g'`
               fi

               dnl * Turn on flag once -lmpi is found in output
               lmpi_p=`echo $casc_this_output | grep "\-lmpi"`
               if test -n "$lmpi_p"; then
                  casc_lmpi_exists=yes
               fi

               casc_mpoutput="$casc_mpoutput $casc_this_output"
               casc_this_output=

            done

            rm -rf mpconftest*

            dnl * little test to identify $CC as IBM's xlc
            echo "main() {}" > cc_conftest.c
            cc_output=`${CC-cc} -v -o cc_conftest cc_conftest.c 2>&1`
            xlc_p=`echo $cc_output | grep xlcentry`
            if test -n "$xlc_p"; then
               casc_compiler_is_xlc=yes
            fi
            rm -rf cc_conftest*

            dnl * $MPICC might not produce '-lmpi', but we still need it.
            dnl * Add -lmpi to $casc_mplibs if it was never found
            if test "$casc_lmpi_exists" = "no"; then
               casc_mplibs="-lmpi"
            else
               casc_mplibs=
            fi

            casc_want_arg=

            dnl * Loop through every word in output to find possible flags.
            dnl * If the word is the absolute path of a library, it is added
            dnl * to $casc_flags.  Any "-llib", "-L/dir", "-R/dir" and
            dnl * "-I/dir" is kept.  If '-l', '-L', '-R', '-I', '-u', or '-Y'
            dnl * appears alone, then the next word is checked.  If the next
            dnl * word is another flag beginning with '-', then the first
            dnl * word is discarded.  If the next word is anything else, then
            dnl * the two words are coupled in the $casc_arg variable.
            dnl * "-binitfini:poe_remote_main" is a flag needed especially
            dnl * for IBM MPI, and it is always kept if it is found.
            dnl * Any other word is discarded.  Also, after a word is found
            dnl * and kept once, it is discarded if it appears again

            for casc_arg in $casc_mpoutput; do

               casc_old_want_arg=$casc_want_arg
               casc_want_arg=

               if test -n "$casc_old_want_arg"; then
                  case "$casc_arg" in
                  [-*)]
                     casc_old_want_arg=
                  ;;
                  esac
               fi

               case "$casc_old_want_arg" in
               ['')]
                  case $casc_arg in
                  [/*.a)]
                     exists=false
                     for f in $casc_flags; do
                        if test x$casc_arg = x$f; then
                           exists=true
                        fi
                     done
                     if $exists; then
                        casc_arg=
                     else
                        casc_flags="$casc_flags $casc_arg"
                     fi
                  ;;
                  [-binitfini:poe_remote_main)]
                     exists=false
                     for f in $casc_flags; do
                        if test x$casc_arg = x$f; then
                           exists=true
                        fi
                     done
                     if $exists; then
                        casc_arg=
                     else
                        casc_flags="$casc_flags $casc_arg"
                     fi
                  ;;
                  [-lang*)]
                     casc_arg=
                  ;;
                  [-[lLR])]
                     casc_want_arg=$casc_arg
                     casc_arg=
                  ;;
                  [-[lLR]*)]
                     exists=false
                     for f in $casc_flags; do
                        if test x$casc_arg = x$f; then
                           exists=true
                        fi
                     done
                     if $exists; then
                        casc_arg=
                     else
                       casc_flags="$casc_flags $casc_arg"
                     fi
                  ;;
                 [-u)]
                     casc_want_arg=$casc_arg
                     casc_arg=
                  ;;
                  [-Y)]
                     casc_want_arg=$casc_arg
                     casc_arg=
                  ;;
                  [-I)]
                     casc_want_arg=$casc_arg
                     casc_arg=
                  ;;
                  [-I*)]
                     exists=false
                     for f in $casc_flags; do
                        if test x$casc_arg = x$f; then
                           exists=true
                        fi
                     done
                     if $exists; then
                        casc_arg=
                     else
                        casc_flags="$casc_flags $casc_arg"
                     fi
                  ;;
                  [*)]
                     casc_arg=
                  ;;
                  esac

               ;;
               [-[lLRI])]
                  casc_arg="casc_old_want_arg $casc_arg"
               ;;
               [-u)]
                  casc_arg="-u $casc_arg"
               ;;
               [-Y)]
                  casc_arg=`echo $casc_arg | sed -e 's%^P,%%'`
                  SAVE_IFS=$IFS
                  IFS=:
                  casc_list=
                  for casc_elt in $casc_arg; do
                     casc_list="$casc_list -L$casc_elt"
                  done
                  IFS=$SAVE_IFS
                  casc_arg="$casc_list"
               ;;
               esac

               dnl * Still inside the big for loop, we separate each flag
               dnl * into includes, libdirs, libs, flags
               if test -n "$casc_arg"; then
                  case $casc_arg in
                  [-I*)]

                     dnl * if the directory given in this flag contains mpi.h
                     dnl * then the flag is assigned to $MPIINCLUDE
                     if test -z "$MPIINCLUDE"; then
                        casc_cppflags="$casc_cppflags $casc_arg"
                        casc_include_dir=`echo "$casc_arg" | sed 's/-I//g'`

                        SAVE_CPPFLAGS="$CPPFLAGS"
                        CPPFLAGS="$casc_cppflags"
                        dnl changequote([, ])dnl

                        unset ac_cv_header_mpi_h
                        AC_CHECK_HEADER(mpi.h,
                                        MPIINCLUDE="$casc_cppflags")

                        dnl changequote(, )dnl
                        CPPFLAGS="$SAVE_CPPFLAGS"

                     else
                        casc_arg=
                     fi
                  ;;
                  [-[LR]*)]

                     dnl * These are the lib directory flags
                     casc_mplibdirs="$casc_mplibdirs $casc_arg"
                  ;;
                  [-l* | /*)]

                     dnl * These are the libraries
                     casc_mplibs="$casc_mplibs $casc_arg"
                  ;;
                  [-binitfini:poe_remote_main)]
                     if test "$casc_compiler_is_xlc" = "yes"; then
                        casc_mpflags="$casc_mpflags $casc_arg"
                     fi
                  ;;
                  [*)]
                     dnl * any other flag that has been kept goes here
                     casc_mpflags="$casc_mpflags $casc_arg"
                  ;;
                  esac

                  dnl * Upcoming test needs $LIBS to contain the flags
                  dnl * we've found
                  LIBS_SAVE=$LIBS
                  LIBS="$MPIINCLUDE $casc_mpflags $casc_mplibdirs $casc_mplibs"

                  if test -n "`echo $LIBS | grep '\-R/'`"; then
                     LIBS=`echo $LIBS | sed 's/-R\//-R \//'`
                  fi

                  dnl changequote([, ])dnl


                  dnl * Test to see if flags found up to this point are
                  dnl * sufficient to compile and link test program.  If not,
                  dnl * the loop keeps going to the next word
                  AC_LANG_PUSH(C)
                  AC_TRY_LINK(
dnl                     ifelse(AC_LANG, [C++],

dnl [#ifdef __cplusplus
dnl extern "C"
dnl #endif
dnl ])dnl
[#include "mpi.h"
], [int rank, size;
   int argc;
   char **argv;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Finalize();
],
                     casc_result=yes)
                  AC_LANG_POP(C)
                  LIBS=$LIBS_SAVE

                  if test "$casc_result" = yes; then
                     casc_result=
                     break
                  fi
               fi
            done

            dnl * After loop is done, set variables to be substituted
            MPILIBS=$casc_mplibs
            MPILIBDIRS=$casc_mplibdirs
            MPIFLAGS="$MPIFLAGS $casc_mpflags"

            dnl * IBM MPI uses /usr/lpp/ppe.poe/libc.a instead of /lib/libc.a
            dnl * so we need to make sure that -L/lib is not part of the
            dnl * linking line when we use IBM MPI.  This only appears in
            dnl * configure when CASC_FIND_MPI is called first.
	    dnl            ifdef([AC_PROVIDE_CASC_FIND_F77LIBS],
            dnl               if test -n "`echo $F77LIBFLAGS | grep '\-L/lib '`"; then
            dnl                  if test -n "`echo $F77LIBFLAGS | grep xlf`"; then
            dnl                     F77LIBFLAGS=`echo $F77LIBFLAGS | sed 's/-L\/lib //g'`
            dnl                  fi
            dnl               fi
            dnl            )

            if test -n "`echo $MPILIBS | grep pmpich`" &&
               test -z "`echo $MPILIBS | grep pthread`"; then
                  LIBS_SAVE=$LIBS
                  LIBS="$MPIINCLUDE $MPIFLAGS $MPILIBDIRS $MPILIBS -lpthread"
                  AC_LANG_PUSH(C)
                  AC_TRY_LINK(
dnl                     ifelse(AC_LANG, [C++],

dnl [#ifdef __cplusplus
dnl extern "C"
dnl #endif
dnl ])dnl
[#include "mpi.h"
], [int rank, size;
   int argc;
   char **argv;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Finalize();
],
                     MPILIBS="$MPILIBS -lpthread")
                  AC_LANG_POP(C)
                  LIBS=$LIBS_SAVE
            fi

            AC_MSG_CHECKING(for MPI include directories)
            AC_MSG_RESULT($MPIINCLUDE)
            AC_MSG_CHECKING(for MPI library directories)
            AC_MSG_RESULT($MPILIBDIRS)
            AC_MSG_CHECKING(for MPI libraries)
            AC_MSG_RESULT($MPILIBS)
            AC_MSG_CHECKING(for other MPI-related flags)
            AC_MSG_RESULT($MPIFLAGS)
         fi
      fi

      AC_CACHE_VAL(casc_cv_mpi_include, casc_cv_mpi_include=$MPIINCLUDE)
      AC_CACHE_VAL(casc_cv_mpi_lib_dirs, casc_cv_mpi_lib_dirs=$MPILIBDIRS)
      AC_CACHE_VAL(casc_cv_mpi_libs, casc_cv_mpi_libs=$MPILIBS)
      AC_CACHE_VAL(casc_cv_mpi_flags, casc_cv_mpi_flags=$MPIFLAGS)
   fi

   AC_SUBST(MPIINCLUDE)
   AC_SUBST(MPILIBDIRS)
   AC_SUBST(MPILIBS)
   AC_SUBST(MPIFLAGS)

])dnl

dnl ********************************************************************
dnl * CASC_FIND_MPI_ALPHA is a special case of CASC_FIND_MPI for the
dnl * compass cluster.  The original CASC_FIND_MPI looks for existence
dnl * of mpCC and mpiCC.  If the former is found it uses native (proprietary)
dnl * mpi and if the latter is found, it uses mpich.  The DECs are a
dnl * special case because mpCC does not exist and mpiCC does, but we want
dnl * to use the native version by default.  Therefore, the original macro
dnl * did not work for this case so I added this one to deal with it.
dnl * AMW 9/00
dnl ********************************************************************

AC_DEFUN([CASC_FIND_MPI_ALPHA],
[

   casc_find_mpi_cache_used=yes

   AC_CACHE_VAL(casc_cv_mpi_include, casc_find_mpi_cache_used=no)
   AC_CACHE_VAL(casc_cv_mpi_libs, casc_find_mpi_cache_used=no)
   AC_CACHE_VAL(casc_cv_mpi_lib_dirs, casc_find_mpi_cache_used=no)
   AC_CACHE_VAL(casc_cv_mpi_flags, casc_find_mpi_cache_used=no)

   if test "$casc_find_mpi_cache_used" = "yes"; then
      AC_MSG_CHECKING(for location of mpi.h)
      MPIINCLUDE=$casc_cv_mpi_include
      AC_MSG_RESULT("\(cached\) $MPIINCLUDE")

      AC_MSG_CHECKING(for MPI library directories)
      MPILIBDIRS=$casc_cv_mpi_lib_dirs
      AC_MSG_RESULT("\(cached\) $MPILIBDIRS")

      AC_MSG_CHECKING(for MPI libraries)
      MPILIBS=$casc_cv_mpi_libs
      AC_MSG_RESULT("\(cached\) $MPILIBS")

      AC_MSG_CHECKING(for other MPI-related flags)
      MPIFLAGS=$casc_cv_mpi_flags
      AC_MSG_RESULT("\(cached\) $MPIFLAGS")
   else


      dnl * Set up user options.  If user uses any of the fist three options,
      dnl * then automatic tests are not run.

      casc_user_chose_mpi=no
      AC_ARG_WITH(mpi-include, [  --with-mpi-include=DIR  mpi.h is in DIR],
                  for mpi_dir in $withval; do
                     MPIINCLUDE="$MPIINCLUDE -I$withval"
                  done; casc_user_chose_mpi=yes)

      AC_ARG_WITH(mpi-libs,
[  --with-mpi-libs=LIBS    LIBS is space-separated list of library names
                          needed for MPI, e.g. \"nsl socket mpi\"],
                  for mpi_lib in $withval; do
                     MPILIBS="$MPILIBS -l$mpi_lib"
                  done; casc_user_chose_mpi=yes)


      AC_ARG_WITH(mpi-lib-dirs,
[  --with-mpi-lib-dirs=DIRS
                          DIRS is space-separated list of directories
                          containing the libraries specified by
                          \`--with-mpi-libs', e.g \"/usr/lib /usr/local/mpi/lib\"],
                  for mpi_lib_dir in $withval; do
                     MPILIBDIRS="-L$mpi_lib_dir $MPILIBDIRS"
                  done; casc_user_chose_mpi=yes)

      dnl * --with-mpi-flags only adds to automatic selections,
      dnl * does not override

      AC_ARG_WITH(mpi-flags,
[  --with-mpi-flags=FLAGS  FLAGS is space-separated list of whatever flags other
                          than -l and -L are needed to link with mpi libraries],
                          MPIFLAGS=$withval)


      if test "$casc_user_chose_mpi" = "no"; then

         dnl * Set defaults for Compass cluster here.  This is the point where
         dnl * we call CASC_SET_MPI in CASC_FIND_MPI macro.

         casc_mpi_include_dir=
         casc_mpi_lib_dirs=
         casc_mpi_libs="mpi rt rpc gs pthread"

         for casc_incl_dir in $casc_mpi_include_dir; do
            MPIINCLUDE="-I$casc_incl_dir $MPIINCLUDE"
         done
         for casc_lib_dir in $casc_mpi_lib_dirs; do
            MPILIBDIRS="-L$casc_lib_dir $MPILIBDIRS"
         done
         for casc_lib in $casc_mpi_libs; do
            MPILIBS="$MPILIBS -l$casc_lib"
         done
      fi


      AC_MSG_CHECKING(for MPI include directories)
      AC_MSG_RESULT($MPIINCLUDE)
      AC_MSG_CHECKING(for MPI library directories)
      AC_MSG_RESULT($MPILIBDIRS)
      AC_MSG_CHECKING(for MPI libraries)
      AC_MSG_RESULT($MPILIBS)
      AC_MSG_CHECKING(for other MPI-related flags)
      AC_MSG_RESULT($MPIFLAGS)

   fi

])dnl



# ===========================================================================
#
# SYNOPSIS
#
#   CHECK_SZLIB()
#
# DESCRIPTION
#
#   This macro searches for an installed szlib library. If nothing was
#   specified when calling configure, it searches first in /usr/local and
#   then in /usr. If the --with-szlib=DIR is specified, it will try to find
#   it in DIR/include/szlib.h and DIR/lib/libsz.a. If --without-szlib is
#   specified, the library is not searched at all.
#
#   If either the header file (szlib.h) or the library (libsz) is not found,
#   the configuration exits on error, asking for a valid szlib installation
#   directory or --without-szlib.
#
#   The macro defines the symbol HAVE_LIBSZ if the library is found. You
#   should use autoheader to include a definition for this symbol in a
#   config.h file. Sample usage in a C/C++ source is as follows:
#
#     #ifdef HAVE_LIBSZ
#     #include <szlib.h>
#     #endif /* HAVE_LIBSZ */
#

AC_DEFUN([CHECK_SZLIB],
#
# DEFINES :
#	        szlib_PREFIX
#		szlib_INCLUDES
#		szlib_LIBS
#
[AC_MSG_CHECKING(if szlib is wanted)
AC_ARG_WITH(szlib,
[  --with-szlib=DIR root directory path of szlib installation [DIR defaults to
                    /usr/local or /usr if not found in /usr/local]
  --without-szlib to disable szlib usage completely [the default]],
[if test "$withval" != no ; then
  AC_MSG_RESULT(yes)
  if test "$withval" == yes ;
  then
     SZLIB_HOME=/usr/local
  else
     SZLIB_HOME="$withval"
  fi
  if test ! -d "$SZLIB_HOME"
  then
    AC_MSG_WARN([Sorry, $SZLIB_HOME does not exist, checking usual places])
    SZLIB_HOME=/usr/local
    if test ! -f "${SZLIB_HOME}/include/szlib.h"
    then
       SZLIB_HOME=/usr
    fi
  fi
else
  AC_MSG_RESULT(no)
fi])

#
# Locate szlib, if wanted
#
if test -n "${SZLIB_HOME}"
then
        SZLIB_OLD_LDFLAGS=$LDFLAGS
        SZLIB_OLD_CPPFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS -L${SZLIB_HOME}/lib"
        CPPFLAGS="$CPPFLAGS -I${SZLIB_HOME}/include"
        AC_LANG_SAVE
        AC_LANG_C
        AC_CHECK_LIB(z, inflateEnd, [szlib_cv_libsz=yes], [szlib_cv_libsz=no])
        AC_CHECK_HEADER(szlib.h, [szlib_cv_szlib_h=yes], [szlib_cv_szlib_h=no])
        AC_LANG_RESTORE
        if test "$szlib_cv_libsz" = "yes" -a "$szlib_cv_szlib_h" = "yes"
        then
	        szlib_PREFIX="${SZLIB_HOME}"
		szlib_INCLUDES="-I${SZLIB_HOME}/include"
		szlib_LIBS="-L${SZLIB_HOME}/lib -lsz"
                #
                # If both library and header were found, use them
                #
                AC_CHECK_LIB(z, inflateEnd)
                AC_MSG_CHECKING(szlib in ${SZLIB_HOME})
                AC_MSG_RESULT(ok)
        else
                #
                # If either header or library was not found, revert and bomb
                #
                AC_MSG_CHECKING(szlib in ${SZLIB_HOME})
                LDFLAGS="$SZLIB_OLD_LDFLAGS"
                CPPFLAGS="$SZLIB_OLD_CPPFLAGS"
                AC_MSG_RESULT(failed)
                AC_MSG_ERROR(either specify a valid szlib installation with --with-szlib=DIR or disable szlib usage with --without-szlib)
        fi
fi

])

# ===========================================================================
#           http://www.nongnu.org/autoconf-archive/check_zlib.html
# ===========================================================================
#
# SYNOPSIS
#
#   CHECK_ZLIB()
#
# DESCRIPTION
#
#   This macro searches for an installed zlib library. If nothing was
#   specified when calling configure, it searches first in /usr/local and
#   then in /usr. If the --with-zlib=DIR is specified, it will try to find
#   it in DIR/include/zlib.h and DIR/lib/libz.a. If --without-zlib is
#   specified, the library is not searched at all.
#
#   If either the header file (zlib.h) or the library (libz) is not found,
#   the configuration exits on error, asking for a valid zlib installation
#   directory or --without-zlib.
#
#   The macro defines the symbol HAVE_LIBZ if the library is found. You
#   should use autoheader to include a definition for this symbol in a
#   config.h file. Sample usage in a C/C++ source is as follows:
#
#     #ifdef HAVE_LIBZ
#     #include <zlib.h>
#     #endif /* HAVE_LIBZ */
#
# LICENSE
#
#   Copyright (c) 2008 Loic Dachary <loic@senga.org>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([CHECK_ZLIB],
#
# DEFINES :
#	        zlib_PREFIX
#		zlib_INCLUDES
#		zlib_LIBS
#
[AC_MSG_CHECKING(if zlib is wanted)
AC_ARG_WITH(zlib,
[  --with-zlib=DIR root directory path of zlib installation [DIR defaults to
                    /usr/local or /usr if not found in /usr/local]
  --without-zlib to disable zlib usage completely [the default]],
[if test "$withval" != no ; then
  AC_MSG_RESULT(yes)
  if test "$withval" == yes ;
  then
     ZLIB_HOME=/usr/local
  else
     ZLIB_HOME="$withval"
  fi
  if test ! -d "$ZLIB_HOME"
  then
    AC_MSG_WARN([Sorry, $ZLIB_HOME does not exist, checking usual places])
    ZLIB_HOME=/usr/local
    if test ! -f "${ZLIB_HOME}/include/zlib.h"
    then
       ZLIB_HOME=/usr
    fi
  fi
else
  AC_MSG_RESULT(no)
fi],
  [AC_MSG_RESULT(no)]
)


#
# Locate zlib, if wanted
#
if test -n "${ZLIB_HOME}"
then
        ZLIB_OLD_LDFLAGS=$LDFLAGS
        ZLIB_OLD_CPPFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS -L${ZLIB_HOME}/lib"
        CPPFLAGS="$CPPFLAGS -I${ZLIB_HOME}/include"
        AC_LANG_SAVE
        AC_LANG_C
        AC_CHECK_LIB(z, inflateEnd, [zlib_cv_libz=yes], [zlib_cv_libz=no])
        AC_CHECK_HEADER(zlib.h, [zlib_cv_zlib_h=yes], [zlib_cv_zlib_h=no])
        AC_LANG_RESTORE
        if test "$zlib_cv_libz" = "yes" -a "$zlib_cv_zlib_h" = "yes"
        then
	        zlib_PREFIX="${ZLIB_HOME}"
		zlib_INCLUDES="-I${ZLIB_HOME}/include"
		zlib_LIBS="-L${ZLIB_HOME}/lib -lz"
                #
                # If both library and header were found, use them
                #
                AC_CHECK_LIB(z, inflateEnd)
                AC_MSG_CHECKING(zlib in ${ZLIB_HOME})
                AC_MSG_RESULT(ok)
        else
                #
                # If either header or library was not found, revert and bomb
                #
                AC_MSG_CHECKING(zlib in ${ZLIB_HOME})
                LDFLAGS="$ZLIB_OLD_LDFLAGS"
                CPPFLAGS="$ZLIB_OLD_CPPFLAGS"
                AC_MSG_RESULT(failed)
                AC_MSG_ERROR(either specify a valid zlib installation with --with-zlib=DIR or disable zlib usage with --without-zlib)
        fi
fi

])

dnl $Id$

dnl Determines which compiler is being used.
dnl This check uses the compiler behavior when possible.
dnl For some compiler, we resort to a best guess,
dnl because we do not know a foolproof way to get the info.

dnl Much of the information used here came from the very
dnl helpful predef project (http://predef.sourceforge.net/).




dnl Simple wrappers to allow using CASC_INFO_CXX_ID_NAMES and
dnl CASC_INFO_CC_ID_NAMES without arguments.
dnl The names CC_ID and CC_VERSION are used for the C compiler id and version.
dnl The names CXX_ID and CXX_VERSION are used for the C++ compiler id and version.
AC_DEFUN([CASC_INFO_CXX_ID],[
  CASC_INFO_CXX_ID_NAMES(CXX_ID,CXX_VERSION)
])
AC_DEFUN([CASC_INFO_CC_ID],[
  CASC_INFO_CC_ID_NAMES(CC_ID,CC_VERSION)
])
AC_DEFUN([CASC_INFO_CC_CXX_ID],[
  AC_REQUIRE([CASC_INFO_CC_ID])
  AC_REQUIRE([CASC_INFO_CXX_ID])
])


dnl CASC_INFO_CXX_ID and CASC_INFO_C_ID determine which C or C++ compiler
dnl is being used.
# Set the variables CXX_ID or C_ID as follows:
# Gnu		-> gnu
# SUNWspro	-> sunpro
# Dec		-> dec
# KCC		-> kai
# Intel		-> intel
# SGI		-> sgi
# IBM xlc	-> xlc


AC_DEFUN([CASC_INFO_CXX_ID_NAMES],
dnl Arguments are:
dnl 1. Name of variable to set to the ID string.
dnl 2. Name of variable to set to the version number.
[
# Start macro CASC_INFO_CXX_ID_NAMES
  AC_REQUIRE([AC_PROG_CXXCPP])
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  CASC_AC_LOG(CXXP is $CXX)
  CASC_AC_LOG(CXXCPP is $CXXCPP)

  $1=unknown
  $2=unknown

dnl Do not change the following chain of if blocks into a case statement.
dnl We may eventually have a compiler that must be tested in a different
dnl method


  # Check if it is a Sun compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CXX is sunpro)
    AC_EGREP_CPP([^0x[0-9]+],__SUNPRO_CC,
      $1=sunpro
      # SUN compiler defines __SUNPRO_CC to the version number.
      echo __SUNPRO_CC > conftest.C
      $2=`${CXXCPP} conftest.C | sed -n 2p`
      rm -f conftest.C
    )
  fi


  # Check if it is a Intel compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CXX is intel)
    AC_EGREP_CPP(^yes,
#ifdef __INTEL_COMPILER
yes;
#endif
,
      $1=intel
      # Intel compiler defines __INTEL_COMPILER to the version number.
      echo __INTEL_COMPILER > conftest.C
      $2=`${CXXCPP} conftest.C | sed -n 2p`
      rm -f conftest.C
    )
  fi


  # Check if it is a GNU compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CXX is gnu)
    AC_EGREP_CPP(^yes,
#ifdef __GNUC__
yes;
#endif
,
    $1=gnu
    # GNU compilers output version number with option --version.
    # Alternatively, it also defines the macros __GNUC__,
    # GNUC_MINOR__ and __GNUC_PATCHLEVEL__
    [[$2=`$CXX --version | sed -e 's/[^0-9]\{0,\}\([^ ]\{1,\}\).\{0,\}/\1/' -e 1q`]]
    )
  fi


  # Check if it is a DEC compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CXX is dec)
    AC_EGREP_CPP(^1,__DECCXX,
      $1=dec
      # DEC compiler defines __DECCXX_VER to the version number.
      echo __DECCXX_VER > conftest.C
      $2=`${CXXCPP} conftest.C | sed -n 2p`
      rm -f conftest.C
    )
  fi


  # Check if it is a KAI compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CXX is kai)
    AC_EGREP_CPP(^1,__KCC,
      $1=kai
      # KCC compiler defines __KCC_VERSION to the version number.
      echo __KCC_VERSION > conftest.C
      $2=`${CXXCPP} conftest.C | sed -n 2p`
      rm -f conftest.C
    )
  fi


  # Check if it is a SGI compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CXX is sgi)
    AC_EGREP_CPP(^1,__sgi,
      $1=sgi
      # SGI compiler defines _COMPILER_VERSION to the version number.
      echo _COMPILER_VERSION > conftest.C
      $2=`${CXXCPP} conftest.C | sed /^\\#/d`
      rm -f conftest.C
    )
  fi


  # Check if it is a IBM compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CXX is xlc)
    AC_EGREP_CPP(^yes,
#ifdef __xlC__
yes;
#endif
,
    $1=xlc
    # IBM compiler defines __xlC__ to the version number.
    echo __xlC__ > conftest.C
    $2=`${CXXCPP} conftest.C | sed /^\\#/d`
    rm -f conftest.C
    )
  fi


  AC_LANG_RESTORE
  CASC_AC_LOG_VAR(CXX_ID CXX_VERSION)
# End macro CASC_INFO_CXX_ID_NAMES
])





AC_DEFUN([CASC_INFO_CC_ID_NAMES],
dnl Arguments are:
dnl 1. Name of variable to set to the ID string.
dnl 2. Name of variable to set to the version number.
[
# Start macro CASC_INFO_CC_ID_NAMES
  AC_REQUIRE([AC_PROG_CPP])
  AC_LANG_SAVE
  AC_LANG_C
  CASC_AC_LOG(CC is $CC)
  CASC_AC_LOG(CPP is $CPP)

  $1=unknown
  $2=unknown

dnl Do not change the following chain of if blocks into a case statement.
dnl We may eventually have a compiler that must be tested in a different
dnl method


  # Check if it is a Sun compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CC is sunpro)
    AC_EGREP_CPP([^ 0x[0-9]+],__SUNPRO_C,
      $1=sunpro
      # SUN compiler defines __SUNPRO_C to the version number.
      echo __SUNPRO_C > conftest.c
      $2=`${CPP} ${CPPFLAGS} conftest.c | sed -n -e 's/^ //' -e 2p`
      rm -f conftest.c
    )
  fi


  # Check if it is a Intel compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CC is intel)
    AC_EGREP_CPP(^yes,
#ifdef __INTEL_COMPILER
yes;
#endif
,
      $1=intel
      # Intel compiler defines __INTEL_COMPILER to the version number.
      echo __INTEL_COMPILER > conftest.C
      $2=`${CPP} conftest.C | sed -n 2p`
      rm -f conftest.C
    )
  fi


  # Check if it is a GNU compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CC is gnu)
    AC_EGREP_CPP(^yes,
#ifdef __GNUC__
yes;
#endif
,
    $1=gnu
    [[$2=`$CC --version | sed -e 's/[^0-9]\{0,\}\([^ ]\{1,\}\).\{0,\}/\1/' -e 1q`]]
    )
  fi


  # Check if it is a DEC compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CC is dec)
    AC_EGREP_CPP(^ 1,__DECC,
      $1=dec
      # DEC compiler defines __DECC_VER to the version number.
      echo __DECC_VER > conftest.c
      $2=`${CPP} ${CPPFLAGS} conftest.c | sed -n -e 's/^ //' -e 2p`
      rm -f conftest.c
    )
  fi


  # Check if it is a KAI compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CC is kai)
    AC_EGREP_CPP(^1,__KCC,
      $1=kai
      # KCC compiler defines __KCC_VERSION to the version number.
      echo __KCC_VERSION > conftest.c
      $2=`${CPP} ${CPPFLAGS} conftest.c | sed -n 2p`
      rm -f conftest.c
    )
  fi


  # Check if it is a SGI compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CC is sgi)
    AC_EGREP_CPP(^1,__sgi,
      $1=sgi
      # SGI compiler defines _COMPILER_VERSION to the version number.
      echo _COMPILER_VERSION > conftest.c
      $2=`${CPP} ${CPPFLAGS} conftest.c | sed /^\\#/d`
      rm -f conftest.c
    )
  fi


  # Check if it is a IBM compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CC is xlc)
    if echo "$host_os" | grep "aix" >/dev/null ; then
      # The wretched IBM shell does not eval correctly,
      # so we have to help it with a pre-eval eval statement.
      ac_cpp=`eval "echo $ac_cpp"`
      save_ac_cpp=$ac_cpp
      CASC_AC_LOG(ac_cpp is temporarily set to $ac_cpp)
    else
      save_ac_cpp=
    fi
    CASC_AC_LOG(ac_cpp is $ac_cpp)
    AC_EGREP_CPP(^yes,
#ifdef __xlC__
yes;
#endif
,
    $1=xlc
    # IBM compiler defines __xlC__ to the version number.
    echo __xlC__ > conftest.C
    $2=`${CPP} conftest.C | sed /^\\#/d`
    rm -f conftest.C
    )
    test "$save_ac_cpp" && ac_cpp=$save_ac_cpp
    CASC_AC_LOG(ac_cpp is restored to $ac_cpp)
  fi


  AC_LANG_RESTORE
  CASC_AC_LOG_VAR(CC_ID CC_VERSION)
# End macro CASC_INFO_CC_ID_NAMES
])

AC_DEFUN([CASC_AC_LOG],[echo "configure:__oline__:" $1 >&AC_FD_CC])

AC_DEFUN([CASC_AC_LOG_VAR],[
dnl arg1 is list of variables to log.
dnl arg2 (optional) is a label.
dnl
dnl This macro makes code that write out at configure time
dnl label: x is '...'
dnl if x is set and
dnl label: x is unset
dnl otherwise.
define([log_label],ifelse($2,,,[$2: ]))
log_vars="$1"
for log_vars_index in $log_vars ; do
  eval "test \"\${${log_vars_index}+set}\" = set"
  if test $? = 0; then
    log_vars_value="'`eval echo \\${$log_vars_index}`'";
  else
    log_vars_value="unset";
  fi
  CASC_AC_LOG("log_label$log_vars_index is $log_vars_value");
dnl
dnl This is a shorter version, but it does not work for some Bourne shells
dnl due to misinterpretation of the multiple backslashes
dnl CASC_AC_LOG("log_label$log_vars_index is `eval if test \\\"\$\{$log_vars_index+set\}\\\"\; then echo \\\""'"\$\{$og_vars_index\}"'"\\\"\; else echo 'unset'\; fi`")
done
undefine([log_label])
])

AC_DEFUN([SPLIT_LIBS_STRING],[
dnl
dnl This macro takes an automake-style LIBS string (arg1) and
dnl splits it into the -L part (arg2, what is usually called
dnl LIB_PATH) and -l part (arg3, what is usually called LIB_NAME).
dnl
# Split $1 into the LIB_PATH part ($2) and the LIB_NAME part ($3)
if test -n "${$1}"; then
  for i in ${$1}; do
    case "$i" in
    -L*) $2="${$2} $i" ;;
    *) $3="${$3} $i" ;;
    esac
  done
fi
])

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

dnl $Id$


AC_DEFUN([CASC_LIBS_ADD_RPATH],[
# Begin macro CASC_LIBS_ADD_RPATH
dnl Support RPATH by going in a LIBS string and, for each -L flag,
dnl add a flag immediately following it to set the RPATH, for
dnl paths that contain shared libraries.
dnl
dnl arg1 is a LIBS string.
dnl arg2 is the name of the variable to set to the new LIBS string.
dnl arg3 is non-empty to use id of the C++ compiler instead of the C compiler.


dnl Determine which compiler is being used, because
dnl the syntax of the RPATH flag depends on the compiler.
dnl Use the C++ compiler and assume the C compiler
dnl is from the same family.
AC_REQUIRE([CASC_INFO_CC_CXX_ID])


AC_ARG_ENABLE(rpath,
[  --enable-rpath=SYNTAX	When linking add syntax for rpath for every
			-L option that points to a directory with .so
			files in it.  If SYNTAX is omitted, an attempt
			is made to find out the correct rpath syntax for
			the compiler being used.]
,,enable_rpath=yes)

if test "$enable_rpath" = yes; then
  # Determine the proper rpath syntax.

  AC_LANG_SAVE

  ifelse([$3],,
  AC_LANG_C
  rpath_compiler_id="$CC_ID",
  AC_LANG_CPLUSPLUS
  rpath_compiler_id="$CXX_ID"
  )


  # Unset the rpath syntax variable so we can check on whether we
  # found a way to set it.
  unset rpath_beginning;

  # Determine, based on the compiler, the syntax for specifying RPATH.
  # It should be of the form "$rpath_beginning$the_path", where
  # rpath_beginning is the compiler-dependent part.
  case "$rpath_compiler_id" in
    gnu)
      # This compiler may use a variable rpath syntax because it may use
      # the native loader.
      CASC_LIBS_FIND_RPATH(rpath_beginning,
	['---bogus-flag-meant-to-cause-error' '-Wl,-rpath ' '-Wl,-R' '-Wl,-R '])
    ;;
    intel)
      # This compiler may use a variable rpath syntax because it may use
      # the native loader.
      CASC_LIBS_FIND_RPATH(rpath_beginning,
	['---bogus-flag-meant-to-cause-error' '-Wl,-rpath ' '-Wl,-R' '-Wl,-R '])
      if test "$rpath_beginning" = "---bogus-flag-meant-to-cause-error"; then
        # Do not rely on the compiler return value to test for syntax
        # Guess the syntax assuming the native loader will be used.
        case "$host_os" in
          linux*) rpath_beginning='-Wl,-rpath ' ;;
          sun*|solaris*) rpath_beginning='-R' ;;
          osf*) rpath_beginning='-rpath ' ;;
          *) rpath_beginning='' ;;
        esac
        AC_MSG_WARN(
  [Your compiler ifelse($3,,$CC,$CXX) returns 0 even when it is
  given a bogus flag.  Therefore, I cannot find the proper syntax
  for the rpath for this compiler.  I have resorted to a guess that
  may not be correct: '$rpath_beginning'.
  You can override this by using --enable-rpath=SYNTAX])
      fi
    ;;
    sunpro)
      # This compiler may use a variable rpath syntax.
      CASC_LIBS_FIND_RPATH(rpath_beginning,['---bogus-flag-meant-to-cause-error' '-R' '-R '])
    ;;
    kai)
      # The KAI compilers use the system native loader.
      #
      # On some platforms (PC/Linux at least), this compiler seems
      # to return 0 even if it encounters error, thus it can return
      # the first guess for the rpath syntax, even if the guess is
      # wrong.  We try to catch this by making the first flag bogus.
      # If the compiler accepts this flag (by returning 0), we know
      # it is wrong and we resort to an alternative method for
      # getting the rpath syntax.
      CASC_LIBS_FIND_RPATH(rpath_beginning,
	['---bogus-flag-meant-to-cause-error' '-R' '-R ' '-rpath ' '-Wl,-rpath ' '-Wl,-R' '-Wl,-R '])
      if test "$rpath_beginning" = "---bogus-flag-meant-to-cause-error"; then
        # Do not rely on the compiler return value to test for syntax
        # Guess the syntax assuming the native loader will be used.
        case "$host_os" in
          linux*) rpath_beginning='-Wl,-rpath ' ;;
          sun*|solaris*) rpath_beginning='-R' ;;
          osf*) rpath_beginning='-rpath ' ;;
          *) rpath_beginning='' ;;
        esac
        AC_MSG_WARN(
  [Your compiler ifelse($3,,$CC,$CXX) returns 0 even when it is
  given a bogus flag.  Therefore, I cannot find the proper syntax
  for the rpath for this compiler.  I have resorted to a guess that
  may not be correct: '$rpath_beginning'.
  You can override this by using --enable-rpath=SYNTAX])
      fi
    ;;
    *)
      CASC_LIBS_FIND_RPATH(rpath_beginning)
    ;;
  esac
  CASC_AC_LOG_VAR(host_os CC_ID CXX_ID rpath_compiler_id rpath_beginning, forming rpaths)

  AC_LANG_RESTORE

  # It is valid to have rpath_beginning be blank.
  # but if it is unset, we could not find a way to set it.
  if test ! "${rpath_beginning+set}" = set; then
    AC_MSG_WARN(I cannot find a working syntax for setting relocatable paths)
  fi

elif test ! "${enable_rpath}" = no; then

  # User has provided the rpath syntax.
  rpath_beginning=$enable_rpath

fi;	# End block determining the proper rpath syntax.


# Use the rpath syntax.
if test "${rpath_beginning+set}" = set	\
  && test -n "${rpath_beginning}" ; then
  # Add the RPATH flags only if we know the syntax for it,
  # and if it is needed as indicated by a non-empty rpath_beginning.

  # Loop through the flags in $1, looking for the -L flag,
  # and append RPATH flag to each one found, if the the
  # path specified by the flag includes shared libraries.
  for i in ${$1}; do
    new_$2="${new_$2} ${i}"
    tmp_addl_string=`echo $i | sed 's/^-L//'`
    test "$tmp_addl_string" = "$i" && continue	# does not contain -L.
    test -d "$tmp_addl_string" || continue;	# directory nonexistent.
    test "`echo $tmp_addl_string/*.so`" = "$tmp_addl_string/*.so" \
      && continue;	# does not contain shared libraries.
    echo "${new_$2}"	\
      | grep ".*${rpath_beginning}[[ 	]]*${tmp_addl_string}"	\
      > /dev/null	\
      && continue	# already contains the flag we want to add.
    new_$2="${new_$2} ${rpath_beginning}${tmp_addl_string}"
  done
  $2="${new_$2}"

fi

dnl Now, arg2 should be similar to arg1, but with the additional RPATH flags.

# End macro CASC_LIBS_ADD_RPATH
])

AC_DEFUN([CASC_LIBS_FIND_RPATH],[
# Begin macro CASC_LIBS_FIND_RPATH
dnl Find the correct rpath syntax from the list given in arg1.
dnl arg1: variable to set to the syntax string
dnl arg2: list of syntaxes to try;
dnl   if blank, a large number of syntaxes will be tried.
dnl
dnl arg1 is list of possible rpath syntaxes to try.
define(possible_rpaths,dnl
[ifelse($2,,['-R ' '-R' '-rpath ' '-Wl,-rpath ' '-Wl,-R ' '-Wl,-R'],[[$2]])])
  save_LIBS="$LIBS";
  for i in possible_rpaths; do
    LIBS="${i}/usr/local"
    AC_TRY_LINK(,,$1="$i", unset $1)
    # Intel compiler does not fail on bad args but warning message is
    # created. If warning is found in the log then continue searching
    # for syntax as the current one is no good.  If warning is not
    # found use return status of link attempt to determine if
    # parameter was accepted by the compiler.
    SEARCH=`echo "ignoring unknown option '${LIBS}'" | sed -e "s/---/-f-/"`
    if ( grep "$SEARCH" config.log ) >/dev/null 2>&1
    then 
	:
    else 
        if test "${$1+set}" = set; then break; fi
    fi
  done
  LIBS="$save_LIBS"
undefine([possible_rpaths])
# End macro CASC_LIBS_FIND_RPATH
])

