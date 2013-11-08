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
