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



