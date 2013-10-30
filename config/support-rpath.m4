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
