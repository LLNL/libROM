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
