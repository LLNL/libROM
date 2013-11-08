dnl Define a macro for supporting HDF5

AC_DEFUN([CASC_SUPPORT_HDF5],[

# Begin CASC_SUPPORT_HDF5
# Defines hdf5_PREFIX hdf5_INCLUDES and hdf5_LIBS if with-hdf5 is specified.
AC_ARG_WITH(hdf5,
[ --with-hdf5[=PATH]  Use HDF5 and optionally specify where HDF5 is installed.],
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

