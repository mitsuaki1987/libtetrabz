#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([libtetrabz], [2.0.0], [mkawamura@issp.u-tokyo.ac.jp])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_HEADERS([config.h])

echo ""
echo "  Check compiler"
echo ""
AC_PROG_CC
AC_PROG_F77
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
                                 [compile with MPI (parallelization) support.
                                  If none is found, MPI is not used. Default: no])],,
            [with_mpi=no])
if test x"$with_mpi" != xno; then
  AX_PROG_FC_MPI(,[FCFLAGS="$FCFLAGS -D__MPI"],[
    AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
    ])
else
  AC_PROG_FC
fi
echo ""
echo "  Check OpenMP"
echo ""
AC_ARG_WITH(openmp, [AS_HELP_STRING([--with-openmp],
                                 [compile with OpenMP (parallelization) support.
                                  If none is found, OpenMP is not used. Default: yes])],,
            [with_openmp=yes])
if test x"$with_openmp" != xno; then
  AX_OPENMP([FCFLAGS="$FCFLAGS $OPENMP_CFLAGS";CFLAGS="$CFLAGS $OPENMP_CFLAGS"],[
      AC_MSG_WARN([OpenMP requested, but couldn't use OpenMP.])
  ])
fi
echo ""
echo "  Check fortran module"
echo ""
# Find out how to add to the Fortran include path.
AX_F90_MODULE_FLAG
# The module flag may be "-I " or so, and some make implementations
# swallow trailing space in variable settings, so we employ a little
# trick.
AC_SUBST([empty], [""])
AC_SUBST([FC_MODINC], [$ax_cv_f90_modflag'${empty}'])
# Find out the file name extensions of Fortran module files.
AX_F90_MODULE_EXTENSION
if test "x$ax_f90_modext" != xunknown; then
  FC_MODEXT=$ax_cv_f90_modext
else
  FC_MODEXT=$ax_cv_f90_modext
fi
AC_SUBST([FC_MODEXT])
echo ""
echo "  Check library generation"
echo ""
AM_PROG_AR
AC_PROG_LIBTOOL
#
echo ""
echo "  Output files"
echo ""
AC_CONFIG_FILES([
Makefile
src/Makefile
test/Makefile
example/Makefile
])

AC_OUTPUT

echo ""
echo "  SUMMARY"
echo ""
echo "             FC = ${FC}"
echo "             CC = ${CC}"
echo "        FCFLAGS = ${FCFLAGS}"
echo "         CFLAGS = ${CFLAGS}"
echo ""
echo "  Configure done."
echo ""
