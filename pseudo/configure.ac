

AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([Pseudo - a gaussian based pseudo-generator],[1.7.5],[Thierry.Deutsch@cea.fr],[pseudo])

dnl Init something, don't know exactly what...
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE(-Wno-portability)

dnl If FC environment variable is not already set, it
dnl looks for a modern Fortran compiler, prefering newest
dnl fortran available (e. g. Fortran95). Then FC variable
dnl is set.
AC_PROG_FC([], [Fortran 90])
AC_FC_SRCEXT(f90)

AX_MPI()
AX_FC_MOD()

AX_LINALG()

AX_LIBXC()

dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
src/Makefile
tests/Makefile
tests/atoms/Makefile
tests/C-pseudo/Makefile
tests/H-pseudo/Makefile
tests/pawpatchtest/Makefile
tests/tutorials/Makefile
tests/tutorials/oxygen.NLCC/Makefile
])

AC_OUTPUT
