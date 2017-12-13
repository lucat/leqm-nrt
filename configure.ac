#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([leqm-nrt], [0.17], [luca.trisciani@dcp-werkstatt.ch])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_FILES([Makefile src/Makefile])
# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h unistd.h stdio.h math.h sndfile.h unistd.h pthread.h string.h time.h ctype.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([clock_gettime memset])

AC_OUTPUT
