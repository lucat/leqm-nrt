#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.69])
AC_INIT([leqm-nrt], [0.20], [luca.trisciani@dcp-werkstatt.ch])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR([src/leqm-nrt.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC


CFLAGS="$CFLAGS -I/usr/include/di"
# CXXFLAGS="$CXXFLAGS -I/opt/local/include"
LDFLAGS="$LDFLAGS -L/usr/lib/di"

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lavcodec':
AC_CHECK_LIB([avcodec], [main])
# FIXME: Replace `main' with a function in `-lavformat':
AC_CHECK_LIB([avformat], [main])
# FIXME: Replace `main' with a function in `-lavutil':
AC_CHECK_LIB([avutil], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [main])

AS_IF([test "x$test" = xfalse],
[AC_CHECK_LIB([sndfile], [main])],
[AC_MSG_RESULT([Linking with ffmpeg])])


# FIXME: Replace `main' with a function in `-lrt':

AC_CHECK_LIB([rt], [main])
# FIXME: Replace `main' with a function in `-ldi':


AC_CHECK_LIB([di], [di_init])


# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h sys/param.h unistd.h di/di.h stdint.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([clock_gettime memset pow])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
