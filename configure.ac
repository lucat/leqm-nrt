AC_INIT([leqm-nrt], [0.18], [luca.trisciani@dcp-werkstatt.ch])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CC
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
	Makefile
	src/Makefile
])
AC_CHECK_HEADERS([stdlib.h nistd.h stdio.h math.h sndfile.h unistd.h pthread.h string.h time.h ctype.h iso646.h])

AC_FUNC_MALLOC
AC_CHECK_FUNCS([clock_gettime memset powf])

AC_OUTPUT
