#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.69])
AC_INIT([leqm-nrt], [0.20], [luca.trisciani@dcp-werkstatt.ch])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR([src/leqm-nrt.c])


# Checks for programs.
AC_PROG_CC


# Determine target operating system
# AC_CANONICAL_HOST is needed to access the 'host_os' variable    
AC_CANONICAL_HOST

build_linux=no
build_windows=no
build_mac=no

# Detect the target system
case "${host_os}" in
    linux*)
        build_linux=yes
        ;;
    cygwin*|mingw*)
        build_windows=yes
        ;;
    darwin*)
        build_mac=yes
        ;;
    *)
        AC_MSG_ERROR(["OS $host_os is not supported"])
        ;;
esac

# Pass the conditionals to automake
AM_CONDITIONAL([LINUX], [test "$build_linux" = "yes"])
AM_CONDITIONAL([WINDOWS], [test "$build_windows" = "yes"])
AM_CONDITIONAL([OSX], [test "$build_mac" = "yes"])





#CFLAGS="$CFLAGS -I/usr/include/ffmpeg"


#CXXFLAGS="$CXXFLAGS -I/opt/local/include"

## Dolby DI



# Checks for libraries.

## FFMPEG
## Add alternative for older FFMPEG versions 


FFMPEG_libs=


# FIXME: Replace `main' with a function in `-lavcodec':
AC_CHECK_LIB([avcodec], [avcodec_send_packet, avcodec_receive_frame], [
FFMPEG_libs=$FFMPEG_libs -lavcodec
AC_CHECK_LIB([avformat], [avformat_open_input, avformat_find_stream_info], [
FFMPEG_libs=$FFMPEG_libs -lavformat
AC_CHECK_LIB([avutil], [av_get_channel_layout], [
FFMPEG_libs=$FFMPEG_libs -lavutil
])
])
])

LIBS=$FFMPEG_libs

# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [sqrt])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [pthread_attr_init])

AS_IF([test x = "x${FFMPEG_libs}"],
[AC_CHECK_LIB([sndfile], [sf_readf_double])],
[AC_MSG_RESULT([Neither ffmpeg nor sndfile were found. Please install at least one of the two libraries.])])


# FIXME: Replace `main' with a function in `-lrt':
if test "$build_linux" == "yes"
then
AC_CHECK_LIB([rt], [clock_gettime])
fi

## Dolby DI
not_found_di=0
#LDFLAGS="$LDFLAGS -L/usr/lib/di"
AC_CHECK_HEADERS([di/di.h], [], [not_found_di=1])

LDFLAGS_save=$LDFLAGS
di_libdir=
found_di_lib=0
for di_libdir in "/usr/lib" "/usr/lib/di" "/usr/local/lib" "/usr/local/lib/di"; do
    LDFLAGS="$LDFLAGS -L$di_libdir"
    AC_CHECK_LIB([di], [di_init], [
    LIBS="${LIBS} -ldi"
    found_di_lib=1
    break
    ], [
    di_libdir=
    ])
    unset ac_cv_lib_di_di_init
    done
if test $found_di_lib == 0
then
	not_found_di=1
fi


## Dolby DI
#LDFLAGS="$LDFLAGS -L/usr/lib/di"
#not_found_di=0
#AC_CHECK_HEADERS([di/di.h], [], [not_found_di=1])
#AC_SEARCH_LIBS([di_init], [di], [], [not_found_di=1])
#
if test $not_found_di == 0
then
	LDFLAGS="${LDFLAGS_save} -L${di_libdir}"
	AC_DEFINE([HAVE_LIBDI],[1],[Define if Dolby libdi found])
fi

#CFLAGS="$CFLAGS -I/usr/include/di"



#AC_CHECK_LIB([di], [di_init])



#AC_CHECK_LIB([di], [di_init])

#
# Checks for header files.
#AC_CHECK_HEADERS([stdlib.h string.h sys/param.h unistd.h di/di.h stdint.h libavformat/format.h libavcodec/avcodec.h libavutil/avutil.h libavutil/assert.h sndfile.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([clock_gettime memset pow])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
