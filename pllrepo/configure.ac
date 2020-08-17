#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([libpll], [1.0.0], [Tomas.Flouri@h-its.org])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE

# AM_MAINTAINER_MODE

LIBPLL_VERSION=1.0.0
LIBPLL_MAJOR=1
LIBPLL_MINOR=0
LIBPLL_REV=0

AC_CONFIG_SRCDIR([src/pll.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_CPP
LT_INIT

AM_PROG_CC_C_O
# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h float.h limits.h malloc.h stddef.h stdint.h stdlib.h string.h strings.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_FUNC_REALLOC
AC_CHECK_FUNCS([atexit clock_gettime getcwd getpagesize gettimeofday memmove memset munmap pow sqrt strcasecmp strchr strdup strndup strstr])

have_pthreads=no
AC_SEARCH_LIBS([pthread_create], [pthread], [have_pthreads=yes])

if test "x${have_pthreads}" = xyes; then
  AC_CHECK_HEADERS([pthread.h], [], [have_pthreads=no])
fi


have_generic=yes
have_avx=no
have_avx_pthreads=no
have_avx_mpi=no
have_sse3=no
have_sse3_pthreads=no
have_sse3_mpi=no

AX_EXT

if test "x${ax_cv_have_avx_ext}" = "xyes"; then
  have_avx=yes
  if test "x${have_pthreads}"="xyes"; then
    have_avx_pthreads=yes
  fi
fi

if test "x${ax_cv_have_sse3_ext}" = "xyes"; then
  have_sse3=yes
  if test "x${have_pthreads}"="xyes"; then
    have_sse3_pthreads=yes
  fi
fi

AC_ARG_ENABLE(generic,AC_HELP_STRING([--enable-generic],[build generic version of the library]),[have_generic=yes])
AC_ARG_ENABLE(generic,AC_HELP_STRING([--disable-generic],[build generic version of the library]),[have_generic=no])
AC_ARG_ENABLE(avx,AC_HELP_STRING([--enable-avx],[build AVX version of the library]),[have_avx=yes])
AC_ARG_ENABLE(avx,AC_HELP_STRING([--disable-avx],[build AVX version of the library]),[have_avx=no])
AC_ARG_ENABLE(avx-pthreads,AC_HELP_STRING([--enable-avx-pthreads],[build pthreads AVX version of the library]),[have_avx_pthreads=yes])
AC_ARG_ENABLE(avx-pthreads,AC_HELP_STRING([--disable-avx-pthreads],[build pthreads AVX version of the library]),[have_avx_pthreads=no])
AC_ARG_ENABLE(avx-mpi,AC_HELP_STRING([--enable-avx-mpi],[build MPI AVX version of the library]),[have_avx_mpi=yes])
AC_ARG_ENABLE(avx-mpi,AC_HELP_STRING([--disable-avx-mpi],[build MPI AVX version of the library]),[have_avx_mpi=no])
AC_ARG_ENABLE(sse3,AC_HELP_STRING([--enable-sse3],[build SSE3 version of the library]),[have_sse3=yes])
AC_ARG_ENABLE(sse3,AC_HELP_STRING([--disable-sse3],[build SSE3 version of the library]),[have_sse3=no])
AC_ARG_ENABLE(sse3-pthreads,AC_HELP_STRING([--enable-sse3-pthreads],[build pthreads SSE3 version of the library]),[have_sse3_pthreads=yes])
AC_ARG_ENABLE(sse3-pthreads,AC_HELP_STRING([--disable-sse3-pthreads],[build pthreads SSE3 version of the library]),[have_sse3_pthreads=no])
AC_ARG_ENABLE(sse3-mpi,AC_HELP_STRING([--enable-sse3-mpi],[build MPI SSE3 version of the library]),[have_sse3_mpi=yes])
AC_ARG_ENABLE(sse3-mpi,AC_HELP_STRING([--disable-sse3-mpi],[build MPI SSE3 version of the library]),[have_sse3_mpi=no])

AX_MPI([have_mpi=yes],[have_mpi=no])

AM_INIT_AUTOMAKE([subdir-objects])
AM_CONDITIONAL(BUILD_GENERIC, test "x${have_generic}" = "xyes")
AM_CONDITIONAL(BUILD_AVX, test "x${have_avx}" = "xyes")
AM_CONDITIONAL(BUILD_AVX_PTHREADS, test "x${have_avx_pthreads}" = "xyes")
AM_CONDITIONAL(BUILD_SSE3, test "x${have_sse3}" = "xyes")
AM_CONDITIONAL(BUILD_SSE3_PTHREADS, test "x${have_sse3_pthreads}" = "xyes")
AM_CONDITIONAL(BUILD_MPI, test "x${have_mpi}" = "xyes")

AC_SUBST(LIBPLL_VERSION)
AC_SUBST(LIBPLL_MAJOR)
AC_SUBST(LIBPLL_MINOR)
AC_SUBST(LIBPLL_REV)
AC_SUBST(MPICC)



AC_CONFIG_FILES([Makefile
                 MPI/Makefile
                 src/Makefile
                 man/Makefile
                 examples/Makefile])
AC_OUTPUT
