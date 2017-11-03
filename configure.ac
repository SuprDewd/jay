
AC_INIT([jay], [0.1], [bjarki.agust@gmail.com])
AM_INIT_AUTOMAKE([subdir-objects -Wall -Werror foreign])
AC_PROG_CXX
AC_CONFIG_FILES([
    Makefile
    src/Makefile
])

CXXFLAGS="$CXXFLAGS -std=c++0x"

AC_CHECK_HEADERS([gmp.h])
AC_CHECK_HEADERS([gmpxx.h])
AC_CHECK_LIB([gmp], [__gmpz_init])
AC_CHECK_LIB([gmpxx], [main])

AC_OUTPUT

