AC_INIT([gaia], [0.1.0])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_AUX_DIR([m4])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE([foreign])
AM_CONFIG_HEADER([config.h])
LT_INIT(disable-shared)

AC_PROG_INSTALL
AC_PROG_CC
AC_PROG_CXX
AC_PROG_LN_S

AC_ARG_WITH(jansson_include,
	[  --with-jansson-include=DIR  jansson include directory ], [
	CPPFLAGS="$CPPFLAGS -I$withval"
	])
AC_ARG_WITH(jansson_lib,
	[  --with-jansson-lib=DIR      jansson lib directory ], [
	LDFLAGS="$LDFLAGS -L$withval"
	jansson_lib_path="$withval"
	])
AC_ARG_WITH(jansson,
	[  --with-jansson=DIR          jansson install directory ], [
	if test -z "$with_jansson_lib" && test -z "$with_jansson_include" ; then
		CPPFLAGS="$CPPFLAGS -I$withval/include"
		LDFLAGS="$LDFLAGS -L$withval/lib"
		jansson_lib_path="$withval/lib"
	else
		AC_MSG_FAILURE([Do not use --with-jansson and --with-jansson-include/--with-jansson-lib options simultaneously.])
	fi
	])

AC_MSG_CHECKING(for jansson.h)
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include "jansson.h"
]])],[AC_MSG_RESULT(yes)],[
	AC_MSG_RESULT(no)
	AC_MSG_ERROR([jansson.h can't be found, or is unusable.])
])

dnl Check if we can link with jansson
AC_CHECK_LIB(jansson, json_object, [LIBS="-ljansson $LIBS"],
	[AC_MSG_ERROR([libjansson not found, or is unusable.])])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
