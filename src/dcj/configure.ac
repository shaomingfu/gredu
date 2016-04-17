#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(dcj, 1.0, shaomingfu@gmail.com)
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/gene.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

AS_IF([test "s$BOOST_HOME" = "s"],[AC_MSG_ERROR([*** Please set environment variable BOOST_HOME ***])]) 
AS_IF([test "s$GUROBI_HOME" = "s"],[AC_MSG_ERROR([*** Please set environment variable GUROBI_HOME ***])]) 
AS_IF([test "s$GRB_LICENSE_FILE" = "s"],[AC_MSG_ERROR([*** Please set environment variable GRB_LICENSE_FILE ***])]) 

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_FUNCS([pow])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
