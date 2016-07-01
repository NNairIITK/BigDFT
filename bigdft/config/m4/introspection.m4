# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_INTROSPECTION],
[dnl Test for GLib and introspection.
  AC_ARG_WITH(gobject, AS_HELP_STRING([--with-gobject],
                       [Build GObject bindings (default is no).]),
              ax_have_glib=$withval, ax_have_glib="no")
  if test x"$ax_have_glib" = x"yes"; then
    PKG_CHECK_MODULES(GLIB, glib-2.0 gobject-2.0 gthread-2.0 gio-2.0 >= 2.22, [ax_have_glib=yes], [ax_have_glib=no])
    PKG_CHECK_MODULES([GOBJECT_INTROSPECTION],
                      [gobject-introspection-1.0 >= 0.9],
                      [ax_have_introspection=yes],
                      [ax_have_introspection=no])
    if test "x$ax_have_introspection" = "xyes" ; then
      dnl girdir=$($PKG_CONFIG --variable=girdir gobject-introspection-1.0)
      dnl typelibsdir=$($PKG_CONFIG --variable=typelibdir gobject-introspection-1.0)
      girdir="$datadir/gir-1.0"
      typelibsdir="$libdir/girepository-1.0"
      AC_SUBST(girdir)
      AC_SUBST(typelibsdir)
      AC_SUBST([G_IR_SCANNER], [$($PKG_CONFIG --variable=g_ir_scanner gobject-introspection-1.0)])
      AC_SUBST([G_IR_COMPILER], [$($PKG_CONFIG --variable=g_ir_compiler gobject-introspection-1.0)])
  
      dnl Add Python support for the PythonGI plug-in.
      AC_LANG_PUSH(C)
      AM_CHECK_PYTHON_HEADERS([AC_DEFINE([HAVE_PYTHON], [], [If set, we can call Python.h])],[AC_MSG_WARN(could not find Python headers)])
      PYTHON_LIBS="`python-config --libs`"
      AC_SUBST(PYTHON_LIBS)
      AC_LANG_POP(C)
    fi
  else
    ax_have_glib=no
    ax_have_introspection=no
  fi
  if test x"$ax_have_glib" = x"yes" ; then
    AC_DEFINE([HAVE_GLIB], [], [If set, we can call glib.h])
    AC_SUBST([GLIB_TRUE], [""])
    AC_SUBST([GLIB_END_TRUE], [""])
    AC_SUBST([GLIB_FALSE], ["/*"])
    AC_SUBST([GLIB_END_FALSE], ["*/"])
  else
    AC_SUBST([GLIB_TRUE], ["/*"])
    AC_SUBST([GLIB_END_TRUE], ["*/"])
    AC_SUBST([GLIB_FALSE], [""])
    AC_SUBST([GLIB_END_FALSE], [""])
  fi

  AM_CONDITIONAL([HAVE_GLIB], [test "$ax_have_glib" = "yes"])
  AM_CONDITIONAL([WITH_GOBJECT_INTROSPECTION],
                 [test "x$enable_introspection" = "xyes"])
])
