## Process this file with automake to produce Makefile.in

SUBDIRS = 

check_PROGRAMS = check_lib_c

check_lib_c_SOURCES = check_lib_c.c
check_lib_c_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/src/modules -I$(top_builddir)/src/modules
check_lib_c_LDADD = @LIBS_DEPENDENCIES@ @FCLIBS@

