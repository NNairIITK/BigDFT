## Process this file with automake to produce Makefile.in

SUBDIRS = 

check_PROGRAMS = check_lib_c

check_lib_c_SOURCES = check_lib_c.c
check_lib_c_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/src/modules -I$(top_builddir)/src/modules @GLIB_CFLAGS@
check_lib_c_LDADD = @LIBS_DEPENDENCIES@ @FCLIBS@

EXTRA_DIST = bindings.out.ref posinp.ascii test.dft
CLEANFILES = bindings.out bindings.report output.ascii

check: bindings.report

#check_lib_c: $(top_builddir)/src/libbigdft.a

$(EXTRA_DIST):
	@cp -p $(srcdir)/$@ .

bindings.report: $(srcdir)/bindings.out.ref bindings.out
	python $(top_srcdir)/tests/fldiff.py --discrepancy=1.5e-9 $^ | tee $@

bindings.out: check_lib_c posinp.ascii $(srcdir)/test.dft
	./check_lib_c $(srcdir) > $@
	cat output.ascii >> $@