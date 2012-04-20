## Process this file with automake to produce Makefile.in

SUBDIRS = 

check_PROGRAMS = check_lib_c check_client_signals

check_lib_c_SOURCES = check_lib_c.c
check_lib_c_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/src/modules -I$(top_builddir)/src/modules @GLIB_CFLAGS@
check_lib_c_LDADD = @LIBS_DEPENDENCIES@ @FCLIBS@

check_client_signals_SOURCES = check_client_signals.c
check_client_signals_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/src/modules -I$(top_srcdir)/src/lib -I$(top_builddir)/src/modules @GLIB_CFLAGS@
check_client_signals_LDADD = @LIBS_DEPENDENCIES@ @FCLIBS@

EXTRA_DIST = bindings.out.ref signals.out.ref posinp.ascii test.dft
CLEANFILES = bindings.out bindings.report output.ascii default.* forces_test.ascii test.perf signals.out signals.report time.yaml

check: bindings.report signals.report

#check_lib_c: $(top_builddir)/src/libbigdft.a

$(EXTRA_DIST):
	@cp -p $(srcdir)/$@ .

bindings.report: $(srcdir)/bindings.out.ref bindings.out
	python $(top_srcdir)/tests/fldiff.py --discrepancy=1.5e-9 $^ | tee $@

bindings.out: check_lib_c posinp.ascii $(srcdir)/test.dft
	./check_lib_c $(srcdir) > $@
	cat output.ascii >> $@

signals.report: $(srcdir)/signals.out.ref signals.out
	python $(top_srcdir)/tests/fldiff.py --discrepancy=1.5e-9 $^ | tee $@

signals.out: check_client_signals posinp.ascii $(srcdir)/test.dft
	@if ! test x"$(srcdir)" = x"." ; then \
	  ln -sf $(srcdir)/posinp.ascii test.ascii; \
	  cp -p $(srcdir)/test.dft .; \
        fi
	@echo "signaling T" > test.perf
	@echo "signalTimeout 15" >> test.perf
	$(run_parallel) $(top_builddir)/src/bigdft test > /dev/null &
	sleep 3s; ./check_client_signals > $@
	@if ! test x"$(srcdir)" = x"." ; then \
	  rm -f test.ascii ; \
	  rm -f test.dft ; \
        fi
