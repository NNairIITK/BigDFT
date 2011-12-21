## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.kpt \
	write.dft \
	read.dft \
	posinp.xyz \
	write.out.ref \
	read.out.ref

EXTRA_CLEAN = wavefunction.etsf

PPFILES = psppar.O

CLEANFILES = write.out write_plain.out read.out *.prc fldiff.*.report geopt.mon default.* input.dft \
	wf_etsf.cube wf_plain.cube exportwf.etsf.log exportwf.plain.log $(EXTRA_CLEAN) $(PPFILES)

check_PROGRAMS = check_lib check_lib_c

check_lib_SOURCES = check_lib.f90
check_lib_FCFLAGS = -I$(top_builddir)/src/modules
check_lib_LDADD = @LIBS_DEPENDENCIES@

check_lib_c_SOURCES = check_lib_c.c
check_lib_c_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/src/modules
check_lib_c_LDADD = @LIBS_DEPENDENCIES@ @FCLIBS@

clean-local: 
	rm -rf data

distclean: clean clean-local
	rm -rf Makefile

#dependencies
$(PPFILES):
	@cp $(top_srcdir)/utils/PSPfiles/HGH-K/$@ .

in:
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    cp -p $(srcdir)/$$i . ; \
	  done ; \
        fi

clean.in:
	@if ! test x"$(srcdir)" = x"." ; then \
	  for i in $(EXTRA_DIST) ; do \
	    rm -f $$i ; \
	  done ; \
        fi

#Check read and write by doing twice run
write.out: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz write.dft 
	$(MAKE) in
	@cp write.dft input.dft
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	@rm -f input.dft
	$(MAKE) clean.in

write_plain.out: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz write.dft 
	$(MAKE) in
	@cp write.dft input.dft
	sed -i "s/^0 3 0/0 1 0/g" input.dft
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	@rm -f input.dft
	$(MAKE) clean.in

wf_etsf.cube: write.out $(top_builddir)/src/bigdft-tool $(PPFILES) input.kpt posinp.xyz read.dft
	$(MAKE) in
	@cp read.dft input.dft
	$(top_builddir)/src/bigdft-tool -a export-wf data/wavefunction.etsf --i-band 3 --i-spin 2 --i-kpt 2 --i-complex 2 > exportwf.etsf.log
	mv wavefunction-k002-DI.b0003.cube $@
	@rm -f wavefunction-k002-DI.b0003_avg_*
	@rm -f input.dft
	$(MAKE) clean.in

wf_plain.cube: write_plain.out $(top_builddir)/src/bigdft-tool $(PPFILES) input.kpt posinp.xyz read.dft
	$(MAKE) in
	@cp read.dft input.dft
	$(top_builddir)/src/bigdft-tool -a export-wf data/wavefunction-k002-DI.b0003 > exportwf.plain.log
	mv wavefunction-k002-DI.b0003.cube $@
	@rm -f wavefunction-k002-DI.b0003_avg_*
	@rm -f input.dft
	$(MAKE) clean.in

read.out: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz read.dft 
	$(MAKE) in
	@cp read.dft input.dft
	$(run_parallel) $(top_builddir)/src/bigdft > $@
	@rm -f input.dft
	$(MAKE) clean.in

check_lib.out: check_lib
	./check_lib > $@

fldiff.write.report: write.out write.out.ref
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft write.out \
		$(srcdir)/write.out.ref | tee $@

fldiff.read.report: read.out read.out.ref
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft read.out \
		$(srcdir)/read.out.ref | tee $@

fldiff.exportwf.report: wf_etsf.cube wf_plain.cube
	python $(top_srcdir)/tests/fldiff.py wf_etsf.cube \
		wf_plain.cube | tee $@

fldiff.check_lib.report: check_lib.out.ref check_lib.out
	python $(top_srcdir)/tests/fldiff.py check_lib.out \
		$(srcdir)/check_lib.out.ref | tee $@

check: fldiff.write.report fldiff.read.report fldiff.exportwf.report fldiff.check_lib.report
