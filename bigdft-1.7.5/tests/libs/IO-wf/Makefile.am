## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.kpt \
	write.dft \
	read.dft \
	posinp.xyz \
	write.out.ref \
	read.out.ref \
        check_lib.out.ref

PPFILES = psppar.O

if HAVE_ETSF_IO
#etsf_tests = write_etsf.report read_etsf.report exportwf_etsf.report check_lib_etsf.report
etsf_tests = check_lib_etsf.report
endif

if BUILD_BINDINGS
bindings_tests = check_lib_bindings.report
bindings_bin = check_lib_c
endif

# Checks to do.
#check: write_plain.report write_binary.report \
#	read_plain.report read_binary.report \
#	exportwf.report check_lib_plain.report check_lib_binary.report \
#	$(etsf_tests) $(bindings_tests)

check: check_lib_plain.report check_lib_binary.report $(etsf_tests) $(bindings_tests)

report:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

check_PROGRAMS = check_lib $(bindings_bin)

AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
else
LIBS += @LIBS_SHORT_DEPS@
endif

AM_FCFLAGS = -I$(top_builddir)/includes -I$(top_builddir)/src/PSolver -I$(top_srcdir)/src/PSolver -I. @LIBABINIT_INCLUDE@ @LIB_XC_CFLAGS@  @MPI_INCLUDE@

check_lib_SOURCES = check_lib.f90
check_lib_FCFLAGS = -I$(top_builddir)/src/ -I$(top_builddir)/includes  -I$(top_builddir)/src/modules @LIBABINIT_INCLUDE@
check_lib_LDADD = -lbigdft-1

check_lib_c_SOURCES = check_lib_c.c
check_lib_c_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/src/modules -I$(top_builddir)/src/modules
check_lib_c_LDADD = -lbigdft-1 @FCLIBS@

clean-local: 
	rm -rf data-*

distclean: clean clean-local
	rm -rf Makefile

CLEANFILES = write_*.out read_*.out *.prc *.report geopt.mon default.* \
	wf_*.cube exportwf_*.log binary.* plain.* etsf.* forces_* \
	check_lib_*.out $(EXTRA_CLEAN) $(PPFILES)


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

if BUILD_DYNAMIC_LIBS
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/src
endif

#Check read and write by doing twice run
write_etsf.out: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz write.dft 
	$(MAKE) in
	@cp write.dft etsf.dft
	@cp input.kpt etsf.kpt
	@cp posinp.xyz etsf.xyz
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	$(run_parallel) $(top_builddir)/src/bigdft etsf > $@
	@rm -f etsf.dft etsf.kpt etsf.xyz
	$(MAKE) clean.in

write_plain.out: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz write.dft 
	$(MAKE) in
	@cp write.dft plain.dft
	@cp input.kpt plain.kpt
	@cp posinp.xyz plain.xyz
	sed -i "s/^0 3 0/0 1 0/g" plain.dft
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	$(run_parallel) $(top_builddir)/src/bigdft plain > $@
	@rm -f plain.dft plain.kpt plain.xyz
	$(MAKE) clean.in

write_binary.out: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz write.dft 
	$(MAKE) in
	@cp write.dft binary.dft
	@cp input.kpt binary.kpt
	@cp posinp.xyz binary.xyz
	sed -i "s/^0 3 0/0 2 0/g" binary.dft
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	$(run_parallel) $(top_builddir)/src/bigdft binary > $@
	@rm -f binary.dft binary.kpt binary.xyz
	$(MAKE) clean.in


read_etsf.out: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz read.dft write_etsf.out
	$(MAKE) in
	@cp read.dft etsf.dft
	@cp input.kpt etsf.kpt
	@cp posinp.xyz etsf.xyz
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	$(run_parallel) $(top_builddir)/src/bigdft etsf > $@
	@rm -f etsf.dft etsf.kpt etsf.xyz
	$(MAKE) clean.in

read_plain.out: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz read.dft write_plain.out
	$(MAKE) in
	@cp read.dft plain.dft
	@cp input.kpt plain.kpt
	@cp posinp.xyz plain.xyz
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	$(run_parallel) $(top_builddir)/src/bigdft plain > $@
	@rm -f plain.dft plain.kpt plain.xyz
	$(MAKE) clean.in

read_binary.out: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz read.dft write_binary.out
	$(MAKE) in
	@cp read.dft binary.dft
	@cp input.kpt binary.kpt
	@cp posinp.xyz binary.xyz
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	$(run_parallel) $(top_builddir)/src/bigdft binary > $@
	@rm -f binary.dft binary.kpt binary.xyz
	$(MAKE) clean.in


wf_etsf.cube: write_etsf.out $(top_builddir)/src/bigdft-tool $(PPFILES) input.kpt posinp.xyz read.dft
	$(MAKE) in
	@cp read.dft etsf.dft
	@cp input.kpt etsf.kpt
	@cp posinp.xyz etsf.xyz
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	$(top_builddir)/src/bigdft-tool --name etsf -a export-wf data-etsf/wavefunction.etsf --i-band 3 --i-spin 2 --i-kpt 2 --i-complex 2 > exportwf_etsf.log
	mv wavefunction-k002-DI.b0003.cube $@
	@rm -f wavefunction-k002-DI.b0003_avg_*
	@rm -f etsf.dft etsf.kpt etsf.xyz
	$(MAKE) clean.in

wf_plain.cube: write_plain.out $(top_builddir)/src/bigdft-tool $(PPFILES) input.kpt posinp.xyz read.dft
	$(MAKE) in
	@cp read.dft plain.dft
	@cp input.kpt plain.kpt
	@cp posinp.xyz plain.xyz
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	$(top_builddir)/src/bigdft-tool --name plain -a export-wf data-plain/wavefunction-k002-DI.b0003 > exportwf_plain.log
	mv wavefunction-k002-DI.b0003.cube $@
	@rm -f wavefunction-k002-DI.b0003_avg_*
	@rm -f plain.dft plain.kpt plain.xyz
	$(MAKE) clean.in

wf_binary.cube: write_binary.out $(top_builddir)/src/bigdft-tool $(PPFILES) input.kpt posinp.xyz read.dft
	$(MAKE) in
	@cp read.dft binary.dft
	@cp input.kpt binary.kpt
	@cp posinp.xyz binary.xyz
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	$(top_builddir)/src/bigdft-tool --name binary -a export-wf data-binary/wavefunction-k002-DI.bin.b0003 > exportwf_binary.log
	mv wavefunction-k002-DI.b0003.cube $@
	@rm -f wavefunction-k002-DI.b0003_avg_*
	@rm -f binary.dft binary.kpt binary.xyz
	$(MAKE) clean.in

check_lib_plain.out: check_lib write_plain.out
	./check_lib data-plain/wavefunction-k002-DI.b0003 > $@

check_lib_binary.out: check_lib write_binary.out
	./check_lib data-binary/wavefunction-k002-DI.bin.b0003 > $@

check_lib_etsf.out: check_lib write_etsf.out
	./check_lib data-etsf/wavefunction.etsf > $@

check_lib_bindings.out: check_lib_c write_binary.out
	./check_lib_c data-binary/wavefunction-k002-DI.bin.b0003 > $@

# Write tests.
write_plain.report: write_plain.out $(srcdir)/write.out.ref
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft $^ | tee $@

write_binary.report: write_binary.out $(srcdir)/write.out.ref
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft $^ | tee $@

write_etsf.report: write_etsf.out $(srcdir)/write.out.ref
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft $^ | tee $@

# Read tests.
read_plain.report: read_plain.out $(srcdir)/read.out.ref
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft $^ | tee $@

read_binary.report: read_binary.out $(srcdir)/read.out.ref
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft $^ | tee $@

read_etsf.report: read_etsf.out $(srcdir)/read.out.ref
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft $^ | tee $@

# Bigdft-tool exportation
exportwf.report: wf_binary.cube wf_plain.cube
	python $(top_srcdir)/tests/fldiff.py --discrepancy=2.e-9 $^ | tee $@

exportwf_etsf.report: wf_binary.cube wf_etsf.cube
	python $(top_srcdir)/tests/fldiff.py --discrepancy=2.e-9 $^ | tee $@

# Library in Fortran
check_lib_plain.report: $(srcdir)/check_lib.out.ref check_lib_plain.out
	python $(top_srcdir)/tests/fldiff.py $^ | tee $@

check_lib_binary.report: $(srcdir)/check_lib.out.ref check_lib_binary.out
	python $(top_srcdir)/tests/fldiff.py $^ | tee $@

check_lib_etsf.report: $(srcdir)/check_lib.out.ref check_lib_etsf.out
	python $(top_srcdir)/tests/fldiff.py $^ | tee $@

# C bindings
check_lib_bindings.report: $(srcdir)/check_lib.out.ref check_lib_bindings.out
	python $(top_srcdir)/tests/fldiff.py $^ | tee $@
