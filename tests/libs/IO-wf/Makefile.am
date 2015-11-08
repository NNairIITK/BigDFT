## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	input.kpt \
	write.dft \
	read.dft \
	posinp.xyz \
	write.ref.yaml \
	read.ref.yaml \
    check_lib_plain.ref.yaml

PPFILES = psppar.O

if HAVE_ETSF_IO
#etsf_tests = write_etsf.report read_etsf.report exportwf_etsf.report check_lib_etsf.report
etsf_tests = check_lib_etsf.report.yaml
endif

if BUILD_BINDINGS
bindings_tests = check_lib_bindings.report.yaml
bindings_bin = check_lib_c
endif

# Checks to do.
#check: write_plain.report write_binary.report \
#	read_plain.report read_binary.report \
#	exportwf.report check_lib_plain.report check_lib_binary.report \
#	$(etsf_tests) $(bindings_tests)

check: check_lib_plain.report.yaml check_lib_binary.report.yaml $(etsf_tests) $(bindings_tests)

report:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

check_PROGRAMS = check_lib $(bindings_bin)

AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
else
LIBS += @LIBS_SHORT_DEPS@
endif

AM_FCFLAGS = -I$(top_builddir)/includes -I$(top_builddir)/src/PSolver -I$(top_srcdir)/src/PSolver -I. @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@  @MPI_INCLUDE@

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
write_etsf.out.yaml: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz write.dft 
	$(MAKE) in
	@cp write.dft etsf.dft
	@cp input.kpt etsf.kpt
	@cp posinp.xyz etsf.xyz
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	$(run_parallel) $(top_builddir)/src/bigdft etsf > $@
	@rm -f etsf.dft etsf.kpt etsf.xyz
	$(MAKE) clean.in

write_plain.out.yaml: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz write.dft 
	$(MAKE) in
	@cp write.dft plain.dft
	@cp input.kpt plain.kpt
	@cp posinp.xyz plain.xyz
	sed -i "s/^0 3 0/0 1 0/g" plain.dft
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	$(run_parallel) $(top_builddir)/src/bigdft plain > $@
	@rm -f plain.dft plain.kpt plain.xyz
	$(MAKE) clean.in

write_binary.out.yaml: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz write.dft 
	$(MAKE) in
	@cp write.dft binary.dft
	@cp input.kpt binary.kpt
	@cp posinp.xyz binary.xyz
	sed -i "s/^0 3 0/0 2 0/g" binary.dft
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	$(run_parallel) $(top_builddir)/src/bigdft binary > $@
	@rm -f binary.dft binary.kpt binary.xyz
	$(MAKE) clean.in


read_etsf.out.yaml: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz read.dft write_etsf.out
	$(MAKE) in
	@cp read.dft etsf.dft
	@cp input.kpt etsf.kpt
	@cp posinp.xyz etsf.xyz
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	$(run_parallel) $(top_builddir)/src/bigdft etsf > $@
	@rm -f etsf.dft etsf.kpt etsf.xyz
	$(MAKE) clean.in

read_plain.out.yaml: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz read.dft write_plain.out
	$(MAKE) in
	@cp read.dft plain.dft
	@cp input.kpt plain.kpt
	@cp posinp.xyz plain.xyz
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	$(run_parallel) $(top_builddir)/src/bigdft plain > $@
	@rm -f plain.dft plain.kpt plain.xyz
	$(MAKE) clean.in

read_binary.out.yaml: $(top_builddir)/src/bigdft $(PPFILES) input.kpt posinp.xyz read.dft write_binary.out
	$(MAKE) in
	@cp read.dft binary.dft
	@cp input.kpt binary.kpt
	@cp posinp.xyz binary.xyz
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	$(run_parallel) $(top_builddir)/src/bigdft binary > $@
	@rm -f binary.dft binary.kpt binary.xyz
	$(MAKE) clean.in


wf_etsf.cube: write_etsf.out.yaml $(top_builddir)/src/bigdft-tool $(PPFILES) input.kpt posinp.xyz read.dft
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

wf_plain.cube: write_plain.out.yaml $(top_builddir)/src/bigdft-tool $(PPFILES) input.kpt posinp.xyz read.dft
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

wf_binary.cube: write_binary.out.yaml $(top_builddir)/src/bigdft-tool $(PPFILES) input.kpt posinp.xyz read.dft
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

check_lib_plain.out.yaml: check_lib write_plain.out.yaml
	./check_lib data-plain/wavefunction-k002-DI.b0003 > $@

check_lib_binary.out.yaml: check_lib write_binary.out.yaml
	./check_lib data-binary/wavefunction-k002-DI.bin.b0003 > $@

check_lib_etsf.out.yaml: check_lib write_etsf.out.yaml
	./check_lib data-etsf/wavefunction.etsf > $@

check_lib_bindings.out.yaml: check_lib_c write_binary.out.yaml
	./check_lib_c data-binary/wavefunction-k002-DI.bin.b0003 > $@

# Write tests.
write_plain.report.yaml: write_plain.out.yaml
	python $(top_srcdir)/tests/fldiff_yaml.py -t $(top_srcdir)/tests/tols-BigDFT.yaml -r $(srcdir)/write.ref.yaml -d $< -o $@

write_binary.report.yaml: write_binary.out.yaml
	python $(top_srcdir)/tests/fldiff_yaml.py -t $(top_srcdir)/tests/tols-BigDFT.yaml -r $(srcdir)/write_binary.ref.yaml -d $< -o $@

write_etsf.report.yaml: write_etsf.out.yaml
	python $(top_srcdir)/tests/fldiff_yaml.py -t $(top_srcdir)/tests/tols-BigDFT.yaml -r $(srcdir)/write_etsf.ref.yaml -d $< -o $@

# Read tests.
read_plain.report.yaml: read_plain.out.yaml
	python $(top_srcdir)/tests/fldiff_yaml.py -t $(top_srcdir)/tests/tols-BigDFT.yaml -r $(srcdir)/read_plain.ref.yaml -d $< -o $@

read_binary.report.yaml: read_binary.out.yaml
	python $(top_srcdir)/tests/fldiff_yaml.py -t $(top_srcdir)/tests/tols-BigDFT.yaml -r $(srcdir)/read_binary.ref.yaml -d $< -o $@

read_etsf.report.yaml: read_etsf.out.yaml
	python $(top_srcdir)/tests/fldiff_yaml.py -t $(top_srcdir)/tests/tols-BigDFT.yaml -r $(srcdir)/read_etsf.ref.yaml -d $< -o $@

# Bigdft-tool exportation
exportwf.report: wf_binary.cube wf_plain.cube
	python $(top_srcdir)/tests/fldiff.py --discrepancy=2.e-9 $^ | tee $@

exportwf_etsf.report: wf_binary.cube wf_etsf.cube
	python $(top_srcdir)/tests/fldiff.py --discrepancy=2.e-9 $^ | tee $@

# Library in Fortran
check_lib_plain.report.yaml: check_lib_plain.out.yaml
	python $(top_srcdir)/tests/fldiff_yaml.py -t $(top_srcdir)/tests/tols-BigDFT.yaml -r $(srcdir)/check_lib_plain.ref.yaml -d $< -o $@

check_lib_binary.report.yaml: check_lib_binary.out.yaml
	python $(top_srcdir)/tests/fldiff_yaml.py -t $(top_srcdir)/tests/tols-BigDFT.yaml -r $(srcdir)/check_lib_binary.ref.yaml -d $< -o $@

check_lib_etsf.report.yaml: check_lib_etsf.out.yaml
	python $(top_srcdir)/tests/fldiff_yaml.py -t $(top_srcdir)/tests/tols-BigDFT.yaml -r $(srcdir)/check_lib_etsf.ref.yaml -d $< -o $@

# C bindings
check_lib_bindings.report.yaml: check_lib_bindings.out.yaml
	python $(top_srcdir)/tests/fldiff_yaml.py -t $(top_srcdir)/tests/tols-BigDFT.yaml -r $(srcdir)/check_lib_bindings.ref.yaml -d $< -o $@
