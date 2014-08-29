# Generic part, of the testing Makefiles.
# Possible calls are:
#  make in: generate all input dirs.
#  make failed-check: run check again on all directories with missing report
#                     or failed report.
#  make X.in: generate input dir for directory X.
#  make X.check: generate a report for directory X (if not already existing).
#  make X.recheck: force the creation of the report in directory X.
#  make X.clean: clean the given directory X.
#  make X.diff: make the difference between the output and the reference (with DIFF envvar)
#  make X.updateref: update the reference with the output (prompt the overwrite)

if USE_MPI
  mpirun_message = mpirun
else
  mpirun_message =
endif
if USE_OCL
oclrun_message = oclrun
accel_in_message = in_message
else
oclrun_message =
accel_in_message =
endif

if BUILD_LIBYAML
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/yaml-0.1.4/src/.libs
PYTHONPATH := ${PYTHONPATH}:`ls -d $(abs_top_builddir)/PyYAML-3.10/build/lib.*`
endif

if BUILD_DYNAMIC_LIBS
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/src
endif

AM_FCFLAGS = -I$(top_builddir)/includes -I$(top_srcdir)/PSolver/src -I. @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@  @MPI_INCLUDE@

PSPS = psppar.H \
       psppar.C \
       psppar.Li \
       psppar.Ca \
       psppar.Mn \
       psppar.N \
       psppar.Si \
       HGH/psppar.H \
       HGH/psppar.Na \
       HGH/psppar.Cl \
       HGH/psppar.O \
       HGH/psppar.Si \
       HGH/psppar.Fe \
       HGH/psppar.Mg \
       HGH/psppar.Ag \
       HGH/psppar.N \
       HGH/psppar.C \
       HGH-K/psppar.H \
       HGH-K/psppar.Si \
       HGH-K/psppar.N \
       HGH-K/psppar.O \
       HGH-K/psppar.Ti \
       extra/psppar.H \
       Xabs/psppar.Fe

ALLDIRS = $(EXTRA_TESTDIRS) $(TESTDIRS)

INS = $(ALLDIRS:=.in)
RUNS = $(ALLDIRS:=.run)
CHECKS = $(TESTDIRS:=.check) $(TESTDIRS:=.yaml-check)
EXTRA_CHECKS = $(EXTRA_TESTDIRS:=.check) $(EXTRA_TESTDIRS:=.yaml-check)
DIFFS = $(ALLDIRS:=.diff)
UPDATES = $(ALLDIRS:=.updateref)
FAILEDCHECKS = $(TESTDIRS:=.recheck)
CLEANS = $(ALLDIRS:=.clean)

EXTRA_DIST += README $(ALLDIRS)

in: $(INS)

check: $(CHECKS) report

complete-check: $(EXTRA_CHECKS) check

diff: $(DIFFS)

update-references: $(UPDATES)

clean: $(CLEANS)

distclean: $(CLEANS)
	rm -rf Makefile

failed-check: $(FAILEDCHECKS) report

report:
	@if test $(MAKELEVEL) = 0 ; then	export PYTHONPATH=${PYTHONPATH}; export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ;python $(top_srcdir)/tests/report.py ; fi

#Binary dependencies
$(abs_top_builddir)/src/BigDFT2Wannier:
	cd $(abs_top_builddir)/src && $(MAKE) BigDFT2Wannier WaCo;

%.memguess.out: $(abs_top_builddir)/src/memguess $(abs_top_builddir)/src/bigdft-tool
	@name=`basename $@ .memguess.out | sed "s/[^_]*_\?\(.*\)$$/\1/"` ; \
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
	echo "$(run_serial) $(abs_top_builddir)/src/bigdft-tool -n 1 > $@"; \
	$(run_serial) $(abs_top_builddir)/src/bigdft-tool -n 1 > $@ ; \
	mv log.yaml log-memguess.yaml ; \
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.out.out: $(abs_top_builddir)/src/bigdft
	@name=`basename $@ .out.out | sed "s/[^_]*_\?\(.*\)$$/\1/"` ; \
	if test -f list_posinp; then \
	   name=`echo '--runs-file=list_posinp --taskgroup-size=1'`; \
	fi; \
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
	echo "Running $(run_parallel) $(abs_top_builddir)/src/bigdft $$name > $@" ; \
	$(run_parallel) $(abs_top_builddir)/src/bigdft $$name > $@ ; \
	if test -f list_posinp; then cat `tail -n $$(($$(wc -l < list_posinp) - 1)) list_posinp | sed "s/^\(.*\)$$/log-\1.yaml/g"` > log.yaml ; fi ; \
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.geopt.mon.out: $(abs_top_builddir)/src/bigdft
	name=`basename $@ .geopt.mon.out | sed "s/[^_]*_\?\(.*\)$$/\1/"` ; \
	if test -n "$$name" ; then datadir="data-"$$name ; else datadir="data" ; fi ; \
	$(MAKE) -f ../Makefile $*.out.out && cp $$datadir/geopt.mon $@
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.bader.out: $(abs_top_builddir)/src/tools/bader/bader %.out.out
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
	$(abs_top_builddir)/src/tools/bader/bader data/electronic_density.cube > $@ && mv dipole.yaml log-dipole.yaml
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.freq.out: $(abs_top_builddir)/src/frequencies
	@name=`basename $@ .freq.out | sed "s/[^_]*_\?\(.*\)$$/\1/"` ; \
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
	echo "Running $(run_parallel) $(abs_top_builddir)/src/frequencies > $@" ; \
	$(run_parallel) $(abs_top_builddir)/src/frequencies > $@
	name=`basename $@ .freq.out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.NEB.out: $(abs_top_builddir)/src/NEB NEB_include.sh NEB_driver.sh
	rm -f neb.it*
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
	$(run_serial) $(abs_top_builddir)/src/NEB | tee $@
	cat neb.NEB.0*/log.yaml | grep -v "Unable to read mpd.hosts" > log.yaml
	echo "---" >> log.yaml
	grep ":" NEB.NEB.out | grep -v "<BigDFT>" >> log.yaml
	grep -v ":" NEB.NEB.out > tmp-neb.out
	mv tmp-neb.out NEB.NEB.out
	rm -rf neb.NEB.0*
	rm -f gen_output_file velocities_file
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.splsad.out: $(abs_top_builddir)/src/splsad
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
	$(run_parallel) $(abs_top_builddir)/src/splsad > $@
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.minhop.out: $(abs_top_builddir)/src/global
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
	$(run_parallel) $(abs_top_builddir)/src/global > $@
#	mv log-mdinput.yaml log.yaml
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.xabs.out: $(abs_top_builddir)/src/abscalc
	name=`basename $@ .xabs.out` ; \
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
	$(run_serial) $(abs_top_builddir)/src/abscalc $$name > $@
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.b2w.out: $(abs_top_builddir)/src/BigDFT2Wannier
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
	$(run_parallel) $(abs_top_builddir)/src/bigdft > $@
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
	$(run_parallel) $(abs_top_builddir)/src/BigDFT2Wannier > $@
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.testforces.out: $(abs_top_builddir)/src/test_forces
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
	$(run_parallel) $(abs_top_builddir)/src/test_forces > $@
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"

$(PSPS):
	cp $(abs_top_srcdir)/utils/PSPfiles/$@ .

%.clean:
	@name=`basename $@ .clean` ; dir=$$name-test ; \
	rm -rf $$dir ; rm -f $$name.* ; \
    echo "Clean $$dir"

%.post-in: ;
%.psp: ;
%.post-clean: ;
%.post-out: ;

in_message:
	@if test -n "$(run_ocl)" ; then \
	  echo "==============================================" ; \
	  echo "Will generate a 'input.perf' file to force OCL" ; \
	    if test -n "$(ocl_platform)" ; then \
	      echo "Forcing use of $(ocl_platform)" ; \
	    fi ; \
	    if test -n "$(ocl_devices)" ; then \
	      echo "Forcing use of $(ocl_devices)" ; \
	    fi ; \
	  echo "==============================================" ; \
	fi

$(INS): in_message
	@name=`basename $@ .in` ; dir=$$name-test ; \
	if test ! -d $$dir ; then mkdir $$dir ; fi ; \
	for i in $(srcdir)/$$name/* ; do cp -f $$i $$dir ; done ; \
	if test -n "$(accel_in_message)" -a -n "$(run_ocl)" ; then \
		if test "$(run_ocl)" = "CPU" ; then \
				echo "ACCEL OCLCPU" > $$dir/check.perf ; \
		elif test "$(run_ocl)" = "ACC" ; then \
				echo "ACCEL OCLACC" > $$dir/check.perf ; \
		else \
				echo "ACCEL OCLGPU" > $$dir/check.perf ; \
		fi ; \
		if test -n "$(ocl_platform)" ; then \
				echo "OCL_PLATFORM $(ocl_platform)" >> $$dir/check.perf ; \
		fi ; \
		if test -n "$(ocl_devices)" ; then \
				echo "OCL_DEVICES $(ocl_devices)" >> $$dir/check.perf ; \
		fi ; \
	else echo -n "" > $$dir/check.perf ; fi ; \
	echo "outdir ./" >> $$dir/check.perf ; \
	chmod u+w $$dir/* ; \
	for i in $$dir/*.out.ref.yaml ; do \
	    base=`basename $$i .out.ref.yaml | sed "s/[^_]*_\?\(.*\)$$/\1/"` ; \
	    if test -n "$$base" ; then cat $$dir/check.perf >> $$dir/$$base.perf ; fi ; \
	done ; \
	cat $$dir/check.perf >> $$dir/input.perf ; \
    cd $$dir && $(MAKE) -f ../Makefile $$name".psp"; \
    $(MAKE) -f ../Makefile $$dir".post-in"; \
    echo "Input prepared in \"$$dir\" directory, make $$name.run available"
	touch $@


run_message:
	@if test -n "$(run_parallel)" ; then \
	  echo "==============================================" ; \
	  echo "Will run tests in parallel with '$$run_parallel'" ; \
	  echo "==============================================" ; \
	fi

%.run: %.in run_message
	@name=`basename $@ .run` ; dir=$$name-test ; \
        runs="$(srcdir)/$$name/*.ref" ; \
	tgts=`for r in $$runs ; do echo $$(basename $$r .ref)".out"; done` ; \
        cd $$dir && $(MAKE) -f ../Makefile $$tgts ; \
        echo "Tests have run in \"$$dir\" directory, make $$name.check available"
	@touch $@

%.check: %.run %.yaml-check
	@name=`basename $@ .check` ; dir=$$name-test ; \
        chks="$(srcdir)/$$name/*.ref" ; \
	tgts=`for c in $$chks ; do echo $$(basename $$c .ref)".report"; done` ; \
        cd $$dir && $(MAKE) -f ../Makefile $$tgts
	touch $@

%.yaml-check: %.run
	@name=`basename $@ .yaml-check` ; dir=$$name-test ; \
        chks="$(srcdir)/$$name/*.ref.yaml" ; \
	tgts=`for c in $$chks ; do echo $$(basename $$c .ref.yaml)".report.yaml"; done` ; \
        cd $$dir && $(MAKE) -f ../Makefile $$tgts
	touch $@


%.diff	: %.run
	@if test -z "$$DIFF" ; then echo "The environment variable DIFF is missing!"; else \
			name=`basename $@ .diff` ; dir=$$name-test ; \
				chks="$(srcdir)/$$name/*.ref" ; \
			for c in $$chks ; do \
			    echo "$$DIFF $$c $$dir/$$(basename $$c .ref).out"; \
				$$DIFF $$c $$dir/$$(basename $$c .ref)".out"; \
			done ; \
				ychks="$(srcdir)/$$name/*.ref.yaml" ; \
			for c in $$ychks ; do base=`basename $$c .out.ref.yaml | sed s/.out// | sed s/.xabs// | sed "s/[^_]*_\?\(.*\)$$/\1/"` ; \
			if test -n "$$base" ; then \
			echo "$$DIFF $$c $$dir/log-$$base.yaml" ; \
			$$DIFF $$c $$dir/log-$$base.yaml; \
			else \
			echo "$$DIFF $$c $$dir/log.yaml" ; \
			$$DIFF $$c $$dir/log.yaml; \
			fi ;\
			done ; \
    fi ; \
	touch $@

%.updateref: #%.run %.diff
	@name=`basename $@ .updateref` ; dir=$$name-test ; \
        chks="$(srcdir)/$$name/*.ref" ; \
	for c in $$chks ; do echo "Update reference with " $$dir/$$(basename $$c .ref)".out"; \
	                     cp -vi $$dir/$$(basename $$c .ref)".out"  $$c;\
	done ; \
        ychks="$(srcdir)/$$name/*.ref.yaml" ; \
	for c in $$ychks ; do base=`basename $$c .out.ref.yaml | sed s/.out// | sed s/.xabs// | sed "s/[^_]*_\?\(.*\)$$/\1/"`  ;\
	if test -n "$$base" ; then \
	echo "Update reference with " $$dir/log-$$base.yaml; \
	                     cp -vi $$dir/log-$$base.yaml $$c;\
	else \
	echo "Update reference with " $$dir/log.yaml; \
	                     cp -vi $$dir/log.yaml $$c;\
	fi ;\
	done ; \
	touch $@

%.recheck: %.in
	@name=`basename $@ .recheck` ; dir=$$name-test ; \
        refs="$$dir/*.ref" ; \
	for r in $$refs ; do \
	  rep=`basename $$r .ref`".report" ; \
	  if ! grep -qs "succeeded\|passed" $$dir/$$rep ; then \
	    target=` basename $$r .ref` ; \
	    rm -f $$dir/$$target".out" $$dir/$$target".report" ; \
	    cd $$dir && $(MAKE) -f ../Makefile $$target".out" $$target".report" && cd - ; \
	  fi \
	done
	touch $*".check"

# Avoid copying in dist the builddir files.
distdir: $(DISTFILES)
	@srcdirstrip=`echo "$(srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	topsrcdirstrip=`echo "$(top_srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	list='$(DISTFILES)'; \
	  dist_files=`for file in $$list; do echo $$file; done | \
	  sed -e "s|^$$srcdirstrip/||;t" \
	      -e "s|^$$topsrcdirstrip/|$(top_builddir)/|;t"`; \
	case $$dist_files in \
	  */*) $(MKDIR_P) `echo "$$dist_files" | \
			   sed '/\//!d;s|^|$(distdir)/|;s,/[^/]*$$,,' | \
			   sort -u` ;; \
	esac; \
	for file in $$dist_files; do \
	  d=$(srcdir); \
	  if test -d $$d/$$file; then \
	    dir=`echo "/$$file" | sed -e 's,/[^/]*$$,,'`; \
	    if test -d "$(distdir)/$$file"; then \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    if test -d $(srcdir)/$$file && test $$d != $(srcdir); then \
	      cp -fpR $(srcdir)/$$file "$(distdir)$$dir" || exit 1; \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    cp -fpR $$d/$$file "$(distdir)$$dir" || exit 1; \
	  else \
	    test -f "$(distdir)/$$file" \
	    || cp -p $$d/$$file "$(distdir)/$$file" \
	    || exit 1; \
	  fi; \
	done

# Doc messages.
all:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) foot_message ; fi

head_message:
	@echo "=============================================================================="
	@echo " This is a directory for tests. Beside the 'make check'"
	@echo " one can use the following commands:"
	@echo "  make in:           generate all input dirs."
	@echo "  make failed-check: run check again on all directories"
	@echo "                     with missing report or failed report."
	@echo "  make complete-check: for developers, makes long and extensive tests."
	@echo "  make X.in:         generate input dir for directory X."
	@echo "  make X.check:      generate a report for directory X"
	@echo "                     (if not already existing)."
	@echo "  make X.recheck:    force the creation of the report in directory X."
	@echo "  make X.clean:      clean the given directory X."
	@echo "  make X.diff:       make the difference between output and the reference"
	@echo "                     (with the environment variable DIFF)"
	@echo "  make X.updateref   update the reference with the output"
	@echo "                     (prompt the overwrite)"	

mpirun: head_message
	@echo ""
	@echo " Use the environment variable run_parallel"
	@echo "     ex: export run_parallel='mpirun -np 2'  "

oclrun: head_message $(mpirun_message)
	@echo ""
	@echo " Use the environment variable run_ocl"
	@echo "     ex: export run_ocl='on' to use OpenCL acceleration"
	@echo "     use run_ocl='CPU' or 'ACC' to force use of hardware different than GPU"
	@echo " and the environment variables ocl_platform and ocl_devices"
	@echo "     ex: export ocl_platform='NVIDIA'"
	@echo "     ex: export ocl_devices='K20'"

foot_message: $(mpirun_message) $(oclrun_message) head_message
	@echo "=============================================================================="

