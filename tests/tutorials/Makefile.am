# Give here the list of existing tests

#SUBDIRS = Xabs/pawpatchgen  Xabs/scfpotgen   Xabs/spectra    

TESTDIRS = CH4 \
	  CH3- \
	  H2O-slab \
	  N2 \
	  H2O-bader \
	  scalability \
	  dev

EXTRA_DIST = $(TESTDIRS) \
     Figures \
     H1A-Basics.html \
     H1B-Basic_surfaces.html \
     H1C-Basic_post_processing.html \
     H2B-scalability.html \
     H2E1-troubleshooting_minima_hopping.html \
     H2E-minima_hopping.html \
     H3B-charge_analysis.html \
     H3D-Xabs.html \
     H3X-hacking_bigdft.html \
     Mg7-minhop \
     minhop.pdf \
     GPU \
     H2A-GPU.html

H2O-bader.post-in :
	cp -f input-all.dft input.dft
H2O-bader.dipole.dat.out: H2O-bader.out.out
	$(abs_top_builddir)/src/tools/bader/bader data/electronic_density.cube > bader.out &&\
	mv dipole.dat $@

dev.check: dev/toy_model

dev/toy_model: $(srcdir)/dev/toy_model.f90
	$(FC) -o $@ $(FCFLAGS) $(AM_FCFLAGS) $(srcdir)/dev/toy_model.f90 -L$(top_builddir)/src -lbigdft  -L$(top_builddir)/src/PSolver -lpoissonsolver $(libs_opt) $(LDFLAGS) $(LIBS)
toy_model.dat.out: dev-N2.out.out toy_model
	./toy_model > toy_model.dat.out

#Flags for optimization libraries
if USE_DGEMMSY
libs_dgemmsy = $(top_builddir)/src/dgemmsy/libdgemmsy.a
else
libs_dgemmsy =
endif
if USE_CUDA_GPU
libs_cuda=$(top_builddir)/src/CUDA/libGPU.a
else
libs_cuda=
endif
if USE_OCL
libs_ocl = $(top_builddir)/src/OpenCL/libOCL.a -lstdc++
else
libs_ocl =
endif


libs_opt=$(libs_cuda)\
 $(libs_ocl)\
 $(libs_dgemmsy)



# Generic part, do not touch.
AM_FCFLAGS = -I$(top_builddir)/src -I$(top_builddir)/src/PSolver -I$(top_builddir)/src/modules @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@

INS = $(TESTDIRS:=.in)
RUNS = $(TESTDIRS:=.run)
CHECKS = $(TESTDIRS:=.check)
FAILEDCHECKS = $(TESTDIRS:=.recheck)
CLEANS = $(TESTDIRS:=.clean)

in: $(INS)

check: $(CHECKS) report

clean: $(CLEANS)

distclean : $(CLEANS)
	rm -rf Makefile

#distclean-check: clean

#	.PHONY: distclean-check

failed-check: $(FAILEDCHECKS) report

report:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

%.memguess.out: $(abs_top_builddir)/src/memguess
	name=`basename $@ .out.out | sed "s/[^_]*_\?\(.*\)$$/\1/"` ; \
	$(abs_top_builddir)/src/memguess 1 $$name > $@
%.out.out: $(abs_top_builddir)/src/bigdft
	name=`basename $@ .out.out | sed "s/[^_]*_\?\(.*\)$$/\1/"` ; \
	$(run_parallel) $(abs_top_builddir)/src/bigdft $$name > $@
%.geopt.mon.out:
	$(MAKE) -f ../Makefile $*.out.out && mv geopt.mon $@


%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref) mode="--mode=bigdft";; \
          *) mode="";; \
        esac ; \
        case $* in \
          CH3-.out*) prec="5.e-7" ;; \
          scalability.out*)  prec="1.e-9" ;; \
          *) prec="1.1e-10" ;; \
        esac ; \
	python $(abs_top_srcdir)/tests/fldiff.py $$mode --discrepancy=$$prec $*".out" $< | tee $@

%.clean:
	@dir=`basename $@ .clean` ; \
        if test x"$(srcdir)" = x"." ; then \
          rm -f $$dir.* $$dir/psppar.* $$dir/*.out $$dir/*.report $$dir/default* ; \
	  rm -fr $$dir/data ; \
        else \
          rm -rf $$dir.* $$dir ; \
        fi ; \
        echo "Clean in "$$dir

%.post-in: ;

in_message:
	@if test -n "$(run_ocl)" ; then \
	  echo "==============================================" ; \
	  echo "Will generate a 'input.perf' file to force OCL" ; \
	  echo "==============================================" ; \
	fi

$(INS): in_message
	@dir=`basename $@ .in` ; \
        if ! test x"$(srcdir)" = x"." ; then \
          if [ ! -d $$dir ] ; then mkdir $$dir ; fi ; \
          for i in $(srcdir)/$$dir/* ; do cp -f $$i $$dir; done ; \
        fi ; \
	if ! test -f $(srcdir)/$$dir/input.perf ; then \
	  rm -f $$dir/input.perf ; \
	  if test -n "$(run_ocl)" ; then \
	    echo "ACCEL OCLGPU" > $$dir/input.perf ; \
	  fi ; \
	fi ; \
        cd $$dir && $(MAKE) -f ../Makefile $$dir".post-in"; \
        echo "Input prepared in "$$dir" dir. make $$dir.run available"
	touch $@

run_message:
	@if test -n "$(run_parallel)" ; then \
	  echo "==============================================" ; \
	  echo "Will run tests in parallel with '$$run_parallel'" ; \
	  echo "==============================================" ; \
	fi

%.run: %.in run_message
	@dir=`basename $@ .run` ; \
        runs="$(srcdir)/$$dir/*.ref" ; \
	tgts=`for r in $$runs ; do echo $$(basename $$r .ref)".out"; done` ; \
        cd $$dir && $(MAKE) -f ../Makefile $$tgts ; \
        echo "Tests have run in "$$dir" dir. make $$dir.check available"
	touch $@

%.check: %.run
	@dir=`basename $@ .check` ; \
        chks="$(srcdir)/$$dir/*.ref" ; \
	tgts=`for c in $$chks ; do echo $$(basename $$c .ref)".report"; done` ; \
        cd $$dir && $(MAKE) -f ../Makefile $$tgts
	touch $@

%.recheck: %.in
	@dir=`basename $@ .recheck` ; \
        reps="$$dir/*.report" ; \
	for r in $$reps ; do \
	  if ! grep -qs "succeeded" $$r ; then \
	    target=` basename $$r .report` ; \
	    $(MAKE) $*".in" && rm -f $*".in" ; \
	    rm -f $$dir/$$target".out" $$dir/$$target".report" ; \
	    cd $$dir && $(MAKE) -f ../Makefile $$target".out" $$target".report" && cd - ; \
	  fi \
	done
	touch $*".check"

if USE_MPI
  mpirun_message=mpirun
else
  mpirun_message=
endif

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
