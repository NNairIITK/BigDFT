# Give here the list of existing tests

#SUBDIRS = Xabs/pawpatchgen  Xabs/scfpotgen   Xabs/spectra    

TESTDIRS = CH4 \
	  CH3- \
	  H2O-slab \
	  H2O-bader \
	  scalability \
	  dev

if USE_OCL
  TESTDIRS +=
else
  TESTDIRS += N2
endif

EXTRA_DIST = Figures \
     H1A_first_run.html \
     H1B_grid_cv.html \
     H1C_basic_surfaces.html \
     H1D_wfn_plot.html \
     H2A_MPI_OMP.html \
     H2B_GPU.html \
     H2C-BART.html \
     H2D-splined_saddle.html \
     H2E-minima_hopping.html \
     H2E1-troubleshooting_minima_hopping.html \
     H3B-charge_analysis.html \
     H3C-Wannier.html \
     H3D-Xabs.html \
     H3X-hacking_bigdft.html \
     Mg7-minhop \
     minhop.pdf \
     bart.pdf \
     GPU \
     C2H6 \
     Mg5 \
     Li6 \
     Xabs \
     Wannier

# Additional post-in commands
# Remove HF.perf in case of acceleration
N2.post-in :
	rm -f input.perf

if BUILD_DYNAMIC_LIBS
bigdft_library = -L$(top_builddir)/src -lbigdft-1
bigdft_library_rpath = -Wl,-rpath=$(abs_top_builddir)/src
else
bigdft_library = $(top_builddir)/src/libbigdft.a @LIBS_SHORT_DEPS@
bigdft_library_rpath =
endif

# Additional rules for the dev tutorial.
dev.check: dev/toy_model
dev/toy_model toy_model: $(abs_top_srcdir)/tests/tutorials/dev/toy_model.f90
	$(FC) -o $@ $(FCFLAGS) $(AM_FCFLAGS) $(abs_top_srcdir)/tests/tutorials/dev/toy_model.f90 \
			$(LDFLAGS) $(bigdft_library_rpath) $(bigdft_library) $(LIBS)
toy_model.dat.out: dev-N2.out.out toy_model 
	./toy_model > toy_model.dat.out
dev.post-clean:
	rm -f toy_model

# Define the precision for specific directories.
%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref) mode="--mode=bigdft";; \
          *) mode="";; \
        esac ; \
        case $* in \
          CH3-.out*) prec="5.2e-7" ;; \
          scalability.out*)  prec="1.e-9" ;; \
          *) prec="1.1e-10" ;; \
        esac ; \
	python $(abs_top_srcdir)/tests/fldiff.py $$mode --discrepancy=$$prec $*".out" $< | tee $@

%.report.yaml: %.ref.yaml
	@export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	name=`basename $< .out.ref.yaml | sed "s/[^_]*_\?\(.*\)$$/\1/"`  ;\
	tols=`basename $< .out.ref.yaml` ;\
	if test -n "$$name" ; then log=log-$$name.yaml ; else log=log.yaml ; fi ;\
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@ 

include $(srcdir)/../check.mk
