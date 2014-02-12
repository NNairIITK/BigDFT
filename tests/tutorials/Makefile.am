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


AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
else
LIBS += @LIBS_SHORT_DEPS@
endif


# Additional rules for the dev tutorial.
dev.check: dev.test/toy_model
dev.test/toy_model toy_model: $(abs_top_srcdir)/tests/tutorials/dev/toy_model.f90
	$(FC) -o $@ $(FCFLAGS) $(AM_FCFLAGS) $(abs_top_srcdir)/tests/tutorials/dev/toy_model.f90 \
			$(AM_LDFLAGS) $(LDFLAGS) -lbigdft-1 $(LIBS)
dev-N2.out.post-out: dev-N2.out.out toy_model 
	rm -f input.perf
	./toy_model > log-model.yaml
dev.post-clean:
	rm -f toy_model

# Define the precision for specific directories.
%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref ) mode="--mode=bigdft";; \
          *) mode="";; \
        esac ; \
        case $* in \
          *) prec="1.1e-10" ;; \
        esac ; \
	python $(abs_top_srcdir)/tests/fldiff.py $$mode --discrepancy=$$prec $*".out" $< | tee $@

%.report.yaml: %.ref.yaml
	@export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	name=`basename $< .out.ref.yaml | sed "s/[^_]*_\?\(.*\)$$/\1/"`  ;\
	tols=`basename $< .out.ref.yaml` ;\
	if test -n "$$name" ; then log=log-$$name.yaml ; else log=log.yaml ; fi ;\
	echo $$log; \
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@ 

include $(srcdir)/../check.mk
