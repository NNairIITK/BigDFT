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
H2O-bader.post-in :
	cp -f input-all.dft input.dft
# Remove HF.perf in case of acceleration
N2.post-in :
	rm -f input.perf


# Additional rules for the dev tutorial.
dev.check: dev/toy_model
dev/toy_model toy_model: $(abs_top_srcdir)/tests/tutorials/dev/toy_model.f90
	$(FC) -o $@ $(FCFLAGS) $(AM_FCFLAGS) $(abs_top_srcdir)/tests/tutorials/dev/toy_model.f90 $(LDFLAGS) @LIBS_DEPENDENCIES@ $(LIBS)
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
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	name=`basename $< .out.ref.yaml | sed "s/[^_]*_\?\(.*\)$$/\1/"`  ;\
        case $* in \
          scalability*)     tols="tols-scalability.yaml" ;; \
          *)                tols="../DFT/cubic/tols.yaml" ;; \
        esac ;\
	if test -n "$$name" ; then \
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d log-$$name.yaml -t $(srcdir)/$$tols -o $@ ; \
	else \
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d log.yaml -t $(srcdir)/$$tols -o $@ ; \
	fi 


include $(srcdir)/../check.mk
