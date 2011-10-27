# Give here the list of existing tests
TESTDIRS = \
	   GEOPT-BFGS \
	   GEOPT-LBFGS \
           H2-freq \
	   NEB \
	   SPLSAD \
	   XabsCheb \
	   XabsLanc \
	   Xabsb2B

#This will run the automatic tests also for minima hopping
if USE_MINIMA_HOPPING
TESTDIRS += MINHOP
endif

EXTRA_DIST = \
        $(TESTDIRS) \
	MINHOP

# Give here the pseudo-potentials used by each test.
H2-freq.psp: psppar.H
NEB.psp: psppar.H
SPLSAD.psp: HGH/psppar.N
MINHOP.psp: psppar.Si psppar.H
GEOPT-BFGS.psp: HGH/psppar.Mg
GEOPT-LBFGS.psp: HGH/psppar.Mg
XabsCheb.psp: Xabs/psppar.Fe
XabsLanc.psp: Xabs/psppar.Fe
Xabsb2B.psp: Xabs/psppar.Fe

# Additional dependencies
H2-freq2.freq.out: H2-freq.freq.out
Xabsb2B.xabs.out: Xabsdel_inputGenPot.xabs.out

# Additional run rules.
Xabsdel_inputGenPot.xabs.out: $(abs_top_builddir)/src/abscalc
	name=`basename $@ .xabs.out | sed "s/[^_]*_\?\(.*\)$$/\1/"` ; \
	$(abs_top_builddir)/src/abscalc $$name > $@
	mv -f local_potentialb2B.cube b2B_xanes.cube
	cp inputGenPot.xyz b2B_xanes.xyz

# Additional clean targets.
SPLSAD.post-clean:
	rm -f anchorpoints* datatime.prc fort.* nogt.* path*.xyz vogt.*
MINHOP.post-clean:
	rm -f *pos* *out *mon global.out latest.pos.force.*.dat MINHOP.out \
              *.prc fort.* default* CPUlimit test
XabsCheb.post-clean:
	rm -f cheb_spectra_* electronic_density.cube
XabsLanc.post-clean:
	rm -f alphabeta* electronic_density.cube
Xabsb2B.post-clean:
	rm -f cheb_spectra_* b2B_xanes.* electronic_density.cube local_potentialb2B*

# Define the precision for specific directories.
%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref | *.freq.ref | *.splsad.ref | *.minhop.ref | *.xabs.ref ) mode="--mode=bigdft";; \
          *.NEB.ref) mode="--mode=neb";; \
          *) mode="";; \
        esac ; \
        case $* in \
          GEOPT-BFGS.out*)  prec="1.e-9" ;; \
          GEOPT-LBFGS.out*) prec="1.e-9" ;; \
          H2-freq*)         prec="1.e-8" ;; \
	  NEB*)             prec="5e-5" ;; \
	  SPLSAD*)          prec="2e-8" ;; \
	  MINHOP*)          prec="5e-9" ;; \
          *)                prec="1.1e-10" ;; \
        esac ; \
	python $(abs_top_srcdir)/tests/fldiff.py $$mode --discrepancy=$$prec $*".out" $< | tee $@

# Specific clean traget for NEB.
NEB.clean:
	@dir=`basename $@ .clean` ; \
	rm -f $$dir.* $$dir/triH.NEB.* $$dir/NEB.NEB.out $$dir/psppar.H $$dir/*.report $$dir/NEB_driver.sh $$dir/malloc.prc ; \
	rm -fr $$dir/data ; \
        echo "Clean in "$$dir

include $(srcdir)/../check.mk