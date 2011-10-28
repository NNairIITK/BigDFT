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

# Additional NEB in.
NEB.post-in:
	sed "s;%%abs_top_builddir%%;$(abs_top_builddir);g" NEB_include.sh.in > NEB_include.sh

# Additional run rules.
Xabsdel_inputGenPot.xabs.post-out:
	mv -f local_potentialb2B.cube b2B_xanes.cube
	cp inputGenPot.xyz b2B_xanes.xyz

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

include $(srcdir)/../check.mk