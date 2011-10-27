# Give here the list of existing tests
TESTDIRS = C \
	  Li+ \
	  H2-spin \
	  SiH4 \
	  NH3 \
	  NaCl \
	  Ca2 \
	  H2O-CP2K \
	  SiBulk \
	  AgBulk \
	  H2O-slab \
	  Davidson-SiH4 \
	  TiMix \
	  O2-Spin

#This will run the automatic tests also for GPU
if USE_CUDA_GPU
  TESTDIRS += GPU
endif
if USE_OCL
  TESTDIRS += OCL
else
#Add H (box too small) and Graphene (different algorithm) only if we do not use OCL (see README)
#OCL not implemented with spin=4 (Fehyb and Mn2)
  TESTDIRS += H Graphene FeHyb Mn2
endif

EXTRA_DIST = \
	README \
	inguess.dat \
        $(TESTDIRS) \
        GPU OCL

# Give here the pseudo-potentials used by each test.
H.psp: psppar.H
C.psp: psppar.C
Li+.psp: psppar.Li
H2-spin.psp: HGH/psppar.H
SiH4.psp: HGH-K/psppar.H HGH-K/psppar.Si
NH3.psp: HGH-K/psppar.N HGH-K/psppar.H
TiMix.psp: HGH-K/psppar.Ti HGH/psppar.Na
NaCl.psp: HGH/psppar.Na HGH/psppar.Cl
Ca2.psp: psppar.Ca
Mn2.psp: psppar.Mn
H2O-CP2K.psp: HGH/psppar.H HGH/psppar.O
SiBulk.psp: HGH/psppar.Si
AgBulk.psp: HGH/psppar.Ag
H2O-slab.psp: HGH-K/psppar.H HGH-K/psppar.O
Davidson-SiH4.psp: psppar.H psppar.Si
Graphene.psp: psppar.C
FeHyb.psp: HGH/psppar.Fe
GPU.psp: HGH/psppar.Si
OCL.psp: HGH/psppar.Si
O2-Spin.psp: HGH-K/psppar.O

# Additional rules.
H2O-slab.post-in:
	cp -f velocities.in velocities.xyz
	chmod a+w velocities.xyz

# Define the precision for specific directories.
%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref) mode="--mode=bigdft";; \
          *) mode="";; \
        esac ; \
        case $* in \
          Ca2.out*)         prec="4.e-10" ;; \
          NaCl.out*)        prec="2.e-9" ;; \
          SiBulk.out*)      prec="2.e-7" ;; \
          AgBulk.out*)      prec="6.e-8" ;; \
          SiH4.out*)        prec="5.e-9" ;; \
          *)                prec="1.1e-10" ;; \
        esac ; \
	python $(abs_top_srcdir)/tests/fldiff.py $$mode --discrepancy=$$prec $*".out" $< | tee $@

include $(srcdir)/../check.mk