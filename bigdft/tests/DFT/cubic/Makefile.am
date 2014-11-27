# Give here the list of existing tests
TESTDIRS = C \
	  Li+ \
	  H2-spin \
	  SiH4 \
	  NaCl \
	  Ca2 \
	  H2O-CP2K \
	  SiBulk \
	  AgBulk \
	  TiMix \
	  TiSmear \
	  O2-Spin \
	  ElectField \
	  ElectronGas \
	  NLCC \
	  taskgroups \
	  CH4-D3 

tests_cpu_only = H Graphene FeHyb Mn2 SIC H-Rydberg
tests_ocl_only = OCL

EXTRA_TESTDIRS =

EXTRA_DIST = 

#This will run the automatic tests also for GPU
if USE_CUDA_GPU
  TESTDIRS += GPU
else
  EXTRA_TESTDIRS += GPU
endif
if USE_OCL
  TESTDIRS += $(tests_ocl_only)
  EXTRA_TESTDIRS += $(tests_cpu_only)
else
#Add H (box too small) and Graphene (different algorithm) only if we do not use OCL (see README)
#OCL not implemented with spin=4 (Fehyb and Mn2)
  TESTDIRS += $(tests_cpu_only)
  EXTRA_TESTDIRS += $(tests_ocl_only)
endif

# Give here the pseudo-potentials used by each test.
H.psp: psppar.H
Li+.psp: psppar.Li
H2-spin.psp: HGH/psppar.H
SiH4.psp: HGH-K/psppar.H HGH-K/psppar.Si
TiMix.psp: HGH-K/psppar.Ti HGH/psppar.Na
TiSmear.psp: HGH-K/psppar.Ti HGH/psppar.Na
NaCl.psp: HGH/psppar.Na HGH/psppar.Cl
Ca2.psp: psppar.Ca
Mn2.psp: psppar.Mn
H2O-CP2K.psp: HGH/psppar.H HGH/psppar.O
SiBulk.psp: HGH-K/psppar.Si
AgBulk.psp: HGH/psppar.Ag
Graphene.psp: psppar.C
FeHyb.psp: HGH/psppar.Fe
GPU.psp: HGH/psppar.Si
OCL.psp: HGH/psppar.Si
O2-Spin.psp: HGH-K/psppar.O
ElectField.psp: HGH-K/psppar.H HGH-K/psppar.O


SIC-test.post-in:
	if [ ! -d data-PZ ] ; then mkdir data-PZ; fi
	mv wavefunctions.gau data-PZ/

# Define the precision for specific directories.
%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref) mode="--mode=bigdft";; \
          *) mode="";; \
        esac ; \
        case $* in \
          *)                    prec="1.1e-10" ;; \
        esac ; \
    echo "python $(abs_top_srcdir)/tests/fldiff.py $$mode --discrepancy=$$prec $*.out $<" ; \
	python $(abs_top_srcdir)/tests/fldiff.py $$mode --discrepancy=$$prec $*".out" $< | tee $@

%.report.yaml: %.ref.yaml $(abs_top_builddir)/tests/fldiff_yaml.py
	@name=`basename $< .out.ref.yaml | sed "s/[^_]*_\?\(.*\)$$/\1/"`; \
	tols=`basename $< .out.ref.yaml` ;\
	if test -n "$$name" ; then log=log-$$name.yaml ; else log=log.yaml ; fi ;\
	echo "python $(abs_top_builddir)/tests/fldiff_yaml.py -r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@" ; \
	python $(abs_top_builddir)/tests/fldiff_yaml.py \
		-r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@ 

H-Rydberg-test.post-in:
	@if test -f "H.perf" ; then \
		echo >> H.yaml ; \
		echo "perf:" >> H.yaml ; \
		awk '{ printf("  %s: %s\n",$$1,$$2); }' H.perf >> H.yaml ; \
		echo "  multipole_preserving: Yes" >> H.yaml ; \
		rm H.perf; \
	fi

include $(srcdir)/../../check.mk


