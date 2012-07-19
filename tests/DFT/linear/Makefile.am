# Give here the list of existing tests
TESTDIRS = Locreg \
	   base \
       directmin \
	   incdec \
	   withder \
	   mix \
	   water \
	   cleanup \
	   moreIGorbs

EXTRA_DIST = README

# Give here the pseudo-potentials used by each test.
base.psp: HGH/psppar.C extra/psppar.H
large.psp: HGH/psppar.C extra/psppar.H
directmin.psp : HGH/psppar.C extra/psppar.H
incdec.psp: HGH/psppar.C extra/psppar.H
standardnew.psp: HGH/psppar.C extra/psppar.H
withder.psp: HGH/psppar.C extra/psppar.H
mix.psp: HGH/psppar.C extra/psppar.H
water.psp: HGH/psppar.C extra/psppar.H
moreIGorbs.psp: HGH/psppar.C extra/psppar.H
Locreg.psp : HGH-K/psppar.H HGH-K/psppar.Si

# Define the precision for specific directories.
%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref) mode="--mode=bigdft";; \
          *) mode="";; \
        esac ; \
        case $* in \
	  base.out*)                  prec="1.e-9" ;; \
          large.out*)                 prec="5.e-10" ;; \
          directmin.out*)             prec="1.1e-10" ;; \
	  incdec_decreaselinear.out*) prec="2.e-10" ;; \
	  incdec_increaselocreg.out*) prec="1.1e-9" ;; \
	  standardnew.out*)           prec="2.e-10" ;; \
	  withder_base.out*)          prec="5.e-10" ;; \
	  mix_base.out*)              prec="2.e-10" ;; \
	  mix_new.out*)               prec="1.1e-9" ;; \
	  mix_potential.out*)         prec="1.e-9" ;; \
	  water.out*)                 prec="2.1e-8" ;; \
	  moreIGorbs.out*)            prec="3.e-10" ;; \
          *)                          prec="1.1e-10" ;; \
        esac ; \
	python $(abs_top_srcdir)/tests/fldiff.py $$mode --discrepancy=$$prec $*".out" $< | tee $@

%.report.yaml: %.ref.yaml
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	name=`basename $< .out.ref.yaml | sed "s/[^_]*_\?\(.*\)$$/\1/"`  ;\
	if test -n "$$name" ; then \
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d log-$$name.yaml -t $(srcdir)/../cubic/tols.yaml -o $@ ; \
	else \
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d log.yaml -t $(srcdir)/../cubic/tols.yaml -o $@ ; \
	fi 

include $(srcdir)/../../check.mk
