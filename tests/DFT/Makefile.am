PSPS = psppar.H \
       psppar.C \
       psppar.Li \
       psppar.Ca \
       psppar.Mn \
       psppar.Si \
       HGH/psppar.H \
       HGH/psppar.Na \
       HGH/psppar.Cl \
       HGH/psppar.O \
       HGH/psppar.Si \
       HGH/psppar.Fe \
       HGH-K/psppar.H \
       HGH-K/psppar.Si \
       HGH-K/psppar.N \
       HGH-K/psppar.O \
       HGH-K/psppar.Ti

# Give here the list of existing tests
TESTDIRS = H \
	  C \
	  Li+ \
	  H2-spin \
	  SiH4 \
	  NH3 \
	  NaCl \
	  Ca2 \
	  Mn2 \
	  H2O-CP2K \
	  SiBulk \
	  H2O-slab \
	  Davidson-SiH4 \
	  Graphene \
	  TiMix \
	  O2-Spin \
	  FeHyb 
#This will run the automatic tests also for GPU
if USE_CUDA_GPU
  TESTDIRS += GPU
endif
if USE_OCL
  TESTDIRS += OCL
endif

# Give here the pseudo-potentials used by each test.
H.psp: psppar.H
C.psp: psppar.C
Li+.psp: psppar.Li
H2-spin.psp: HGH/psppar.H
SiH4.psp: HGH-K/psppar.H HGH-K/psppar.Si
NH3.psp: HGH-K/psppar.N HGH-K/psppar.H
TiMix.psp: HGH-K/psppar.Ti
NaCl.psp: HGH/psppar.Na HGH/psppar.Cl
Ca2.psp: psppar.Ca
Mn2.psp: psppar.Mn
H2O-CP2K.psp: HGH/psppar.H HGH/psppar.O
SiBulk.psp: HGH/psppar.Si
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





# Generic part, do not touch.
INS = $(TESTDIRS:=.in)
CHECKS = $(TESTDIRS:=.check)
FAILEDCHECKS = $(TESTDIRS:=.recheck)
CLEANS = $(TESTDIRS:=.clean)

in: $(INS)

check: $(CHECKS) report

clean: $(CLEANS)

failed-check: $(FAILEDCHECKS) report

report:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

%.memguess.out: $(abs_top_builddir)/src/memguess
	$(abs_top_builddir)/src/memguess 1 > $@
%.out.out: $(abs_top_builddir)/src/bigdft
	$(run_parallel) $(abs_top_builddir)/src/bigdft > $@
%.geopt.mon.out:
	$(MAKE) -f ../Makefile $*.out.out && mv geopt.mon $@

%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref) mode="--mode=bigdft";; \
          *) mode="";; \
        esac ; \
        case $* in \
          Ca2.out*) prec="4.e-10" ;; \
          Davidson-SiH4.out*) prec="4.e-9" ;; \
          NaCl.out*) prec="2.e-8" ;; \
          SiBulk.out*) prec="5.e-10" ;; \
          SiH4.out*) prec="1.3e-9" ;; \
          *) prec="1.1e-10" ;; \
        esac ; \
	python $(abs_top_srcdir)/tests/fldiff.py $$mode --discrepancy=$$prec $*".out" $< | tee $@

$(PSPS):
	ln -fs $(abs_top_srcdir)/utils/PSPfiles/$@ 

%.clean:
	@dir=`basename $@ .clean` ; \
        if test x"$(srcdir)" = x"." ; then \
          rm -f $$dir.* $$dir/psppar.* $$dir/*.out $$dir/*.report ; \
        else \
          rm -rf $$dir.* $$dir ; \
        fi ; \
        echo "Clean in "$$dir

%.post-in: ;

%.in:
	@dir=`basename $@ .in` ; \
        if ! test x"$(srcdir)" = x"." ; then \
          if [ ! -d $$dir ] ; then mkdir $$dir ; fi ; \
          for i in $(srcdir)/$$dir/* ; do cp -f $$i $$dir; done ; \
        fi ; \
        cd $$dir && $(MAKE) -f ../Makefile $$dir".psp" $$dir".post-in"; \
        echo "Input prepared in "$$dir" dir. make $$dir.run available"
	touch $@

%.run: %.in
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

EXTRA_DIST = \
	README \
	inguess.dat \
        $(TESTDIRS) \
        GPU OCL
