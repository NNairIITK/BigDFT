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
       HGH/psppar.Mg \
       HGH/psppar.Ag \
       HGH-K/psppar.H \
       HGH-K/psppar.Si \
       HGH-K/psppar.N \
       HGH-K/psppar.O \
       HGH-K/psppar.Ti

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
	  O2-Spin \
	  GEOPT-BFGS \
	  GEOPT-LBFGS 

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
GEOPT-BFGS.psp: HGH/psppar.Mg
GEOPT-LBFGS.psp: HGH/psppar.Mg

# Additional rules.
H2O-slab.post-in:
	cp -f velocities.in velocities.xyz
	chmod a+w velocities.xyz

# Generic part, do not touch.
INS = $(TESTDIRS:=.in)
CHECKS = $(TESTDIRS:=.check)
FAILEDCHECKS = $(TESTDIRS:=.recheck)
CLEANS = $(TESTDIRS:=.clean)

in: $(INS)

check: $(CHECKS) report

clean: $(CLEANS)

distclean : $(CLEANS)
	rm -rf Makefile

failed-check: $(FAILEDCHECKS) report

report:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

%.memguess.out: $(abs_top_builddir)/src/memguess
	$(abs_top_builddir)/src/memguess 1 > $@
%.out.out: $(abs_top_builddir)/src/bigdft
	$(run_parallel) $(abs_top_builddir)/src/bigdft > $@
%.geopt.mon.out: $(abs_top_builddir)/src/bigdft
	$(MAKE) -f ../Makefile $*.out.out && mv geopt.mon $@

%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref) mode="--mode=bigdft";; \
          *) mode="";; \
        esac ; \
        case $* in \
          Ca2.out*)         prec="4.e-10" ;; \
          GEOPT-BFGS.out*)  prec="1.e-9" ;; \
          GEOPT-LBFGS.out*) prec="1.e-9" ;; \
          NaCl.out*)        prec="2.e-9" ;; \
          SiBulk.out*)      prec="2.e-7" ;; \
          AgBulk.out*)      prec="6.e-8" ;; \
          SiH4.out*)        prec="5.e-9" ;; \
          *)                prec="1.1e-10" ;; \
        esac ; \
	python $(abs_top_srcdir)/tests/fldiff.py $$mode --discrepancy=$$prec $*".out" $< | tee $@

$(PSPS):
	ln -fs $(abs_top_srcdir)/utils/PSPfiles/$@ 

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

%.in: in_message
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
        cd $$dir && $(MAKE) -f ../Makefile $$dir".psp" $$dir".post-in"; \
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
        refs="$$dir/*.ref" ; \
	for r in $$refs ; do \
	  rep=`basename $$r .ref`".report" ; \
	  if ! grep -qs "succeeded\|passed" $$dir/$$rep ; then \
	    target=` basename $$r .ref` ; \
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
