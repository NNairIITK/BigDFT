# Give here the list of existing tests
TESTDIRS = \
	   C6H6-freq \
	   GEOPT-all \
	   GEOPT-BFGS \
	   GEOPT-LBFGS \
	   GEOPT-newIG \
	   MINHOP \
	   NH3 \
	   H2O-slab \
	   H2-freq \
	   NEB \
	   SPLSAD \
	   XabsCheb \
	   XabsLanc \
	   Xabsb2B \
	   BigDFT2Wannier \
       testforces \
	   IO-wf

EXTRA_DIST = MH_taskgroups

if HAVE_ETSF_IO
TESTDIRS += IO-wf-ETSF
endif


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
NH3.psp: HGH-K/psppar.N HGH-K/psppar.H
H2O-slab.psp: HGH-K/psppar.H HGH-K/psppar.O

# Additional dependencies
H2-freq2.freq.out: H2-freq.freq.out 
Xabsb2B.xabs.out: InputGenPot.xabs.out
IO-wf_plain-read.out.out: IO-wf_plain.out.out
IO-wf_binary-read.out.out: IO-wf_binary.out.out
IO-wf-ETSF_read.out.out: IO-wf-ETSF_write.out.out


# Additional freq post-in to avoid a fail when the run is done again.
H2-freq-test.post-in:
	rm -f frequencies.res
C6H6-freq-test.post-in: 
	rm -f frequencies.res

# Additional MINHOP in;
MINHOP-test.post-in:
	rm -f posloc* posimed* data/global.mon

# Additional NEB in.
NEB-test.post-in:
	if ! test -f NEB_include.sh ; then sed "s;%%abs_top_builddir%%;$(abs_top_builddir);g" NEB_include.sh.in > NEB_include.sh ; fi
	if ! test -L NEB_driver.sh ; then ln -fs $(abs_top_srcdir)/src/NEB_driver.sh ; fi
	rm -f neb.NEB.restart
NEB-test.post-clean:
	rm -f NEB_include.sh NEB_driver.sh

# Additional H2O-slab in.
H2O-slab-test.post-in:
	cp -f velocities.in velocities.xyz
	chmod a+w velocities.xyz

# Additional reference for geopt.mon test
global.mon.out : test.minhop.out
	cp -f data/global.mon $@


# Additional outs
H2-freq.post-out:
	cp -f log.yaml log-freq.yaml

H2-freq2.post-out:
	cat log.yaml >> log-freq.yaml
	mv log-freq.yaml log.yaml

IO-wf_plain.out.post-out:
	if test -e data-plain-read; then rm -rf data-plain-read; fi
	mv -f data-plain data-plain-read
	$(MAKE) -f ../Makefile IO-wf_plain.cube

IO-wf_binary.out.post-out:
	if test -e data-binary-read; then rm -rf data-binary-read; fi
	mv -f data-binary data-binary-read
	$(MAKE) -f ../Makefile IO-wf_binary.cube

IO-wf-ETSF_write.out.post-out:
	mv -f data-write data-read
	$(MAKE) -f ../Makefile IO-wf_ETSF.cube

NH3.out.post-out:
	mv -f data/geopt.mon NH3.geopt.mon.out

InputGenPot.xabs.post-out:
	mv -f local_potentialb2B.cube b2B_xanes.cube
	cp InputGenPot.xyz b2B_xanes.xyz


# Introduce specific targets for export-wavefunctions
IO-wf_plain.cube: 
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
    rm plain-read.perf ;\
	$(run_serial) $(abs_top_builddir)/src/bigdft-tool --name plain-read -a export-wf data-plain-read/wavefunction-k002-DI.b0003 >> log-plain.yaml
	cp check.perf plain-read.perf
	mv wavefunction-k002-DI.b0003.cube wf_plain.cube
	rm -f wavefunction-k002-DI.b0003_avg_*

IO-wf_binary.cube: 
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
    rm binary-read.perf ;\
	$(run_serial) $(abs_top_builddir)/src/bigdft-tool --name binary-read -a export-wf data-binary-read/wavefunction-k002-DI.bin.b0003 >> log-binary.yaml
	cp check.perf binary-read.perf
	mv wavefunction-k002-DI.b0003.cube wf_binary.cube
	rm -f wavefunction-k002-DI.b0003_avg_*

IO-wf_ETSF.cube: 
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
	$(run_serial) $(abs_top_builddir)/src/bigdft-tool --name read -a export-wf data-read/wavefunction.etsf \
	--i-band 3 --i-spin 2 --i-kpt 2 --i-complex 2 >> IO-wf-ETSF_write.out.out
	mv wavefunction-k002-DI.b0003.cube wf_etsf.cube
	rm -f wavefunction-k002-DI.b0003_avg_*

# Define the precision for specific directories.
%.report: %.ref
	@case $< in \
          *.memguess.ref | *.out.ref | *.freq.ref | *.splsad.ref | *.minhop.ref | *.xabs.ref | *.b2w.ref ) mode="--mode=bigdft";; \
          *.NEB.ref) mode="--mode=neb";; \
          *) mode="";; \
        esac ; \
        case $* in \
          GEOPT-BFGS.out*)  prec="1.e-9" ;; \
          GEOPT-LBFGS.out*) prec="1.e-9" ;; \
	  NEB*)             prec="5e-5" ;; \
	  SPLSAD*)          prec="2e-8" ;; \
	  MINHOP*)          prec="5e-9" ;; \
          *)                prec="1.1e-10" ;; \
        esac ; \
	python $(abs_top_srcdir)/tests/fldiff.py $$mode --discrepancy=$$prec $*".out" $< | tee $@

%.report.yaml: %.ref.yaml
	@export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	name2=`basename $< .xabs.ref.yaml` ; \
	if test "$$name2" = "$<" ; then \
	name=`basename $< .out.ref.yaml | sed "s/[^_]*_\?\(.*\)$$/\1/"`  ;\
	else \
	name=`basename $< .xabs.ref.yaml`  ; fi;\
	tols=`basename $< .out.ref.yaml` ;\
	if test -n "$$name" ; then log=log-$$name.yaml ; else log=log.yaml ; fi ; \
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@ 

include $(srcdir)/../check.mk
