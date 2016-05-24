## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	inputRef.dft \
	inputRef.abscalc \
	posinpRef.xyz \
	XabsCheb.out.ref \
	psppar.Fe_patched_plus_initialwave	\
	inputKpt.kpt \
	inputRef.abscalc \
	input.abscalcLanc \
	input.abscalcLanc \
	inputGenPot.dft \
	posinpGenPot.xyz \
	inputb2B.abscalc \
	inputCG.abscalc \
	XabsCheb.out.ref \
	XabsLanc.out.ref  \
	Xabsb2B.out.ref \
	XabsCG.out.ref

PPFILES = psppar.Fe

CLEANFILES = $(PPFILES) Xabs*.out posout* *.pot *.prc *.cube *.report \
	coeffs_chebychev cheb_spectra_199 cheb_spectra_399 \
	fldiff.report  cheb_spectra_* alpha_beta* \
	alphabeta* local_potential* potx.dat  input* posi* default*

clean.local : 
	rm -rf data

distclean : clean clean.local
	rm -rf Makefile

#dependencies
$(PPFILES):
	@cp $(srcdir)/psppar.Fe_patched_plus_initialwave psppar.Fe

XabsCheb.out: $(top_builddir)/src/abscalc $(PPFILES)
	@ln -fs $(srcdir)/inputRef.abscalc input.abscalc
	@ln -fs $(srcdir)/inputRef.dft input.dft
	@ln -fs $(srcdir)/posinpRef.xyz posinp.xyz
	$(top_builddir)/src/abscalc > $@
	@rm -f input.dft
	@rm -f posinp.xyz
	@rm -f input.abscalc

XabsLanc.out: $(top_builddir)/src/abscalc $(PPFILES)
	@ln -fs $(srcdir)/inputRef.dft input.dft
	@ln -fs $(srcdir)/posinpRef.xyz posinp.xyz
	@ln -fs $(srcdir)/input.abscalcLanc input.abscalc
	$(top_builddir)/src/abscalc > $@
	@rm -f input.dft
	@rm -f posinp.xyz
	@rm -f input.abscalc


Xabsb2B.out: $(top_builddir)/src/abscalc $(PPFILES)
	@ln -fs $(srcdir)/inputGenPot.dft input.dft
	@ln -fs $(srcdir)/posinpGenPot.xyz posinp.xyz
	@ln -fs $(srcdir)/inputRef.abscalc input.abscalc
	$(top_builddir)/src/abscalc > Xabsdel.out
	mv -f local_potentialb2B.cube  b2B_xanes.cube
	cp posinp.xyz  b2B_xanes.xyz
	@rm -f input.dft
	@rm -f posinp.xyz
	@rm -f input.abscalc
	@ln -fs $(srcdir)/inputRef.dft input.dft
	@ln -fs $(srcdir)/posinpRef.xyz posinp.xyz
	@ln -fs $(srcdir)/inputb2B.abscalc input.abscalc	
	@ln -fs $(srcdir)/inputKpt.kpt input.kpt
	$(top_builddir)/src/abscalc > $@
	rm -f  b2B_xanes.xyz
	rm -f  b2B_xanes.cube
	@rm -f input.dft
	@rm -f posinp.xyz
	@rm -f input.abscalc
	@rm -f input.kpt

XabsCG.out: $(top_builddir)/src/abscalc $(PPFILES)
	@ln -fs $(srcdir)/inputCG.abscalc input.abscalc
	@ln -fs $(srcdir)/inputRef.dft input.dft
	@ln -fs $(srcdir)/posinpRef.xyz posinp.xyz
	$(top_builddir)/src/abscalc > $@
	@rm -f input.dft
	@rm -f posinp.xyz
	@rm -f input.abscalc


#dependencies for make check in parallel
#Xabsb2B.out: XabsLanc.out
#XabsLanc.out: XabsCheb.out


check: XabsCheb.out  XabsLanc.out Xabsb2B.out
	python $(pythondir)/fldiff.py --mode=bigdft XabsCheb.out \
		$(top_srcdir)/tests/XabsCheb/XabsCheb.out.ref | tee XabsCheb.out.report
	python $(pythondir)/fldiff.py --mode=bigdft XabsLanc.out \
		$(top_srcdir)/tests/XabsCheb/XabsLanc.out.ref | tee XabsLanc.out.report
	python $(pythondir)/fldiff.py --mode=bigdft Xabsb2B.out \
		$(top_srcdir)/tests/XabsCheb/Xabsb2B.out.ref | tee Xabsb2B.out.report


checkCG: XabsCG.out
	python $(pythondir)/fldiff.py --mode=bigdft XabsCG.out \
		$(top_srcdir)/tests/XabsCheb/XabsCG.out.ref | tee XabsCG.out.report



#	python $(pythondir)/fldiff.py --mode=bigdft Xabsdel.out \
#		$(top_srcdir)/tests/XabsCheb/Xabsdel.out.ref | tee Xabsdel.out.report
