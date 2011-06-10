## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	inputRef.dft \
	inputRef.abscalc \
	posinpRef.xyz \
	XabsCheb.out.ref \
	psppar.Fe_patched_plus_initialwave	

PPFILES = \
	psppar.Fe

CLEANFILES = $(PPFILES) Xabs*.out posout* *.pot *.prc \
	coeffs_chebychev cheb_spectra_199 cheb_spectra_399 \
	fldiff.report  del.del cheb_spectra_* alpha_beta* \
	alphabeta* local_potential* potx.dat

#dependencies
$(PPFILES):
	@ln -fs $(srcdir)/psppar.Fe_patched_plus_initialwave   psppar.Fe

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
	$(top_builddir)/src/abscalc > del.del
	mv -f local_potentialb2B.cube  b2B_xanes.cube
	cp posinp.xyz  b2B_xanes.xyz
	@rm -f input.dft
	@rm -f posinp.xyz
	@rm -f input.abscalc
	@ln -fs $(srcdir)/inputRef.dft input.dft
	@ln -fs $(srcdir)/posinpRef.xyz posinp.xyz
	@ln -fs $(srcdir)/inputb2B.abscalc input.abscalc
	$(top_builddir)/src/abscalc > $@
	rm -f  b2B_xanes.xyz
	rm -f  b2B_xanes.cube
	@rm -f input.dft
	@rm -f posinp.xyz
	@rm -f input.abscalc


XabsCG.out: $(top_builddir)/src/abscalc $(PPFILES)
	@ln -fs $(srcdir)/inputCG.abscalc input.abscalc
	@ln -fs $(srcdir)/inputRef.dft input.dft
	@ln -fs $(srcdir)/posinpRef.xyz posinp.xyz
	$(top_builddir)/src/abscalc > $@
	@rm -f input.dft
	@rm -f posinp.xyz
	@rm -f input.abscalc




check: XabsCheb.out  XabsLanc.out Xabsb2B.out  XabsCG.out
	cat XabsCheb.out  XabsLanc.out  Xabsb2B.out XabsCG.out >  XabsChebLanc.out 
	python $(top_srcdir)/tests/fldiff.py --mode=bigdft XabsChebLanc.out \
		$(top_srcdir)/tests/XabsCheb/XabsCheb.out.ref | tee fldiff.report
