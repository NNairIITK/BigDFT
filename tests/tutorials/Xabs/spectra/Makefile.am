## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXE_FILES = \
input.abscalc 	input.dft 	input.kpt 	input.kpt_20 \
posinp.xyz 	posinp.xyz_18 	posinp.xyz_72p 	psppar.O \
psppar.P 		psppar.Si



# CLEANFILES = $(PPFILES) Xabs*.out posout* *.pot *.prc *.cube *.report \
# 	coeffs_chebychev cheb_spectra_199 cheb_spectra_399 \
# 	fldiff.report  cheb_spectra_* alpha_beta* \
# 	alphabeta* local_potential* potx.dat  input* posi*
#dependencies

$(EXE_FILES):
	cp  $(srcdir)/input.abscalc      .
	cp  $(srcdir)/input.dft     .
	cp  $(srcdir)/input.kpt     .
	cp  $(srcdir)/posinp.xyz     .
	cp  $(srcdir)/psppar.O     .
	cp  $(srcdir)/psppar.Si     .


exercise:  $(EXE_FILES)
