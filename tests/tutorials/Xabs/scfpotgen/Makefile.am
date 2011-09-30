## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXE_FILES = input.dft posinp.xyz  psppar.O psppar.Si  psppar.P



# CLEANFILES = $(PPFILES) Xabs*.out posout* *.pot *.prc *.cube *.report \
# 	coeffs_chebychev cheb_spectra_199 cheb_spectra_399 \
# 	fldiff.report  cheb_spectra_* alpha_beta* \
# 	alphabeta* local_potential* potx.dat  input* posi*
#dependencies

$(EXE_FILES): 
	cp  $(srcdir)/input.dft   .
	cp  $(srcdir)/posinp.xyz   .
	cp  $(srcdir)/psppar.O   .
	cp  $(srcdir)/psppar.Si   .
	cp  $(srcdir)/psppar.P   .


exercise:  $(EXE_FILES)