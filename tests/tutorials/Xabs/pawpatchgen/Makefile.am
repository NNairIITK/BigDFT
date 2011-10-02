## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXE_FILES = atom.dat 	input.dat 	psppar psppar.P  


%: %.src
	cp $(srcdir)/$(@).src $(@)


exercise: $(EXE_FILES)    $(top_srcdir)/pseudo/src/pawpatchfit.py
	cp $(top_srcdir)/pseudo/src/pawpatchfit.py .