## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXE_FILES = \
input.abscalc 	input.dft 	input.kpt 	input.kpt_20 \
posinp.xyz 	posinp.xyz_18 	posinp.xyz_72p 	psppar.O \
psppar.P 		psppar.Si posinp_72X27p.xyz  	posinp_72X8p.xyz


%: %.src
	cp $(srcdir)/$(@).src $(@)



exercise:  $(EXE_FILES)
