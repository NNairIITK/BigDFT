## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	CPUlimit \
	README \
	earr.dat \
	input.dat \
	ioput \
	mdinput.dat \
	poscur.xyz \
	rand.inp

PPFILES = \
	psppar.H \
	psppar.Si

CLEANFILES = $(PPFILES) inguess.new inguess.dat \
			 posloc* posmd* posout* poslow* posbest* \
			 global.out global.mon latest.pos.force.*.dat
			 *.prc fldiff.report $(EXTRA_DIST)

#dependencies
$(PPFILES):
	cp $(top_srcdir)/utils/PSPfiles/$@ .

check: $(PPFILES) 
	cp -f $(top_srcdir)/tests/MINHOP/CPUlimit .
	cp -f $(top_srcdir)/tests/MINHOP/earr.dat .
	cp -f $(top_srcdir)/tests/MINHOP/input.dat .
	cp -f $(top_srcdir)/tests/MINHOP/mdinput.dat .
	cp -f $(top_srcdir)/tests/MINHOP/ioput .
	cp -f $(top_srcdir)/tests/MINHOP/rand.inp .
	cp -f $(top_srcdir)/tests/MINHOP/poscur.xyz .
	rm -f global.mon
	$(run_parallel) $(top_builddir)/src/global > MINHOP.out
	python $(top_srcdir)/tests/fldiff.py global.mon $(top_srcdir)/tests/Refs/global.mon.ref
	python $(top_srcdir)/tests/fldiff.py --bigdft MINHOP.out \
		   $(top_srcdir)/tests/Refs/MINHOP.out.ref | tee fldiff.report
