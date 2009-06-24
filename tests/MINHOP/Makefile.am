## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	CPUlimit \
	README \
	earr.dat \
	input.dft \
	input.geopt \
	ioput \
	mdinput.dft \
	mdinput.geopt \
	poscur.xyz \
	rand.inp

PPFILES = \
	psppar.H \
	psppar.Si

CLEANFILES = $(PPFILES) inguess.new inguess.dat \
			 posloc* posmd* posout* poslow* posbest* \
			 global.out global.mon latest.pos.force.*.dat
			 *.prc fort.* fldiff.report

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/$@

MINHOP.out: $(top_builddir)/src/cluster $(PPFILES)
	if ! test x"$srcdir" = x"." ; then \
	  ln -fs $(srcdir)/CPUlimit ; \
	  ln -fs $(srcdir)/earr.dat ; \
	  ln -fs $(srcdir)/input.dft ; \
	  ln -fs $(srcdir)/mdinput.dft ; \
	  ln -fs $(srcdir)/input.geopt ; \
	  ln -fs $(srcdir)/mdinput.geopt ; \
	  ln -fs $(srcdir)/ioput ; \
	  ln -fs $(srcdir)/rand.inp ; \
	  ln -fs $(srcdir)/poscur.xyz ; \
        fi
	rm -f global.mon
	$(run_parallel) $(top_builddir)/src/global > $@
	if ! test x"$srcdir" = x"." ; then \
	  rm -f CPUlimit ; \
	  rm -f earr.dat ; \
	  rm -f input.dft ; \
	  rm -f mdinput.dft ; \
	  rm -f input.geopt ; \
	  rm -f mdinput.geopt ; \
	  rm -f ioput ; \
	  rm -f rand.inp ; \
	  rm -f poscur.xyz ; \
        fi

check: MINHOP.out
	python $(top_srcdir)/tests/fldiff.py global.mon $(top_srcdir)/tests/Refs/global.mon.ref
	python $(top_srcdir)/tests/fldiff.py --discrepancy=5.e-9 --bigdft MINHOP.out \
		   $(top_srcdir)/tests/Refs/MINHOP.out.ref | tee fldiff.report
