## Process this file with automake to produce Makefile.in

SUBDIRS = 

EXTRA_DIST = \
	CPUlimit_global \
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
			 posloc* posmd* posout* poslow* posbest* possoft*\
			 global.out latest.pos.force.*.dat MINHOP.out
			 *.prc fort.* fldiff.report

#dependencies
$(PPFILES):
	ln -fs $(top_srcdir)/utils/PSPfiles/$@

MINHOP.out: $(top_builddir)/src/global $(PPFILES)
	if ! test x"$(srcdir)" = x"." ; then \
	  cp -f $(srcdir)/CPUlimit_global .; \
	  cp -f $(srcdir)/earr.dat .; \
	  cp -f $(srcdir)/input.dft .; \
	  cp -f $(srcdir)/mdinput.dft .; \
	  cp -f $(srcdir)/input.geopt .; \
	  cp -f $(srcdir)/mdinput.geopt .; \
	  cp -f $(srcdir)/ioput .; \
	  cp -f $(srcdir)/rand.inp .; \
	  cp -f $(srcdir)/poscur.xyz .; \
        fi
	$(run_parallel) $(top_builddir)/src/global > $@
	if ! test x"$(srcdir)" = x"." ; then \
	  rm -f CPUlimit_global ; \
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
	python $(top_srcdir)/tests/fldiff.py --discrepancy=5.e-9 --mode=bigdft MINHOP.out \
	$(top_srcdir)/tests/Refs/MINHOP.out.ref | tee fldiff.report
