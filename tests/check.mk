# Generic part, of the testing Makefiles.
# Possible calls are:
#  make in: generate all input dirs.
#  make failed-check: run check again on all directories with missing report
#                     or failed report.
#  make X.in: generate input dir for directory X.
#  make X.check: generate a report for directory X (if not already existing).
#  make X.recheck: force the creation of the report in directory X.
#  make X.clean: clean the given directroy X.

AM_FCFLAGS = -I$(top_builddir)/src -I$(top_builddir)/src/PSolver -I$(top_builddir)/src/modules @LIBABINIT_INCLUDE@ @LIBXC_INCLUDE@

INS = $(TESTDIRS:=.in)
RUNS = $(TESTDIRS:=.run)
CHECKS = $(TESTDIRS:=.check)
FAILEDCHECKS = $(TESTDIRS:=.recheck)
CLEANS = $(TESTDIRS:=.clean)

in: $(INS)

check: $(CHECKS) report

clean: $(CLEANS)

distclean : $(CLEANS)
	rm -rf Makefile

failed-check: $(FAILEDCHECKS) report

report:
	@if test $(MAKELEVEL) = 0 ; then python $(top_srcdir)/tests/report.py ; fi

%.memguess.out: $(abs_top_builddir)/src/memguess
	name=`basename $@ .out.out | sed "s/[^_]*_\?\(.*\)$$/\1/"` ; \
	$(abs_top_builddir)/src/memguess 1 $$name > $@
%.out.out: $(abs_top_builddir)/src/bigdft
	name=`basename $@ .out.out | sed "s/[^_]*_\?\(.*\)$$/\1/"` ; \
	$(run_parallel) $(abs_top_builddir)/src/bigdft $$name > $@
%.geopt.mon.out:
	$(MAKE) -f ../Makefile $*.out.out && mv geopt.mon $@
%.dipole.dat.out: %.out.out
	$(abs_top_builddir)/src/tools/bader/bader data/electronic_density.cube > bader.out && mv dipole.dat $@

$(PSPS):
	ln -fs $(abs_top_srcdir)/utils/PSPfiles/$@ 

%.clean:
	@dir=`basename $@ .clean` ; \
        if test x"$(srcdir)" = x"." ; then \
          rm -f $$dir.* $$dir/psppar.* $$dir/*.out $$dir/*.report $$dir/default* ; \
	  rm -fr $$dir/data ; \
        else \
          rm -rf $$dir.* $$dir ; \
        fi ; \
        echo "Clean in "$$dir

%.post-in: ;
%.psp: ;

in_message:
	@if test -n "$(run_ocl)" ; then \
	  echo "==============================================" ; \
	  echo "Will generate a 'input.perf' file to force OCL" ; \
	  echo "==============================================" ; \
	fi

$(INS): in_message
	@dir=`basename $@ .in` ; \
        if ! test x"$(srcdir)" = x"." ; then \
          if [ ! -d $$dir ] ; then mkdir $$dir ; fi ; \
          for i in $(srcdir)/$$dir/* ; do cp -f $$i $$dir; done ; \
        fi ; \
	if ! test -f $(srcdir)/$$dir/input.perf ; then \
	  rm -f $$dir/input.perf ; \
	  if test -n "$(run_ocl)" ; then \
	    echo "ACCEL OCLGPU" > $$dir/input.perf ; \
	  fi ; \
	fi ; \
        cd $$dir && $(MAKE) -f ../Makefile $$dir".psp"; \
        $(MAKE) -f ../Makefile $$dir".post-in"; \
        echo "Input prepared in "$$dir" dir. make $$dir.run available"
	touch $@

run_message:
	@if test -n "$(run_parallel)" ; then \
	  echo "==============================================" ; \
	  echo "Will run tests in parallel with '$$run_parallel'" ; \
	  echo "==============================================" ; \
	fi

%.run: %.in run_message
	@dir=`basename $@ .run` ; \
        runs="$(srcdir)/$$dir/*.ref" ; \
	tgts=`for r in $$runs ; do echo $$(basename $$r .ref)".out"; done` ; \
        cd $$dir && $(MAKE) -f ../Makefile $$tgts ; \
        echo "Tests have run in "$$dir" dir. make $$dir.check available"
	touch $@

%.check: %.run
	@dir=`basename $@ .check` ; \
        chks="$(srcdir)/$$dir/*.ref" ; \
	tgts=`for c in $$chks ; do echo $$(basename $$c .ref)".report"; done` ; \
        cd $$dir && $(MAKE) -f ../Makefile $$tgts
	touch $@

%.recheck: %.in
	@dir=`basename $@ .recheck` ; \
        refs="$$dir/*.ref" ; \
	for r in $$refs ; do \
	  rep=`basename $$r .ref`".report" ; \
	  if ! grep -qs "succeeded\|passed" $$dir/$$rep ; then \
	    target=` basename $$r .ref` ; \
	    $(MAKE) $*".in" && rm -f $*".in" ; \
	    rm -f $$dir/$$target".out" $$dir/$$target".report" ; \
	    cd $$dir && $(MAKE) -f ../Makefile $$target".out" $$target".report" && cd - ; \
	  fi \
	done
	touch $*".check"

if USE_MPI
  mpirun_message=mpirun
else
  mpirun_message=
endif

# Avoid copying in dist the builddir files.
distdir: $(DISTFILES)
	@srcdirstrip=`echo "$(srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	topsrcdirstrip=`echo "$(top_srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	list='$(DISTFILES)'; \
	  dist_files=`for file in $$list; do echo $$file; done | \
	  sed -e "s|^$$srcdirstrip/||;t" \
	      -e "s|^$$topsrcdirstrip/|$(top_builddir)/|;t"`; \
	case $$dist_files in \
	  */*) $(MKDIR_P) `echo "$$dist_files" | \
			   sed '/\//!d;s|^|$(distdir)/|;s,/[^/]*$$,,' | \
			   sort -u` ;; \
	esac; \
	for file in $$dist_files; do \
	  d=$(srcdir); \
	  if test -d $$d/$$file; then \
	    dir=`echo "/$$file" | sed -e 's,/[^/]*$$,,'`; \
	    if test -d "$(distdir)/$$file"; then \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    if test -d $(srcdir)/$$file && test $$d != $(srcdir); then \
	      cp -fpR $(srcdir)/$$file "$(distdir)$$dir" || exit 1; \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    cp -fpR $$d/$$file "$(distdir)$$dir" || exit 1; \
	  else \
	    test -f "$(distdir)/$$file" \
	    || cp -p $$d/$$file "$(distdir)/$$file" \
	    || exit 1; \
	  fi; \
	done
