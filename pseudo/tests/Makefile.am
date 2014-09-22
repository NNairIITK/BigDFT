## Process this file with automake to produce Makefile.in

SUBDIRS = \
	  atoms \
	  H-pseudo \
	  C-pseudo \
	  pawpatchtest \
	  tutorials

EXTRA_DIST = \
	README \
	example.O.new \
	pawpatchtest \
	report.py.in

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report.py: report.py.in
	sed							\
		-e s!\@prefix\@!$(prefix)!			\
		-e s!\@exec_prefix\@!$(exec_prefix)!		\
		-e s!\@pyexecdir\@!$(pyexecdir)!		\
		-e s!\@AX_PYYAML_PATH\@!$(AX_PYYAML_PATH)!	\
		< $< > $@

report: 
	python $(top_builddir)/tests/report.py
