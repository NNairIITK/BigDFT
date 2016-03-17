## Process this file with automake to produce Makefile.in

SUBDIRS = fft

EXTRA_DIST = \
	tols.yaml \
	fldiff_yaml.py.in \
	report.py.in \
	yaml_test.ref.yaml $(OTHERS)

OTHERS = yaml_test.f90

check_PROGRAMS = yaml_test
python_PYTHON = report.py fldiff.py fldiff_yaml.py yaml_hl.py
python_DATA = yaml_hl.cfg yaml_hl_remarks.cfg

AM_FCFLAGS = -I../flib
AM_LDFLAGS = -L$(top_builddir)/src
#LIBS += -lflib -lrt @LIB_YAML_LIBS@


if BUILD_DYNAMIC_LIBS
flib_libraries = -lfutile-1
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
else
flib_libraries = $(top_builddir)/src/libfutile-1.a
LIBS += @LIB_YAML_LIBS@
endif

yaml_test_SOURCES = yaml_test.f90 errs.f90 yamlout.f90 dicts.f90 dynmem.f90 yaml_invoice.f90 utls.f90
yaml_test_LDADD = $(flib_libraries)

CLEANFILES = yaml_test.out.yaml yaml_test.report.yaml fldiff_yaml.py report.py

yaml_test.o: $(abs_top_builddir)/src/libfutile-1.a

yaml_test.out.yaml: yaml_test
	./yaml_test > yaml_test.out.yaml

yaml_test.report.yaml: fldiff_yaml.py tols.yaml yaml_test.ref.yaml yaml_test.out.yaml
	python fldiff_yaml.py --label yaml_test -t $(srcdir)/tols.yaml 	\
		-r $(srcdir)/yaml_test.ref.yaml -d yaml_test.out.yaml --output $@

check: yaml_test.report.yaml

report.py: report.py.in
	@sed							\
		-e s!\@prefix\@!$(prefix)!			\
		-e s!\@exec_prefix\@!$(exec_prefix)!		\
		-e s!\@pyexecdir\@!$(pyexecdir)!		\
		-e s!\@AX_PYYAML_PATH\@!$(AX_PYYAML_PATH)!	\
		-e s!\@abs_srcdir\@!$(abs_srcdir)!		\
		< $< > $@

fldiff_yaml.py: fldiff_yaml.py.in
	@sed							\
		-e s!\@prefix\@!$(prefix)!			\
		-e s!\@exec_prefix\@!$(exec_prefix)!		\
		-e s!\@pyexecdir\@!$(pyexecdir)!		\
		-e s!\@AX_PYYAML_PATH\@!$(AX_PYYAML_PATH)!	\
		-e s!\@abs_srcdir\@!$(abs_srcdir)!		\
		< $< > $@
