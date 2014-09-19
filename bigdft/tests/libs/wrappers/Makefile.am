## Process this file with automake to produce Makefile.in

SUBDIRS =

check_PROGRAMS = check_wrappers

AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
LIBS += @FCLIBS@
else
LIBS += @LIBS_SHORT_DEPS@ @FCLIBS@
endif

check_wrappers_SOURCES = check_wrappers.c
check_wrappers_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/src/modules -I$(top_builddir)/src/modules @GLIB_CFLAGS@
check_wrappers_LDADD = -lbigdft-1

EXTRA_DIST = wrappers.out.ref.yaml posinp.ref.yaml pywrappers.out.ref.yaml
CLEANFILES = default.* posinp.yaml time.yaml wrappers.out.yaml pywrappers.out.yaml *report.yaml input_minimal.yaml input.yaml

if WITH_GOBJECT_INTROSPECTION
py_report = pywrappers.out.report.yaml
endif

check: wrappers.out.report.yaml posinp.report.yaml $(py_report)

if BUILD_DYNAMIC_LIBS
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/src
endif

#check_lib_c: $(top_builddir)/src/libbigdft.a

$(EXTRA_DIST):
	@cp -p $(srcdir)/$@ .

posinp.yaml: wrappers.out.yaml

wrappers.out.yaml: check_wrappers
	-$(run_parallel) ./check_wrappers | grep -v "Unable to read mpd.hosts" > $@

pywrappers.out.yaml: check_wrappers.py
	-$(run_parallel) python check_wrappers.py | grep -v "Unable to read mpd.hosts" > $@

%.report.yaml: %.ref.yaml %.yaml
	name=`basename $< .ref.yaml | sed "s/[^_]*_\?\(.*\)[^.].*/\1/"` ;\
	tols=`basename $< .ref.yaml` ;\
	if test -n "$$name" ; then log=log-$$name.yaml ; else log=$$tols.yaml ; fi ;\
	echo "python $(abs_top_builddir)/tests/fldiff_yaml.py -r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@" ; \
	python $(abs_top_builddir)/tests/fldiff_yaml.py \
		-r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@ 
