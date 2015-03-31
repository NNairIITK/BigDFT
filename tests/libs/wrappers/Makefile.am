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

EXTRA_DIST = wrappers.out.ref.yaml posinp.ref.yaml
CLEANFILES = default.* posinp.yaml time.yaml wrappers.out.yaml *report.yaml input_minimal.yaml input.yaml \
			 report report_remarks

check: wrappers.out.report.yaml posinp.report.yaml

if BUILD_DYNAMIC_LIBS
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/src
endif

#check_lib_c: $(top_builddir)/src/libbigdft.a

$(EXTRA_DIST):
	@cp -p $(srcdir)/$@ .

posinp.yaml: wrappers.out.yaml

wrappers.out.yaml: check_wrappers
	-$(run_parallel) ./check_wrappers | grep -v "Unable to read mpd.hosts" > $@

if BUILD_LIBYAML
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/yaml-0.1.4/src/.libs
PYTHONPATH := ${PYTHONPATH}:`ls -d $(abs_top_builddir)/PyYAML-3.10/build/lib.*`
endif

%.report.yaml: %.ref.yaml %.yaml
	@export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	name=`basename $< .ref.yaml | sed "s/[^_]*_\?\(.*\)[^.].*/\1/"` ;\
	tols=`basename $< .ref.yaml` ;\
	if test -n "$$name" ; then log=log-$$name.yaml ; else log=$$tols.yaml ; fi ;\
	echo "python $(abs_top_srcdir)/tests/fldiff_yaml.py -r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@" ; \
	python $(abs_top_srcdir)/tests/fldiff_yaml.py \
		-r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@ 
