## Process this file with automake to produce Makefile.in
chks=wrappers.out.report.yaml posinp.report.yaml $(py_report) $(babel_report)

if BUILD_BINDINGS
SHORT_CHECK = $(chks)
LONG_CHECK = $(SHORT_CHECK)
else
SHORT_CHECK = 
LONG_CHECK = $(chks)
endif

include $(srcdir)/../check.mk


SUBDIRS =

check_PROGRAMS = check_wrappers hooks libconv $(babel_prg)

AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src -Wl,-rpath=$(abs_top_builddir)/src/modules -Wl,-rpath=$(DESTDIR)$(libdir)
bigdft_library = $(top_builddir)/src/libbigdft-1.so.@BIGDFT_MINOR_VERSION@
else
bigdft_library = $(top_builddir)/src/libbigdft-1.a @LIB_BIGDFT_DEPS@
endif

if HAVE_OPENBABEL
babel_prg = babel
babel_report = babel.out.report.yaml
bigdft_library += $(top_builddir)/src/modules/libbabel.so.@BIGDFT_MINOR_VERSION@
if BUILD_DYNAMIC_LIBS
else
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src/modules
endif
endif


check_wrappers_SOURCES = check_wrappers.c
check_wrappers_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/src/bindings -I$(top_builddir)/src/bindings @GLIB_CFLAGS@
check_wrappers_LDADD = $(bigdft_library) @FCLIBS@

hooks_SOURCES = hooks.f90
hooks_FCFLAGS = -I$(top_builddir)/src/orbitals -I$(top_builddir)/src/modules @LIB_FUTILE_CFLAGS@
hooks_LDADD = $(bigdft_library) @LIB_FUTILE_LIBS@

libconv_SOURCES = libconv.f90
libconv_FCFLAGS = -I$(top_builddir)/src/orbitals -I$(top_builddir)/src/modules @LIB_FUTILE_CFLAGS@
libconv_LDADD = $(bigdft_library) @LIB_FUTILE_LIBS@


babel_SOURCES = babel.f90
babel_FCFLAGS = @LIB_FUTILE_CFLAGS@
babel_LDADD = $(bigdft_library) @LIB_FUTILE_LIBS@


EXTRA_DIST = wrappers.out.ref.yaml posinp.ref.yaml pywrappers.out.ref.yaml check_wrappers.py.in babel.out.ref.yaml TiO2.cif
CLEANFILES = default.* posinp.yaml time.yaml wrappers.out.yaml pywrappers.out.yaml *report.yaml input_minimal.yaml input.yaml babel.out.* outfile.xyz \
			 report report_remarks
if WITH_GOBJECT_INTROSPECTION
py_report = pywrappers.out.report.yaml
endif

if BUILD_DYNAMIC_LIBS
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/src
endif

#check_lib_c: $(top_builddir)/src/libbigdft.a

$(EXTRA_DIST):
	@cp -p $(srcdir)/$@ .

posinp.yaml: wrappers.out.yaml

wrappers.out.yaml: check_wrappers
	-$(run_parallel) ./check_wrappers | grep -v "Unable to read mpd.hosts" > $@

check_wrappers.py: check_wrappers.py.in
	sed							\
		-e s!\@prefix\@!$(prefix)!			\
		-e s!\@exec_prefix\@!$(exec_prefix)!		\
		-e s!\@pyexecdir\@!$(pyexecdir)!		\
		< $< > $@

pywrappers.out.yaml: check_wrappers.py
	-$(run_parallel) python check_wrappers.py | grep -v "Unable to read mpd.hosts" > $@

babel.out.yaml: babel TiO2.cif
	-./babel -i TiO2.cif > $@

%.report.yaml: %.ref.yaml %.yaml
	name=`basename $< .ref.yaml | sed "s/[^_]*_\?\(.*\)[^.].*/\1/"` ;\
	tols=`basename $< .ref.yaml` ;\
	if test -n "$$name" ; then log=log-$$name.yaml ; else log=$$tols.yaml ; fi ;\
	echo "python $(pythondir)/fldiff_yaml.py -r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@" ; \
	python $(pythondir)/fldiff_yaml.py \
		-r $< -d $$log -t $(abs_top_srcdir)/tests/tols-BigDFT.yaml --label=$$tols -o $@ 
