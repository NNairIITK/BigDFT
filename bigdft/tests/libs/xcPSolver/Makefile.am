## Process this file with automake to produce Makefile.in

SHORT_CHECK = PS_Check.out.yaml PS_Check.out.yaml.report
LONG_CHECK = $(SHORT_CHECK)

include $(srcdir)/../check.mk


EXTRA_DIST = PS_Check.ref.yaml tols.yaml

check_PROGRAMS = PS_Check PS_Integral gaussian gaussian2D gaussian3D

AM_FCFLAGS = -I$(top_builddir)/includes -I. @LIBABINIT_INCLUDE@ @LIB_XC_CFLAGS@  @MPI_INCLUDE@
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
endif

PS_Check_SOURCES = PS_Check.f90
PS_Check_LDADD = @LIB_BIGDFT_LIBS@ @LINALG_LIBS@

PS_Integral_SOURCES = PS_Integral.f90 fisf.f90
PS_Integral_LDADD = @LIB_BIGDFT_LIBS@ @LINALG_LIBS@

gaussian_SOURCES = gaussian.f90
gaussian_LDADD = @LIB_BIGDFT_LIBS@ @LINALG_LIBS@

gaussian2D_SOURCES = gaussian2D.f90
gaussian2D_LDADD = @LIB_BIGDFT_LIBS@ @LINALG_LIBS@

gaussian3D_SOURCES = gaussian3D.f90
gaussian3D_LDADD = @LIB_BIGDFT_LIBS@ @LINALG_LIBS@


CLEANFILES = PS_Check.out.yaml fldiff.report.yaml malloc.prc time.prc

PS_Check.out.yaml: PS_Check
	$(run_parallel) ./PS_Check --ndim=[57,48,63] --ixc=0  --geocode=F  > PS_Check.out.yaml
	$(run_parallel) ./PS_Check --ndim=[57,48,63] --ixc=1  --geocode=F >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check --ndim=[57,48,63] --ixc=11 --geocode=F >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check --ndim=[57,48,63] --ixc=13 --geocode=F >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check --ndim=[64,64,64] --ixc=0  --geocode=P >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check --ndim=[64,64,64] --ixc=1  --geocode=P >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check --ndim=[64,64,64] --ixc=11 --geocode=P >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check --ndim=[64,64,64] --ixc=13 --geocode=P >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check --ndim=[32,64,48] --ixc=0  --geocode=S >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check --ndim=[32,64,48] --ixc=1  --geocode=S >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check --ndim=[32,64,48] --ixc=11 --geocode=S >> PS_Check.out.yaml
	$(run_parallel) ./PS_Check --ndim=[32,64,48] --ixc=13 --geocode=S >> PS_Check.out.yaml
	cp -f PS_Check.out.yaml tmp
	cat tmp | grep -v 'Unable to read mpd.hosts' > PS_Check.out.yaml
	rm -f tmp

PS_Check.out.yaml.report: PS_Check.out.yaml
	python $(pythondir)/fldiff_yaml.py \
		-r $(abs_srcdir)/PS_Check.ref.yaml -d PS_Check.out.yaml -t $(abs_srcdir)/tols.yaml \
	        --output fldiff.report.yaml

report: 
	python $(pythondir)/report.py
