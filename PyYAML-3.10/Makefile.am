if BUILD_LIBYAML
pyyaml_clean = pyyaml_clean
pyyaml_build = pyyaml_build
else
pyyaml_clean =
pyyaml_build =
endif

pyyaml_build:
	rm -rf build ; \
	echo "[build_ext]" > setup.cfg ; \
	echo "define=YAML_DECLARE_STATIC" >> setup.cfg ; \
	echo "include_dirs=$(abs_top_srcdir)/yaml-0.1.4/include" >> setup.cfg ; \
	echo "library_dirs=$(abs_top_builddir)/yaml-0.1.4/src/.libs" >> setup.cfg ; \
	python setup.py build

pyyaml_clean:
	python setup.py clean

all: $(pyyaml_build)

clean: $(pyyaml_clean)
