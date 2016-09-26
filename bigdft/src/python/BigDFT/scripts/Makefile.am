scrp_PYTHON= \
	__init__.py \
	DeltaTest.py

scrpdir= $(pythondir)/BigDFT/scripts

EXTRA_DIST= CIFs/

install-data-hook: CIFs
	cp -rpL CIFs/ '$(scrpdir)'

uninstall-local:
	rm -rf '$(scrpdir)'/CIFs
