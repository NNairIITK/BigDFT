scrp_PYTHON= \
	__init__.py \
	DeltaTest.py

scrpdir= $(pythondir)/BigDFT/scripts

EXTRA_DIST= CIFs/

install-data-hook: CIFs
	cp -rpL CIFs/ '$(DESTDIR)$(scrpdir)'

uninstall-local:
	rm -rf '$(DESTDIR)$(scrpdir)'/CIFs
