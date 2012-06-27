#!/bin/sh

# Script used to generate configure script from directives.
echo "Listing known macro with 'aclocal'."
aclocal -I config/m4
echo "Creating configure script  with 'autoconf'."
autoconf
echo "Creating config.h.in with 'autoheader'."
autoheader
echo "Creating required files for autotools."
automake --add-missing --copy
echo "Generating PSP database."
python config/pspconf.py > src/init/psp.inc
sed '/!!PSP_TABLE!!/r src/init/psp.inc' src/init/pspconf.in.f90 > src/init/pspconf.f90
rm -f src/init/psp.inc
echo "Autotoolize the libXC source tree."
cd libxc-1.1.0; libtoolize; autoreconf; cd -
echo "Autotoolize the libyaml source tree."
cd yaml-0.1.4; libtoolize -f; autoreconf -i; cd -
echo "Autotoolize the S_GPU source tree."
cd S_GPU; ./autogen.sh; cd -

