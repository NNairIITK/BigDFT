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
echo "Autotoolize the libXC source tree."
cd libXC; autoreconf; cd -
#echo "Autotoolize the S_GPU source tree."
#cd S_GPU; ./autogen.sh; cd -

