#!/bin/sh

# Script used to generate configure script from directives.
echo "Listing known macro with 'aclocal'."
aclocal -I config/m4
echo "Creating configure script  with 'autoconf'."
autoconf
echo "Creating required files for autotools."
automake --add-missing --copy
echo "Autotoolize the libXC source tree."
cd libXC; autoreconf; cd -
