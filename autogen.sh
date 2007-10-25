#!/bin/sh

# Script used to generate configure script from directives.
echo "Listing known macro with 'aclocal'."
aclocal -I config/m4
echo "Creating configure script  with 'autoconf'."
autoconf
echo "Create the config.h.in file with 'autoheader'."
autoheader
echo "Creating required files for autotools."
automake --add-missing --copy
