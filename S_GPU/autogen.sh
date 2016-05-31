#!/bin/sh

# Script used to generate configure script from directives.
echo "Listing known macro with 'aclocal' (may take some time)"
aclocal -I m4
echo "Generate config.h.in"
autoheader
echo "Libtool generation"
case `uname` in
Darwin*) glibtoolize --automake ;;
*) libtoolize --automake ;;
esac
echo "Creating configure script  with 'autoconf'."
autoconf
echo "Creating required files for autotools."
automake --add-missing --copy --gnu
