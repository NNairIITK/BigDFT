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
if test -f config/pspconf.py ; then
  python config/pspconf.py > src/init/psp.inc
  sed '/!!PSP_TABLE!!/r src/init/psp.inc' src/init/pspconf.in.f90 > src/init/pspconf.f90
  rm -f src/init/psp.inc
else
  echo " WARNING, missing config/pspconf.py, cannot generate built-in pseudo-potentials."
fi
GDBUS_CODEGEN=`which gdbus-codegen`
if test -n "$GDBUS_CODEGEN" ; then
  echo "Generate Dbus bindings (obsolete)."
  cd src/bindings
  $GDBUS_CODEGEN --interface-prefix eu.etsf.bigdft.dbus. --generate-c-code bindings_dbus --c-namespace BigdftDBus --c-generate-object-manager bindings_dbus.xml
  cd -
fi
echo "Autotoolize the libXC source tree."
cd libxc-2.0.x; glibtoolize -fc; autoreconf -fi; cd -
echo "Autotoolize the libyaml source tree."
cd yaml-0.1.4; glibtoolize -fc; autoreconf -fi; cd -
echo "Autotoolize the S_GPU source tree."
cd S_GPU; ./autogen.sh; cd -

