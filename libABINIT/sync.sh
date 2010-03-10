#!/bin/sh

#ABINIT=/local/caliste/abinit-6.0.0-caliste
ABINIT=/local/deutsch/L_Sim/ABINIT/abinit-deutsch-6.1.0-private

for file in `find .` ; do
    tmp=`basename $file`
    abfile=`find $ABINIT -name $tmp`
    if test -n "$abfile" -a "$tmp" != "moldyn.F90" ; then
        echo $file "<-" $abfile
        cp $abfile $file
    fi
done

patch -p0 < patch
