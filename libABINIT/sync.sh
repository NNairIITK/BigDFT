#!/bin/sh

ABINIT=/local/caliste/abinit-6.0.0-caliste

for file in `find .` ; do
    tmp=`basename $file`
    abfile=`find $ABINIT -name $tmp`
    if test -n "$abfile" ; then
	echo $file "<-" $abfile
	cp $abfile $file
    fi
done

patch -p0 < patch
