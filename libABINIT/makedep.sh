#!/bin/sh

except="BigDFT_API xc_f90_types_m libxc_funcs_m xc_f90_lib_m"

for file in `find . -name \*.[fF]90` ; do
    uses=`grep -e "^ *use " $file | sed "s/^ *use *\([a-zA-Z0-9_]*\).*/\1/g" | sort | uniq`
    includes=`grep -e "^ *include " $file | sed "s/^ *include [\"']*\([a-zA-Z0-9_.]*\)[\"'].*/\1/g" | sort | uniq`
    if test -n "$uses" -o -n "$includes" ; then
	file=${file##./}
        echo -n $file":"
	if test -n "$uses" ; then
	    for use in $uses ; do
		out="True"
		for e in $except ; do
		    if test x"$use" == x"$e" ; then
			out="False"
			break
		    fi
		done
		if test x"$out" = x"True" ; then
		    if test x"$use" = x"libxc_functionals" ; then
			use="m_"$use
		    fi
		    echo " \\"
		    echo -n "	"$use".o"
		fi
	    done
	fi
	if test -n "$includes" ; then
	    for incl in $includes ; do
		if test "$incl" = "mpif.h" ; then
		    incl='$(mpi_include)'
		else
		    incl=${file%/*}'/'$incl
		fi
		echo " \\"
		echo -n "	"$incl
	    done
	fi
	echo
	echo
    fi
done

	
