#!/bin/sh

except="BigDFT_API xc_f90_types_m libxc_funcs_m xc_f90_lib_m"

incls=""

function adddep
{
    file=$2
    uses=`grep -e "^ *use " $file | sed "s/^ *use *\([a-zA-Z0-9_]*\).*/\1/g" | sort | uniq`
    includes=`grep -e "^ *include " $file | sed "s/^ *include [\"']*\([a-zA-Z0-9_.]*\)[\"'].*/\1/g" | sort | uniq`
    if test -n "$uses" -o -n "$includes" ; then
	file=${file##./}
	base=`basename $file`
	obj=${base%.[fF]90}
	ext=${base#$obj}
	obj=$obj".o"
	if test "$ext" == ".F90" ; then
	    compile='$(PPFCCOMPILE)'
	else
	    compile='$(FCCOMPILE)'
	fi
	if test $1 == "build" ; then
	    echo -n $obj": "$file
	else
	    echo -n $file":"
	fi
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
		    if test $1 == "build" ; then
			incls=$incls" "$incl
		    fi
		fi
		echo " \\"
		echo -n "	"$incl
	    done
	fi
	echo
	if test $1 == "build" ; then
	    echo "	$compile -c -o $obj"' `test -f '"'$file' || echo '\$(srcdir)/'"'`'$file
	fi
	echo
    fi
}

for file in `find . -name \*.[fF]90` ; do
    adddep "build" $file
done

for incl in $incls ; do
    adddep "dep" $incl
done
