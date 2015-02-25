#!/bin/bash

#These modules/includes do not have to be processed
uses_except="BigDFT_API xc_f90_types_m libxc_funcs_m xc_f90_lib_m \
             m_libpaw_libxc_funcs libxc_functionals \
             mpi omp_lib ifcore f90_unix_proc \
             fox_sax ieee_exceptions memory_profiling \
             netcdf etsf_io etsf_io_low_level yaml_output"
includes_except="fexcp.h"

#Some "use module" statements might be expressed as cpp macros
n_macros=3
macros_in=( "USE_DEFS" "USE_MPI_WRAPPERS" "USE_MSG_HANDLING" )
macros_out=( "m_libpaw_defs" "m_libpaw_mpi" "m_libpaw_tools" )

function adddep
{
 file=$1
 uses=`grep -e "^ *use " $file | sed "s/^ *use *\([a-zA-Z0-9_]*\).*/\1/g" | sort | uniq`
 for ((i=0 ; ${n_macros} - $i ; i++));do
  madd=`grep -e "^ *"${macros_in[$i]} ${file} | sed "s/${macros_in[$i]}/${macros_out[$i]}/g"`
  if test -n "${madd}";then uses=`echo -e ${uses}" \n"${madd}`;fi
 done
 includes=`grep -e "^ *include " $file | sed "s/^ *include [\"']*\([a-zA-Z0-9_.]*\)[\"'].*/\1/g" | sort | uniq`
 if test -n "$uses" -o -n "$includes" ; then
	file=${file##./}
	if test -n "$uses" ; then
		for use in $uses ; do
			out="True"
			for e in $uses_except ; do
				if test x"$use" == x"$e" ; then
					out="False"
					break
				fi
			done
			if test x"$out" = x"True" ; then
				if test x"$use" = x"libxc_functionals" ; then
				use="m_"$use
				fi
				echo -n " "$use".o"
			fi
		done
	fi
	if test -n "$includes" ; then
		for incl in $includes ; do
			out="True"
			for e in $includes_except ; do
				if test x"$incl" == x"$e" ; then
					out="False"
					break
				fi
			done
			if test x"$out" = x"True" ; then
				if test "$incl" = "mpif.h" ; then
					incl='$(mpi_include)'
				else
					incl=${file%/*}'/'$incl
					incls=`adddep $incl`
					if test -n "$incls" ; then
						echo -n " "$incls
					fi
				fi
				echo -n " "$incl
			fi
		done
	fi
    fi
}

for file in `find . -name \*.[fF]90` ; do
    incls=`adddep $file`
    file=${file##./}
    base=`basename $file`
    obj=${base%.[fF]90}
    ext=${base#$obj}
    obj=$obj".o"
    echo -n $obj": "$file
    if test "$ext" == ".F90" ; then
	compile='$(PPFCCOMPILE)'
    else
	compile='$(FCCOMPILE)'
    fi
    for incl in $incls ; do
	echo " \\"
	echo -n "	"$incl
    done
    echo
    echo "	$compile -c -o $obj"' `test -f '"'$file' || echo '\$(srcdir)/'"'`'$file
    echo
done
