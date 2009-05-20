#!/bin/bash --noprofile

#############################################
# Copyright (C) 2009 Damien Caliste (DC)    #
# This file is distributed under the terms  #
# of the GNU General Public License.        #
# See http://www.gnu.org/copyleft/gpl.txt . #
#############################################

# Compute the forces for all the current replica.
# INPUTS:
#   [all free_only] [job_name] [workdir] [first_config]
# OUTPUTS:
#   [gen_output_file]
DEBUG=

source NEB_include.sh

type=$1
job_name=$2
restart_file=$job_name".NEB.restart"
n_images=`grep Replica $restart_file | wc -l`
datadir=$PWD
case $3 in
     /*) workdir=$3 ;;
     *) workdir=$PWD/$3 ;;
esac
n_nodes=$((`wc -l $restart_file | cut -d' ' -f1` / $n_images - 2))
case $4 in
     /*) first_config=$4 ;;
     *) first_config=$PWD/$4 ;;
esac

if test x"$DEBUG" != x ; then
    echo "External computation of forces for replicas."
    echo " type     = "$type
    echo " job_name = "$job_name
    echo " n_images = "$n_images
    echo " n_nodes  = "$n_nodes
    echo " datadir  = "$datadir
    echo " workdir  = "$workdir
fi

# Consistency check.
if test "$n_images" -lt 0 ; then
    exit
fi
if [ ! -d "$workdir" ] ; then
    exit
fi

# Set the ids of replica to move.
if [ $type == "all" ]; then
    min=1
    max=$n_images
elif [ $type == "free_only" ]; then
    min=2
    max=$(( $n_images - 1 ))
fi

# Make the computation directories if they don't exist
# and set up the flags.
for ((count=${min};count<=${max};count++)) ; do
    ch=`printf "%02d" $count`
    dir="$workdir/$job_name.NEB.${ch}"
    if [ ! -d $dir ] ;then
	mkdir $dir
    fi
    # Empty the working directory.
    if [ ! -f $dir/RESTART ] ; then
        rm -f $dir/*
        touch $dir/START
    fi
    neb_dir[$count]=$dir
done

# Main loop. Start run in each directories, erasing the START keyword.
jobs_done=0
while [ ${jobs_done} -lt $((${max} - ${min} + 1)) ] ; do

    # Do things in each directories depending on keywords.
    for ((count=${min};count<=${max};count++)) ; do

	cd ${neb_dir[$count]}
	# keyword START
	if [ -f START ] ; then
	    if test x"$DEBUG" != x ; then
		echo "Start job "${count}"."
	    fi

	    # Create the input
	    make_input $job_name $datadir $count $n_nodes $first_config

	    # Erase the keyword
	    rm -f START

	    # Start it
	    run_job $job_name $count
	fi

	cd ${neb_dir[$count]}
	# keyword RESTART
	if [ -f RESTART ] ; then
	    rm -f RESTART
	fi

	cd ${neb_dir[$count]}
	res=`check_job $job_name`
	# Job has finished
	if test $res -gt 0 && [ ! -f OK ] && [ ! -f FAILED ] ; then
	    jobs_done=$(($jobs_done + 1))
	    if test x"$DEBUG" != x ; then
		echo "Job "${count}" finished ("$jobs_done"/"$((${max} - ${min} + 1))")."
	    fi

	    if test x"$res" == x"1" ; then
		touch OK
	    fi
	    if test x"$res" == x"2" ; then
		touch FAILED
	    fi
	fi

    done

    if [ ${jobs_done} -lt $((${max} - ${min} + 1)) ] ; then
	sleep 10s
    fi
done

# Compression of the data
cd $workdir
for ((i = 0; i < 256; i++)) ; do
    ch=`printf "%03d" $i`
    if ! [ -f $job_name.NEB.it${ch}.tar.bz2 ] ; then
	if test x"$DEBUG" != x ; then
	    echo "Compression of replica calculations into '$job_name.NEB.it${ch}.tar.bz2'"
	fi
	tar --exclude \*.bz2 -cjf $job_name.NEB.it${ch}.tar.bz2 $job_name.NEB.*
	break
    fi
done

# Generate the force file.
rm -f $datadir/gen_output_file
for ((count=${min};count<=${max};count++)) ; do
    cd ${neb_dir[$count]}

    forces=`grep_forces $n_nodes`
    
    if [ -f OK ] && test -n "$forces" ; then
	grep_forces $n_nodes >> $datadir/gen_output_file
    else
	echo " WARNING: image ${count} not converged" 
	echo " THE PROGRAM IS GOING TO STOP !!!"
	echo "  9999999.99999999" >> $datadir/gen_output_file
    fi
done
