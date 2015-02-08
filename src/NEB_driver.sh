#!/bin/bash --noprofile

###############################################
# Copyright (C) 2009-2011 Damien Caliste (DC) #
# This file is distributed under the terms    #
# of the GNU General Public License.          #
# See http://www.gnu.org/copyleft/gpl.txt .   #
###############################################

# Compute the forces for all the current replica.
# INPUTS:
#   [all free_only] [job_name] [workdir] [first_config]
# OUTPUTS:
#   [gen_output_file]
DEBUG="yes"

INCLUDED="yes"
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
if [ ! -d $workdir ] ; then
  mkdir -p $wordir
fi
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

# Try to find the iteration number.
cd $datadir
if [ -f $job_name.NEB.tar ] ; then
    i=`tar -tf $job_name.NEB.tar | wc -l`
    i=$(($i / 2))
else
    i=0
fi
neb_iter=`printf "%04d" $i`
if test x"$DEBUG" != x ; then
    echo "Current iter is "${neb_iter}"."
fi

# Save and remove the possible output from previous run.
cd $datadir
if [ -f gen_output_file ] ; then
    j=`printf "%04d" $((${i}-1))`
    cp -f -p gen_output_file $job_name.NEB.it${j}.forces
    tar -rf $job_name.NEB.tar --remove-files $job_name.NEB.it${j}.forces
fi
rm -f gen_output_file

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
    rm -f $dir/OK $dir/FAILED $dir/FAILURES
    if [ ! -f $dir/RESTART ] ; then
        touch $dir/START
    fi
    neb_dir[$count]=$dir
done

# Run an init function if necessary.
cd $workdir
outfile=`init_jobs $job_name $datadir $count $n_nodes $first_config | tail -n 1`

# Main loop. Start run in each directories, erasing the START keyword.
jobs_done=0
iter=0
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
	    run_job $job_name $datadir $count
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
	    if test x"$DEBUG" != x ; then
		echo "Job "${count}" finished ("$(($jobs_done + 1))"/"$((${max} - ${min} + 1))")."
	    fi

	    if test x"$res" == x"1" ; then
		jobs_done=$(($jobs_done + 1))
		touch OK
	    fi
	    if test x"$res" == x"2" ; then
		jobs_done=$(($jobs_done + 1))
		touch FAILED
	    fi
	    if test x"$res" == x"3" ; then
                if [ ! -f FAILURES ] ; then
                    echo -n "1" > FAILURES
                else
                    echo -n $((`cat FAILURES` + 1)) > FAILURES
                fi
                if test -f FAILURES -a `cat FAILURES` -gt 3 ; then
		    jobs_done=$(($jobs_done + 1))
		    touch FAILED
                else
		    touch START
                fi
	    fi
	fi

    done

    cd $workdir
    wait_jobs $iter

    if [ ${jobs_done} -lt $((${max} - ${min} + 1)) ] ; then
	sleep 10s 2> /dev/null
    fi
    iter=$(($iter + 1))
done

# Call a finalise script.
cd $workdir
finalise_jobs $job_name

# Generate the force file.
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
    rm -f OK FAILED
done

# Copy of the outputs
for ((count=${min};count<=${max};count++)) ; do
    ch=`printf "%02d" $count`
    cd $workdir/$job_name.NEB.$ch
    cp -f -p $outfile $job_name.NEB.it${neb_iter}.out
done

# Copy of previous dat file.
cd $datadir
if [ -f $job_name.NEB.dat ] ; then
    cp -f -p $job_name.NEB.dat $job_name.NEB.it${neb_iter}.dat
    tar -rf $job_name.NEB.tar --remove-files $job_name.NEB.it${neb_iter}.dat
fi

