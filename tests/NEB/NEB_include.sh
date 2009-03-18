#!/bin/bash

make_input()
{
    # We create a posinp file from the .NEB.restart file.
    # Set also the input.dat and psp files.
    # INPUTS:
    # $1 the job name
    # $2 the datadir
    # $3 the replica id
    # $4 the number of nodes
    datadir=$2
    restart_file=$datadir/$1".NEB.restart"
    position_file=$datadir/posinp
    n_nodes=$4
    id=$3

    # Create the posinp file
    rm -f posinp
    
    echo $n_nodes" atomic" > posinp
    head -n2 $position_file | tail -n 1 >> posinp
    head -n $(($n_nodes + 2)) $position_file | tail -n $n_nodes | sed "s/^ *//g" | cut -d" " -f1 > atoms
    grep -A $(($n_nodes + 1)) -e "Replica: *$id" $restart_file | tail -n $n_nodes | sed "s/  */ /g" | cut -d' ' -f2,3,4 > positions
    paste -d' ' atoms positions >> posinp
    rm -f atoms positions

    # Copy the input.dat
    cp -f -p $datadir/input.dat .

    # Copy the pseudo-potential files
    cp -f -p $datadir/psppar.* .
}

run_job()
{
    # Create the go.pbs file
    # INPUTS:
    # $1 name of the job
    # $2 the replica id

    cat > go.pbs <<EOF
#!/bin/sh
#PBS -S /bin/bash
#PBS -N $1_$2
#PBS -e $1.qsub.err
#PBS -o $1.qsub.out
#PBS -l select=6:ncpus=8:mpiprocs=8
#PBS -l walltime=24:00:00

cd $PBS_O_WORKDIR
module purge
module load intel
module load intelmpi

## Génération de la hostlist dans un job PBS
cat $PBS_NODEFILE|uniq > mpd.hosts

## récupération du nombre de processus mpi
NPROCS=$(cat $PBS_NODEFILE|wc -l)

##nettoyage
mpdcleanup -f mpd.hosts --rsh=ssh

##Démarrage du demon mpd
mpdboot --rsh=ssh -v -n `cat mpd.hosts|wc -l` -f mpd.hosts

## lancement de l'exécutable
mpiexec -n $NPROCS /work/pochet/bin/cluster > log
echo "job done"

##Arrêt des démons mpd
mpdallexit

EOF

qsub go.pbs
}

check_job()
{
    # INPUTS:
    # $1 name of the job
    # OUPUTS:
    #  0: job running
    # -1: job not start yet
    #  1: job finished and succeeded
    #  2: job finished but failed

    if [ ! -f log ] ; then
	# Job not started yet.
	echo "-1"
	return
    fi

    if grep -qs "job done" $1.qsub.out ; then

	if grep -qs "Final values of the Forces for each atom" log ; then
	    # Success case
	    echo "1"
	    return
	else
	    # Failure
	    echo "2"
	    return
	fi
    fi
    
    # Job not finished
    echo "0"
}

grep_forces()
{
    # Grep the forces from the output file
    # INPUTS:
    # $1 number of nodes to grep forces from
    # OUTPUTS:
    #  the total energy
    #  the list of forces for each atoms
    grep "FINAL iter" log | sed 's/  */ /g' | cut -d' ' -f6
    grep -A $1 "Final values of the Forces for each atom" log | sed 's/  */ /g' | cut -d' ' -f4,5,6 | tail -n $1
}
