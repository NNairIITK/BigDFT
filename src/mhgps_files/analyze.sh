#!/bin/bash

echo "Average number of energy/force calls per minimization:"
grep "energy calls" $1 | awk '{sum+=$10} END { print "Average = ",sum/NR,NR}'

echo "Average number of energy/force calls per saddle computation:"
grep convergence $1 | awk '{sum+=$5} END { print "Average = ",sum/NR,NR}'

echo "Average number of energy/force calls per saddle computation including input guess:"
enersaddle=`grep convergence $1 | awk '{sum+=$5} END { print sum}'`
energuess=`grep "Energy evaluations for TS guess" $1 | awk '{sum+=$8} END { print sum}'`
nsaddle=`grep convergence $1 | awk '{sum+=$5} END { print NR}'`
printf '%0.f\n' `echo "( $enersaddle + $energuess )/ $nsaddle " | bc -l`

echo "Average number of intermediate TS:"
grep "(MHGPS) succesfully connected, intermediate transition states" $1  | awk '{sum+=$7} END { print "Average = ",sum/NR,NR}'

echo "Number of successfully connected input minima pairs: "
grep connected $1 | wc -l

echo "Number of failed connections of input minima pairs:: "
grep establish $1 | wc -l
