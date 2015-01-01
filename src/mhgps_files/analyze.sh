#!/bin/bash


file=$1_tmp
grep "(MHGPS)" $1 | grep -v genv > $file


echo "Average number of energy/force calls per minimization:"
grep "Needed energy calls" $file | awk '{sum+=$10} END { print "Average = ",sum/NR,NR}'

echo "Average number of energy/force calls per saddle computation:"
grep convergence $file | awk '{sum+=$5} END { print "Average = ",sum/NR,NR}'

echo "Average number of energy/force calls per saddle computation including input guess:"
enersaddle=`grep convergence $file | awk '{sum+=$5} END { print sum}'`
energuess=`grep "Energy evaluations for TS guess" $file | awk '{sum+=$8} END { print sum}'`
nsaddle=`grep convergence $file | awk '{sum+=$5} END { print NR}'`
printf '%0.f\n' `echo "( $enersaddle + $energuess )/ $nsaddle " | bc -l`

echo "Average number of intermediate TS:"
grep "(MHGPS) succesfully connected, intermediate transition states" $file  | awk '{sum+=$7} END { print "Average = ",sum/NR,NR}'

echo "Number of successfully connected input minima pairs: "
grep "succesfully connected" $file | wc -l

echo "Number of failed connections of input minima pairs: "
grep establish $file | wc -l
