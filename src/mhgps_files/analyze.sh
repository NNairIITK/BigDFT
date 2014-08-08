#!/bin/bash

echo "Average number of energy/force calls per minimization:"
grep "energy calls" $1 | awk '{sum+=$10} END { print "Average = ",sum/NR,NR}'

echo "Average number of energy/force calls per saddle computation:"
grep convergence $1 | awk '{sum+=$5} END { print "Average = ",sum/NR,NR}'

echo "Number of successfully connected input minima: "
grep connected $1 | wc -l

echo "Number of failed connections: "
grep establish $1 | wc -l
