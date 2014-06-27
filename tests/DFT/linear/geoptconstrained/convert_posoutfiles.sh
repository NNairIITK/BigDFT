#!/bin/bash

cd data
for i in *.int
do
	file=`echo "$i" | sed 's/_/-/g'`
	mv $i $file
	outfile=$(basename $file .int).xyz
	../../../../../src/bigdft-tool -a transform-coordinates --infile=$file --outfile=$outfile
	file=`echo "$outfile" | sed 's/-/_/g'`
	mv $outfile $file
done	
cd  -
