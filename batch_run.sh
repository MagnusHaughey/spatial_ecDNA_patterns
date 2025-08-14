#!/bin/bash

# Set variables from external file
source batch_run_parameters.txt

# Compile simulation script and run batch of simulations
if g++ -o ./spatial_ecDNA_patterns_2D ./spatial_ecDNA_patterns_2D.cpp
then
	for k in $all_k
	do
		for s in $all_s
		do
			for q in $all_q
			do
				printf "Running %s simulations with (n, k, s, q) = (%s, %s, %s, %s)\n" $repeats $n $k $s $q
				for x in $( seq 1 $repeats )
				do
					./spatial_ecDNA_patterns_2D -n $n -k $k -s $s -q $q -x $x
				done
			done
		done
	done
fi

