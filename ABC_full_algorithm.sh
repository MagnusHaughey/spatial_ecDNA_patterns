#!/bin/bash


seed=0
dseed=1
N=10000

# Parameters for parameter priors
# k_range_low=1
# k_range_high=151
# s_range_low=-0.5
# s_range_high=5
# q_values='1 2 5 10 50 1000'
source ./ABC_input_parameters.txt


# Ensure directory for priors exists
if [ ! -d ./ABC_prior_samples/ ] 
then
	mkdir ./ABC_prior_samples
fi


# First, sample parameters from the prior distribution
paramFile='./ABC_prior_samples/prior_'"$seed"'.dat'
PRIOR_GENERATOR_SCRIPT=$PWD'/ABC_generate_prior_samples.py'
if [ -f $PRIOR_GENERATOR_SCRIPT ]
then
	python3 $PRIOR_GENERATOR_SCRIPT --outfile $paramFile --k_range $k_range_low $k_range_high --s_range $s_range_low $s_range_high --q_values $q_values
else
	>&2 echo "ERROR: ABC_generate_prior_samples.py could not be found."
fi



# Compile simulation script and run batch of simulations
SIMULATION_SCRIPT=$PWD'/spatial_ecDNA_patterns_2D.cpp'
SIMULATION_SCRIPT_EXE=${SIMULATION_SCRIPT%.cpp}
PATIENT_DATA_PATH=$PWD'/GB_UK_ecDNA_copy_number_distributions/*/'
OUTPUT_FILE='./ABC_raw_results/raw_results_'"$seed"'.dat'

if [ ! -d ./ABC_raw_results/ ] 
then
	mkdir ./ABC_raw_results/
fi


# Check simulation file can be found
if [ -f $SIMULATION_SCRIPT ]
then
	# Check simulation compiles
	if g++ -o $SIMULATION_SCRIPT_EXE $SIMULATION_SCRIPT
	then
		# Loop over parameter values & simulate spatial data
		while IFS=, read -r s k q
		do
		 
			# Simulate spatial data
			$SIMULATION_SCRIPT_EXE -n $N -s $s -k $k -q $q -x $seed


			# Simulation data path
			parentPath=$PWD'/results/n='"$N"'_k='"$k"'_q='"$q"'_s='"$s"'/seed='"$seed"'/'
			dataPath=$PWD'/results/n='"$N"'_k='"$k"'_q='"$q"'_s='"$s"'/seed='"$seed"'/tissue.csv'



			# Sample and analyse spatial data, writing simulation-patient distance to an output file
			if [ -f $dataPath ]
			then
				for patient in $PATIENT_DATA_PATH
				do
					python3 ABC_rejectionSampling.py --simulation $dataPath --patient $patient >> $OUTPUT_FILE
				done
			fi


			# Remove simulation data
			rm -r $parentPath


			# Increase seed
			seed=$(( seed + dseed ))

		done < $paramFile
	else
		>&2 echo "ERROR: spatial_ecDNA_patterns_2D.cpp failed compilation."
	fi
else
	>&2 echo "ERROR: spatial_ecDNA_patterns_2D.cpp could not be found."
fi






