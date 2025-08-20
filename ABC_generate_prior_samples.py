

import numpy as np 
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--outfile' , help = 'Path to output file where prior samples will be written', type = str , required = True)
#parser.add_argument("--k_range" , help = 'Range of prior distribution for k. Enter integer values' , nargs = 2 , type = int , required = True)
#parser.add_argument("--s_range" , help = 'Range of prior distribution for s. Enter value with precision of 0.1 e.g. -0.4, 0.0, 1.1 etc.' , nargs = 2 , type = float , required = True)
#parser.add_argument("--q_values" , help = 'Values of q to choose from. Enter integer values' , nargs = "*" , type = int , required = True)
args = parser.parse_args()



# Read in parameter prior data from ./ABC_input_parameters.txt
with open("./ABC_input_parameters.txt" , 'r') as priors_file:
	for line in priors_file:
		if ('k_range_low' in line):
			k_range_low = int(line.split("=")[-1])

		if ('k_range_high' in line):
			k_range_high = int(line.split("=")[-1])

		if ('s_range_low' in line):
			s_range_low = float(line.split("=")[-1])

		if ('s_range_high' in line):
			s_range_high = float(line.split("=")[-1])

		if ('q_values' in line):
			q_values = [int(val) for val in (line.split("=")[-1])[1:-1].split(" ")]




### Checks on input parameter values
if (k_range_low >= k_range_high):
	print("Invalid range for k.")
	exit(0)

if (s_range_low >= s_range_high):
	print("Invalid range for s.")
	exit(0)

if (len(q_values) == 0) or (len([val for val in q_values if (val <= 0)]) > 0):
	print("Invalid range for q.")
	exit(0)



np.random.seed()

# Parameter ranges
number_of_parameter_sets = 1000

### Ranges for prior distributions (assume uniform priors for all)
# Initial ecDNA copy number, k
k_choices = np.arange(k_range_low , k_range_high+1 , 1)


# ecDNA selection strength, s
s_choices = np.arange(round(s_range_low , 1) , round(s_range_high , 1)+0.1 , 0.1)
s_choices = [round(val , 1) for val in s_choices]


# Cell pushing strength, q
q_choices = q_values





# Sample parameter triplets from their respective prior ditsributions
s_params = np.random.choice(s_choices , size = number_of_parameter_sets)
k_params = np.random.choice(k_choices , size = number_of_parameter_sets)
q_params = np.random.choice(q_choices , size = number_of_parameter_sets)

for i in range(len(s_params)):
	if (s_params[i].is_integer()):
		s_params[i] = int(s_params[i])

# Write these parameter triplets to an output file 
with open(args.outfile , 'w') as outFile:
	for i in range(number_of_parameter_sets):
		if (s_params[i].is_integer()):
			outFile.write("{},{},{}\n".format(int(s_params[i]) , int(k_params[i]) , q_params[i]))
		else:
			outFile.write("{},{},{}\n".format(s_params[i] , int(k_params[i]) , q_params[i]))
