

import numpy as np 
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--outfile' , help = 'Path to output file where prior samples will be written', type = str , required = True)
parser.add_argument("--k_range" , help = 'Range of prior distribution for k. Enter integer values' , nargs = 2 , type = int , required = True)
parser.add_argument("--s_range" , help = 'Range of prior distribution for s. Enter value with precision of 0.1 e.g. -0.4, 0.0, 1.1 etc.' , nargs = 2 , type = float , required = True)
parser.add_argument("--q_values" , help = 'Values of q to choose from. Enter integer values' , nargs = 1 , type = str , required = True)
args = parser.parse_args()



# Parse q_values from str to int, as inputting int at command line was causing an error
args.q_values = [int(val) for val in args.q_values[0].split(" ")]



### Checks on input parameter values
if (args.k_range[0] >= args.k_range[1]):
	print("Invalid range for k.")
	exit(0)

if (args.s_range[0] >= args.s_range[1]):
	print("Invalid range for s.")
	exit(0)

if (len(args.q_values) == 0) or (len([val for val in args.q_values if (val <= 0)]) > 0):
	print("Invalid range for q.")
	exit(0)



np.random.seed()

# Parameter ranges
number_of_parameter_sets = 1000

### Ranges for prior distributions (assume uniform priors for all)
# Initial ecDNA copy number, k
k_choices = np.arange(args.k_range[0] , args.k_range[1]+1 , 1)


# ecDNA selection strength, s
s_choices = np.arange(round(args.s_range[0] , 1) , round(args.s_range[1] , 1)+0.1 , 0.1)
s_choices = [round(val , 1) for val in s_choices]


# Cell pushing strength, q
q_choices = args.q_values





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
