


import numpy as np 
import sys 
import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt 
from collections import Counter
import argparse
import glob 
import scipy.stats as stats
from pathlib import Path
from scipy.stats import mode




q_mapping = {
	1: 1,
	2: 2,
	5: 3,
	10: 4,
	50: 5,
	1000: 6
}




# Parse command line argument
parser = argparse.ArgumentParser()
parser.add_argument('--epsilon_core', help = 'Epsilon for core sample' , required = False , type = float , default = 0.25)
parser.add_argument('--epsilon_margin', help = 'Epsilon for margin sample' , required = False , type = float , default = 0.25)
parser.add_argument('--epsilon_leadingEdge', help = 'Epsilon for leading edge sample' , required = False , type = float , default = 0.25)
parser.add_argument('--posterior_size', help = 'Desired size of posterior sample ()' , required = False , type = float , default = 5)
args = parser.parse_args()





### Import data 
all_patient_IDs = []
all_ABC_particle_measurements = []

for file in [(i, raw_ABC_file) for i, raw_ABC_file in enumerate(glob.glob("./ABC_raw_results/*.dat"))]:

	if (file[0] == 0):
		all_patient_IDs = np.genfromtxt(file[1] , unpack = True , delimiter = "," , usecols = (0) , dtype = str)
		all_ABC_particle_measurements = np.loadtxt(file[1] , unpack = True , delimiter = "," , usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14))
	else:
		all_patient_IDs += np.genfromtxt(file[1] , unpack = True , delimiter = "," , usecols = (0) , dtype = str)
		all_ABC_particle_measurements += np.loadtxt(file[1] , unpack = True , delimiter = "," , usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14))

all_ABC_particle_measurements_transpose = np.array(all_ABC_particle_measurements).T






# Loop over all patient IDs and summarise ABC results
for patient in np.unique(all_patient_IDs):


	###-----------------------------------------------------------------------------
	### Pre-process raw ABC particle data
	###-----------------------------------------------------------------------------

	# Define and create  path for output data
	outDir = "./ABC_processed_results/" + patient + "/"
	Path(outDir).mkdir(parents=True, exist_ok=True)


	# Find indices in list which match patient 
	index_filter = [i for i, ID in enumerate(all_patient_IDs) if (ID == patient)]


	# Filter other data on these indices
	all_ABC_particle_measurements_patient_specific = all_ABC_particle_measurements_transpose[index_filter]


	# Unpack all_ABC_particle_measurements array into separate lists
	s_unfiltered = all_ABC_particle_measurements_patient_specific.T[0]
	k_unfiltered = all_ABC_particle_measurements_patient_specific.T[1]
	q_unfiltered = all_ABC_particle_measurements_patient_specific.T[2]
	seed_unfiltered = all_ABC_particle_measurements_patient_specific.T[3]
	core_distances = all_ABC_particle_measurements_patient_specific.T[4]
	core_sample_radius = all_ABC_particle_measurements_patient_specific.T[5]
	margin_distances = all_ABC_particle_measurements_patient_specific.T[6]
	margin_sample_x = all_ABC_particle_measurements_patient_specific.T[7]
	margin_sample_y = all_ABC_particle_measurements_patient_specific.T[8]
	margin_sample_radius = all_ABC_particle_measurements_patient_specific.T[9]
	leadingEdge_distances = all_ABC_particle_measurements_patient_specific.T[10]
	leadingEdge_sample_x = all_ABC_particle_measurements_patient_specific.T[11]
	leadingEdge_sample_y = all_ABC_particle_measurements_patient_specific.T[12]
	leadingEdge_sample_radius = all_ABC_particle_measurements_patient_specific.T[13]


	# Parse some of these to ints
	k_unfiltered = [int(val) for val in k_unfiltered]
	q_unfiltered = [int(val) for val in q_unfiltered]
	seed_unfiltered = [int(val) for val in seed_unfiltered]
	margin_sample_x = [int(val) for val in margin_sample_x]
	margin_sample_y = [int(val) for val in margin_sample_y]
	if (np.isnan(leadingEdge_sample_x[0]) == False):
		leadingEdge_sample_x = [int(val) for val in leadingEdge_sample_x]
		leadingEdge_sample_y = [int(val) for val in leadingEdge_sample_y]


	# Round the distances to 2 decimal places
	core_distances = [round(val , 2) for val in core_distances]
	margin_distances = [round(val , 2) for val in margin_distances]
	if (np.isnan(leadingEdge_sample_x[0]) == False):
		leadingEdge_distances = [round(val , 2) for val in leadingEdge_distances]


	# Sort particle distances for core, margin and leading edge for use later on
	core_distances_ascending = np.sort(core_distances)
	margin_distances_ascending = np.sort(margin_distances)
	leadingEdge_distances_ascending = np.sort(leadingEdge_distances)


	# Do the same for the summation of core, margin (and leading edge) particle distances
	if (np.isnan(leadingEdge_sample_x[0]) == False):
		combined_distances = [a + b for a, b in zip(core_distances , margin_distances)]
	else:
		combined_distances = [a + b + c for a, b, c in zip(core_distances , margin_distances , leadingEdge_distances)]






		


	###-----------------------------------------------------------------------------
	### Filter particle data based on cutoff values 
	###-----------------------------------------------------------------------------

	# Repeat until size of posterior sample set is around the size set by user
	posterior_size = round(args.posterior_size*len(s_unfiltered)/100.0)
	attempts = 0
	while (True):

		attempts += 1

		s, k, q = [], [], []

		# Find absolute value of cutoff for core, margin and leading edge distances based on the input epsilon values
		core_cutoff = core_distances_ascending[int(args.epsilon_core*len(core_distances_ascending) - 1)]
		margin_cutoff = margin_distances_ascending[int(args.epsilon_margin*len(margin_distances_ascending) - 1)]
		leadingEdge_cutoff = leadingEdge_distances_ascending[int(args.epsilon_leadingEdge*len(leadingEdge_distances_ascending) - 1)]


		# Construct a posterior sample set consisting of particles which are closer than the region cutoff values
		if (np.isnan(leadingEdge_sample_x[0]) == False):
			for i in range(len(s_unfiltered)):
				# Check if distance to experimental data is below input thresholds
				if (core_distances[i] <= core_cutoff) and (margin_distances[i] <= margin_cutoff) and (leadingEdge_distances[i] <= leadingEdge_cutoff):
					s.append(s_unfiltered[i])
					k.append(k_unfiltered[i])
					q.append(q_unfiltered[i])

		else:
			for i in range(len(s_unfiltered)):
				# Check if distance to experimental data is below input thresholds
				if (core_distances[i] <= core_cutoff) and (margin_distances[i] <= margin_cutoff):
					s.append(s_unfiltered[i])
					k.append(k_unfiltered[i])
					q.append(q_unfiltered[i])




		# If previously constructed posterior sample set is not close enough to the specificied desired size (args.posterior_size), then alter the epsilon values accordingly and repeat above step
		if (len(s) > round((args.posterior_size+2)*len(s_unfiltered)/100.0)):
			if (args.epsilon_core < 0.001):
				if (args.epsilon_core < args.epsilon_margin):
					args.epsilon_margin = np.max([0.0001 , args.epsilon_margin-0.0002])
				else:
					args.epsilon_core = np.max([0.0001 , args.epsilon_core-0.0002])
			else:
				args.epsilon_core = np.max([0.0001 , args.epsilon_core-0.0002])

		elif (len(s) < round(np.max([0 , (args.posterior_size-2)])*len(s_unfiltered)/100.0)):
			if (args.epsilon_margin > 0.2):
				if (args.epsilon_margin > args.epsilon_core):
					args.epsilon_core = np.min([1.0 , args.epsilon_core+0.0015])
				else:
					args.epsilon_margin = np.min([1.0 , args.epsilon_margin+0.0015])
			else:
				args.epsilon_margin = np.min([1.0 , args.epsilon_margin+0.0015])
		else:
			break




		# If algorithm gets stuck (because epsilon increments too large) then choose posterior based on combined core + margin distance
		if (attempts >= 1e4):
			print("\tReverted to sorting particles by combined epsilon")

			allSamples = [(s, k, q, dist) for s, k, q, dist in zip(s_unfiltered , k_unfiltered , q_unfiltered , combined_distances)]
			allSamples_sorted = sorted(allSamples, key = lambda tup: tup[3])
		
			# Take the closest 1000 particles
			s = [line[0] for line in allSamples_sorted[:args.posterior_size]]
			k = [line[1] for line in allSamples_sorted[:args.posterior_size]]
			q = [line[2] for line in allSamples_sorted[:args.posterior_size]]
			break



	# Posterior sample set
	combined_full_posterior = [(a , b , c) for a, b, c in zip(s , k , q)]






	###-----------------------------------------------------------------------------
	### Compute point estimators and errors from newly constructed posterior sample set
	###-----------------------------------------------------------------------------


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




		# Initial ecDNA copy number, k
		k_choices = np.arange(k_range_low , k_range_high+1 , 1)

		# ecDNA selection strength, s
		s_choices = np.arange(round(s_range_low , 1) , round(s_range_high , 1)+0.1 , 0.1)
		s_choices = [round(val , 1) for val in s_choices]
		
		# Cell pushing strength, q
		q_choices = q_values






	### Bootstrapping errors for estimates on s and n
	numBootstraps = 1000
	bootstrapped_s_k_q = []

	for i in range(numBootstraps):

		# Bootstrap sample (sample size set to 20% of posterior sample set size)
		subSample_indices = np.random.choice( [int(val) for val in np.linspace(0 , len(s)-1 , len(s))], size = round(0.2*len(s)), replace = True )

		posterior_subSample = []
		for ind in subSample_indices:
			posterior_subSample.append(combined_full_posterior[ind])


		# Predict s and n from this subsample 
		bootstrapped_s_k_q.append([np.median([line[0] for line in posterior_subSample]) , np.median([line[1] for line in posterior_subSample]) , np.median([line[2] for line in posterior_subSample])])




	### Compute point estimator and errors
	inferred_params = [stats.mode([line[0] for line in bootstrapped_s_k_q])[0] , stats.mode([line[1] for line in bootstrapped_s_k_q])[0] , stats.mode([line[2] for line in bootstrapped_s_k_q])[0]]


	# Errors on estimate for s
	s_error = np.std([line[0] for line in combined_full_posterior])
	s_CI_lower = np.max([inferred_params[0] - s_error , s_choices[0]])
	s_CI_upper = np.min([inferred_params[0] + s_error , s_choices[-1]])

	# Errors on estimate for k
	k_error = np.std([line[1] for line in combined_full_posterior])
	k_CI_lower = inferred_params[1] - k_error
	k_CI_upper = inferred_params[1] + k_error

	# Errors on estimate for q
	q_error = round(np.std([q_mapping[line[2]] for line in combined_full_posterior]))
	q_CI_lower = list(q_mapping.keys())[list(q_mapping.values()).index(np.max([q_mapping[inferred_params[2]] - q_error , 1]))]
	q_CI_upper = list(q_mapping.keys())[list(q_mapping.values()).index(np.min([q_mapping[inferred_params[2]] + q_error , 6]))]


	inferred_params_df = pd.DataFrame([inferred_params[0] , inferred_params[1] , q_mapping[inferred_params[2]]]).transpose()
	inferred_params_df.columns = ['s' , 'k' , 'q']











	###-----------------------------------------------------------------------------
	### Find posterior particle with best-fit parameters which is clostest to patient data
	###-----------------------------------------------------------------------------

	#combined_distances = [a + b for a, b in zip(core_distances , margin_distances)]
	all_params_and_probability = [[a, b, c, d, e, f, g, h, i, j, k, l] for a, b, c, d, e, f, g, h, i, j, k, l in zip(s_unfiltered, k_unfiltered, q_unfiltered, seed_unfiltered, combined_distances, core_sample_radius, margin_sample_x, margin_sample_y, margin_sample_radius, leadingEdge_sample_x, leadingEdge_sample_y, leadingEdge_sample_radius )]
	all_params_and_probability_best_fit_combination_only = [combo for combo in all_params_and_probability if (combo[0] == inferred_params[0]) and (combo[1] == inferred_params[1]) and (combo[2] == inferred_params[2])]
	all_params_and_probability_best_fit_combination_only_sorted = sorted(all_params_and_probability_best_fit_combination_only, key = lambda tup: tup[4])

	with open(outDir + "closest_simulations.csv" , 'w') as f:
		
		f.write("s,k,q,seed,combined_distance,core_radius,margin_x,margin_y,margin_radius,leadingEdge_x,leadingEdge_y,leadingEdge_radius\n")
		for line in all_params_and_probability_best_fit_combination_only_sorted[:3]:
			for element in line[:-1]:
				f.write("{},".format(element))
			f.write("{}\n".format(line[-1]))





	# Write inferred parameters to file 
	with open(outDir + "inferred_params.csv" , 'w') as f:
		f.write(",s,k,q\n")
		f.write("Point_estimator,{},{},{}\n".format(inferred_params[0] , round(inferred_params[1]) , inferred_params[2]))
		f.write("errorLowerBound,{},{},{}\n".format(round(s_CI_lower , 1) , round(k_CI_lower) , q_CI_lower))
		f.write("errorUpperBound,{},{},{}\n".format(round(s_CI_upper , 1) , round(k_CI_upper) , q_CI_upper))












	###-----------------------------------------------------------------------------
	### Make plot of point estimates
	###-----------------------------------------------------------------------------

	s_choices = [round(val , 1) for val in s_choices]
	s_bar_width = s_choices[1]-s_choices[0]
	s_posterior_normalised = [list(s).count(val)/float(len(s)) for val in s_choices]


	k_bar_width = k_choices[1]-k_choices[0]
	k_posterior_normalised = [list(k).count(val)/float(len(k)) for val in k_choices]


	q_bar_width = q_choices[1]-q_choices[0]
	q_posterior_normalised = [list(q).count(val)/float(len(q)) for val in q_choices]


	fig , axes = plt.subplot_mosaic("ABC" , width_ratios = [1,1,1] , figsize = (8 , 4) )

	#print(inferred_params[0],s_CI_lower,s_CI_upper)

	axes["A"].scatter(1 , inferred_params[0] , color = "black" , marker = 's' , s = 50 , zorder = 0)
	axes["A"].errorbar(1 , inferred_params[0], yerr = np.array([[abs(inferred_params[0] - s_CI_lower) , s_CI_upper - inferred_params[0]]]).T , xerr = None , ecolor = "red", linestyle = "dashed" , elinewidth = 2 , capthick = 2 , capsize = 5 , zorder = -10)

	axes["B"].scatter(1 , inferred_params[1] , color = "black" , s = 50 , marker = 's')
	axes["B"].errorbar(1 , inferred_params[1] , yerr = np.array([[np.min([inferred_params[1] - k_CI_lower , inferred_params[1]]), k_CI_upper - inferred_params[1]]]).T , xerr = None , ecolor = "red", elinewidth = 2 , capthick = 2 , capsize = 5 , zorder = -10)

	axes["C"].scatter(1 , q_mapping[inferred_params[2]] , color = "black" , s = 50 , marker = 's')
	axes["C"].errorbar(1 , q_mapping[inferred_params[2]], yerr = np.array([[np.min([q_mapping[inferred_params[2]] - q_mapping[q_CI_lower] , q_mapping[inferred_params[2]]]), q_mapping[q_CI_upper] - q_mapping[inferred_params[2]]]]).T , xerr = None , ecolor = "red", elinewidth = 2 , capthick = 2 , capsize = 5 , zorder = -10)



	axes["A"].set_ylim([s_choices[0] - 0.1 , s_choices[-1] + 0.3])
	axes["B"].set_ylim([k_choices[0] - 5 , k_choices[-1] + 5])
	axes["C"].set_ylim([0 , len(q_choices) + 1])


	axes["A"].tick_params(bottom = False , width = 1.5 , labelsize = 13)
	axes["A"].set_xticks([])

	axes["B"].tick_params(bottom = False , width = 1.5 , labelsize = 13)
	axes["B"].set_xticks([])

	axes["C"].tick_params(bottom = False , width = 1.5 , labelsize = 13)
	axes["C"].set_xticks([])



	axes["A"].set_ylabel(r"$s$" , fontsize = 18 )
	axes["B"].set_ylabel(r"$k$" , fontsize = 18 )
	axes["C"].set_ylabel(r"$q$" , fontsize = 18 )

	axes["C"].set_yticks(np.linspace(1 , len(q_choices) , len(q_choices)))
	axes["C"].set_yticklabels([str(val) for val in q_choices])


	[x.set_linewidth(1.5) for x in axes["A"].spines.values()]
	[x.set_linewidth(1.5) for x in axes["B"].spines.values()]
	[x.set_linewidth(1.5) for x in axes["C"].spines.values()]


	plt.tight_layout()
	plt.subplots_adjust(wspace = 0.75)
	plt.savefig(outDir + "inferred_params.png" , dpi = 600 , format = 'png')

	plt.close()















	###-----------------------------------------------------------------------------
	### Plot posterior distributions
	###-----------------------------------------------------------------------------


	plt.clf()
	fig , axes = plt.subplot_mosaic("ABC;DEF;GHI" , width_ratios = [1,1,1] , figsize = (12 , 12) )


	# Diagonal panels
	axes["A"].bar(k_choices , k_posterior_normalised , color = "#D392BF" , edgecolor = None , align = 'center' , width = k_bar_width)

	axes["E"].bar(s_choices , s_posterior_normalised , color = "#D392BF" , edgecolor = None , align = 'center' , width = s_bar_width)

	axes["I"].bar(np.linspace(1 , len(q_choices) , len(q_choices)) , q_posterior_normalised , color = "#D392BF" , edgecolor = None , align = 'center' , width = q_bar_width)
	axes["I"].set_xticks(np.linspace(1 , len(q_choices) , len(q_choices)))
	axes["I"].set_xticklabels([str(val) for val in q_choices])





	# Draw the 90% credible regions for s and n
	# First, take not of the current x and y limits for the bar plot. Plotting the CI will cause these to change, but we can reset these limits after to their current value
	s_plot_xlim = axes["E"].get_xlim()
	s_plot_ylim = axes["E"].get_ylim()

	k_plot_xlim = axes["A"].get_xlim()
	k_plot_ylim = axes["A"].get_ylim()

	q_plot_xlim = axes["I"].get_xlim()
	q_plot_ylim = axes["I"].get_ylim()



	# Plot 90% credible region for s
	axes["E"].plot([s_CI_lower , s_CI_lower] , [0 , 1.1] , color = "red" , linestyle = 'dashed')
	axes["E"].plot([s_CI_upper , s_CI_upper] , [0 , 1.1] , color = "red" , linestyle = 'dashed')
	axes["E"].axvspan(s_CI_lower, s_CI_upper, alpha = 0.25, color = "#f5b105")

	axes["E"].set_xlim(s_plot_xlim)
	axes["E"].set_ylim(s_plot_ylim)


	# Plot 90% credible region for n
	axes["A"].plot([k_CI_lower , k_CI_lower] , [0 , 1.1] , color = "red" , linestyle = 'dashed')
	axes["A"].plot([k_CI_upper , k_CI_upper] , [0 , 1.1] , color = "red" , linestyle = 'dashed')
	axes["A"].axvspan(k_CI_lower, k_CI_upper, alpha = 0.25, color = "#f5b105")

	axes["A"].set_xlim(k_plot_xlim)
	axes["A"].set_ylim(k_plot_ylim)


	# Plot credible region for q
	axes["I"].plot([q_mapping[q_CI_lower]-(0.5*q_bar_width) , q_mapping[q_CI_lower]-(0.5*q_bar_width)] , [0 , 1.1] , color = "red" , linestyle = 'dashed')
	axes["I"].plot([q_mapping[q_CI_upper]+(0.5*q_bar_width) , q_mapping[q_CI_upper]+(0.5*q_bar_width)] , [0 , 1.1] , color = "red" , linestyle = 'dashed')
	axes["I"].axvspan(q_mapping[q_CI_lower]-(0.5*q_bar_width), q_mapping[q_CI_upper]+(0.5*q_bar_width), alpha = 0.25, color = "#f5b105")

	axes["I"].set_xlim(q_plot_xlim)
	axes["I"].set_ylim(q_plot_ylim)




	# Lower off-diagonal panels
	all_posterior_samples_df = pd.DataFrame([s , k , [q_mapping[val] for val in q]]).transpose()
	all_posterior_samples_df.columns = ['s' , 'k' , 'q']

	sns.kdeplot(data = all_posterior_samples_df, x = "k", y = "s" , color = "#E8CADD" , thresh = 0.2 , fill = True , ax = axes["D"] , warn_singular = False) 
	sns.kdeplot(data = all_posterior_samples_df, x = "k", y = "q" , bw_adjust = 2 , thresh = 0.2 , color = "#E8CADD" , fill = True , ax = axes["G"] , warn_singular = False) 
	sns.kdeplot(data = all_posterior_samples_df, x = "s", y = "q" , bw_adjust = 2 , thresh = 0.2 , color = "#E8CADD" , fill = True , ax = axes["H"] , warn_singular = False) 

	sns.scatterplot(data = inferred_params_df , x = "k" , y = "s" , marker = "X" , color = "#FFCC29" , linewidth = 1.5 , edgecolor = "black" , s = 250 , ax = axes["D"] , label = "Inferred value")
	sns.scatterplot(data = inferred_params_df , x = "k" , y = "q" , marker = "X" , color = "#FFCC29" , linewidth = 1.5 , edgecolor = "black" , s = 250 , ax = axes["G"])
	sns.scatterplot(data = inferred_params_df , x = "s" , y = "q" , marker = "X" , color = "#FFCC29" , linewidth = 1.5 , edgecolor = "black" , s = 250 , ax = axes["H"])







	# More formatting
	axes["H"].set_xlabel(r"$s$" , fontsize = 30 , labelpad = 10)
	axes["G"].set_xlabel(r"$k$" , fontsize = 30 , labelpad = 10)
	axes["I"].set_xlabel(r"$q$" , fontsize = 30 , labelpad = 10)

	axes["D"].set_ylabel(r"$s$" , fontsize = 30 , labelpad = 40)
	axes["A"].set_ylabel(r"$k$" , fontsize = 30 , labelpad = 10)
	axes["G"].set_ylabel(r"$q$" , fontsize = 30 , labelpad = 5)

	axes["A"].set_xticks([])
	axes["A"].set_yticks([])
	axes["D"].set_xticks([])
	axes["E"].set_xticks([])
	axes["E"].set_yticks([])
	axes["H"].set_yticks([])
	axes["I"].set_yticks([])


	axes["A"].tick_params(labelsize = 20)
	axes["D"].tick_params(labelsize = 20)
	axes["G"].tick_params(labelsize = 20)
	axes["H"].tick_params(labelsize = 20)
	axes["I"].tick_params(labelsize = 20)

	axes["E"].set_xlim([s_choices[0] - 0.3 , s_choices[-1] + 0.3])
	axes["H"].set_xlim([s_choices[0] - 0.3 , s_choices[-1] + 0.3])
	axes["D"].set_ylim([s_choices[0] - 0.3 , s_choices[-1] + 0.3])

	axes["A"].set_xlim([k_choices[0] - 5 , k_choices[-1] + 15])
	axes["D"].set_xlim([k_choices[0] - 5 , k_choices[-1] + 15])
	axes["G"].set_xlim([k_choices[0] - 5 , k_choices[-1] + 15])

	axes["G"].set_ylim([0 , len(q_choices) + 1])
	axes["H"].set_ylim([0 , len(q_choices) + 1])

	axes["D"].legend(fontsize = 17)


	axes["G"].set_yticks(np.linspace(1 , len(q_choices) , len(q_choices)))
	axes["G"].set_yticklabels([str(val) for val in q_choices])




	axes["D"].set_yticks(np.arange(round(s_choices[0])-1 , round(s_choices[-1])+1 , 1))
	axes["D"].set_yticklabels([str(val) for val in np.arange(round(s_choices[0])-1 , round(s_choices[-1])+1 , 1)])

	axes["H"].set_xticks(np.arange(round(s_choices[0])-1 , round(s_choices[-1])+1 , 1))
	axes["H"].set_xticklabels([str(val) for val in np.arange(round(s_choices[0])-1 , round(s_choices[-1])+1 , 1)])

	axes["G"].set_xticks(np.arange(0 , k_choices[-1]+1 , 50))
	axes["G"].set_xticklabels([str(val) for val in np.arange(0 , k_choices[-1]+1 , 50)])

	# axes["I"].set_xticks([1 , 2 , 3 , 4 , 5 , 6])
	# axes["I"].set_xticklabels(["1" , "2" , "5" , "10" , "50" , "    1000"])



	axes["A"].set(ylabel = None)
	axes["D"].set(xlabel = None)
	axes["H"].set(ylabel = None)



	fig.delaxes(axes["B"])
	fig.delaxes(axes["C"])
	fig.delaxes(axes["F"])



	#for ax in axes:
	#	[x.set_linewidth(1.5) for x in ax.spines.values()]
	#	ax.tick_params(width = 1.5)

	for ax in ["A" , "D" , "E" , "G" , "H" , "I"]:
		[x.set_linewidth(1.5) for x in axes[ax].spines.values()]
		axes[ax].tick_params(width = 1.5)




	axes["A"].text(0.025 , 1.42 , "Patient: " + patient , weight = "bold" , fontsize = 30 , transform = axes["A"].transAxes , horizontalalignment = 'left' , verticalalignment = 'top')
	axes["A"].text(0.025 , 1.24 , "n={}\nAcceptance rate = {}%".format(len(s) , round(100.0*float(len(s))/float(len(s_unfiltered)) , 2)) , fontsize = 20 , transform = axes["A"].transAxes , horizontalalignment = 'left' , verticalalignment = 'top')



	plt.subplots_adjust(wspace = 0.05 , hspace = 0.05)

	plt.savefig(outDir + "posteriorDistributions.png" , dpi = 300 , format = 'png')
	plt.close()

	print("Summarized ABC data for patient " + patient)




