

import numpy as np 
import argparse
import glob
from scipy import stats




def take_sample(sample_centre_x , sample_centre_y , R):

	x_lower_bound = 0
	y_lower_bound = 0
	x_upper_bound = tissue.shape[0]
	y_upper_bound = tissue.shape[1]


	sample = []
	for dx in range(0 , int(round(R))):

		### Loop through y
		for dy in range(0 , int(round(R))):

			### Compute radial distance of cell
			distance = np.sqrt( np.power( dx , 2 ) + np.power( dy , 2 ) )

			### If cell is within current zone, record copy number
			if (distance <= R):

				if (sample_centre_x - dx < x_lower_bound) or (sample_centre_x + dx >= x_upper_bound) or (sample_centre_y - dy < y_lower_bound) or (sample_centre_y + dy >= y_upper_bound):
					return []

				if (dx == 0) and (dy == 0):
					sample.append([sample_centre_x + dx , sample_centre_y + dy , tissue[sample_centre_x + dx][sample_centre_y + dy]])

				elif (dx == 0):
					sample.append([sample_centre_x + dx , sample_centre_y + dy , tissue[sample_centre_x + dx][sample_centre_y + dy]])
					sample.append([sample_centre_x + dx , sample_centre_y - dy , tissue[sample_centre_x + dx][sample_centre_y - dy]])

				elif (dy == 0):
					sample.append([sample_centre_x + dx , sample_centre_y + dy , tissue[sample_centre_x + dx][sample_centre_y + dy]])
					sample.append([sample_centre_x - dx , sample_centre_y + dy , tissue[sample_centre_x - dx][sample_centre_y + dy]])

				else:
					sample.append([sample_centre_x + dx , sample_centre_y + dy , tissue[sample_centre_x + dx][sample_centre_y + dy]])
					sample.append([sample_centre_x + dx , sample_centre_y - dy , tissue[sample_centre_x + dx][sample_centre_y - dy]])
					sample.append([sample_centre_x - dx , sample_centre_y + dy , tissue[sample_centre_x - dx][sample_centre_y + dy]])
					sample.append([sample_centre_x - dx , sample_centre_y - dy , tissue[sample_centre_x - dx][sample_centre_y - dy]])



	return sample









#========================================================================================================================







# Parse command line argument
parser = argparse.ArgumentParser()
parser.add_argument('--simulation', help = 'File path to cells.csv file to be sampled', type = str , required = True)
parser.add_argument('--patient', help = 'Path to directory containing core and margin samples for patient', type = str , required = True)
args = parser.parse_args()




numMarginSamples = 5
all_patient_measurements = np.zeros(3, dtype = object) # each patient (i.e. array element) in this array contains the following: core_sample; margin_sample
all_patient_measurements[0] = []
all_patient_measurements[1] = []
all_patient_measurements[2] = []



# Find and load patient multi-region samples
LEADING_EDGE_SAMPLE_DETECTED = False
leadingEdge_sample_size = 0
for file in glob.glob(args.patient + "/*.dat"):

	sample_name = file.split("/")[-1]

	if (sample_name == "core.dat"):

		### Trim data to remove dead cells (i.e. cells with -1 copy number) and cells with < 3 ecDNA and round to nearest integer
		exp_data = np.loadtxt(file , unpack = True , delimiter = "," , skiprows = 1)
		all_patient_measurements[0] = exp_data

		### Store number of cells for patient sample 
		core_sample_size = len(exp_data)



	if (sample_name == "margin.dat"):
		#print("Margin sample: " + file)

		### Trim data to remove dead cells (i.e. cells with -1 copy number) and cells with < 3 ecDNA and round to nearest integer
		exp_data = np.loadtxt(file , unpack = True , delimiter = "," , skiprows = 1)
		all_patient_measurements[1] = exp_data

		### Store number of cells for patient sample 
		margin_sample_size = len(exp_data)




	if (sample_name == "leading_edge.dat"):

		LEADING_EDGE_SAMPLE_DETECTED = True

		### Trim data to remove dead cells (i.e. cells with -1 copy number) and cells with < 3 ecDNA and round to nearest integer
		exp_data = np.loadtxt(file , unpack = True , delimiter = "," , skiprows = 1)
		all_patient_measurements[2] = exp_data

		### Store number of cells for patient sample 
		leadingEdge_sample_size = len(exp_data)












### Get simulation parameters from input file path 
params = [word for word in args.simulation.split("/") if "n=" in word][0]
path_split = params.split("_")

s_value = [word for word in path_split if "s=" in word][0]
s_value = s_value.split("=")[-1]

k_value = [word for word in path_split if "k=" in word][0]
k_value = k_value.split("=")[-1]

q_value = [word for word in path_split if "q=" in word][0]
q_value = q_value.split("=")[-1]

seed_value = [word for word in args.simulation.split("/") if "seed" in word][0]
seed_value = seed_value.split("=")[-1]




### Load input tumour data
x , y , copyNumber = np.loadtxt(args.simulation , delimiter = "," , unpack = True , skiprows = 1)



### Group all alive tumour cells 
all_cells = [[a, b, c] for a, b, c in zip(x , y , copyNumber) if (c != -1)]



# matplotlib.pyplot.pcolor() function requires the input data to be square. So first find if x or y is the larger dimension.
minX = int(np.min(x))
minY = int(np.min(y))
maxX = int(np.max(x))
maxY = int(np.max(y))

largest_dim = np.max([maxX-minX , maxY-minY])



### Find system centre (centre of mass)
centre_x = round(np.mean(x))
centre_y = round(np.mean(y))



### Rescale centre_x and centre_y
centre_x -= minX
centre_y -= minY

centre_x = int(centre_x)
centre_y = int(centre_y)




### Parse data to a numpy array 
lengthX = int(maxX - minX + 1)
lengthY = int(maxY - minY + 1)
tissue = np.ones( (lengthX , lengthY) , dtype = int )
tissue = -tissue 	# Set all elements to -1 initially, prior to populating with cell ecDNA copy number data 





### Populate array 
for i in range(len(x)):
	if (copyNumber[i] > -1):
		coordX = int(x[i] - minX)
		coordY = int(y[i] - minY)
		tissue[coordX][coordY] = copyNumber[i]






### Estimate rough radius of tumour according to number of cells
tumour_radius = np.power( (len(all_cells)/np.pi) , (1.0/2.0) )





# Estimate simulated sampled core and margin sample sizes based on UNTRIMMED patient data from corresponding regions
core_sample_radius = np.power( (core_sample_size/np.pi) , (1.0/2.0) )
infiltratingMargin_sample_radius = np.power( (margin_sample_size/np.pi) , (1.0/2.0) )
leadingEdge_sample_radius = np.power( (leadingEdge_sample_size/np.pi) , (1.0/2.0) )









#####################
### Core sample
#####################

### Compute corresponding core sample centre coordinates 
core_sample_centreCoordinates = (int(centre_x) , int(centre_y))


### Take core sample 
SIM_core_sample = take_sample(core_sample_centreCoordinates[0] , core_sample_centreCoordinates[1], core_sample_radius)



# Trimmed core data to remove any cells with <3 ecDNA, as is done with the patient data
SIM_core_sample_copyNumberOnly_trimmed = [cell[-1] for cell in SIM_core_sample if ( cell[-1] >= 3)]


# Append a single cell (with the same copy number) to all distributions, to avoid error with Wasserstein metric if trimmed distributions are empty
np.append(all_patient_measurements[0] , 5)
SIM_core_sample_copyNumberOnly_trimmed.append(5)









#####################
### Infiltrating margin and leading edge samples
#####################

### Take margin samples
for i in range(numMarginSamples):


	### Generate random angles for each margin sample to be taken (number of which should match the number of margin samples in patiet)
	margin_sample_theta = np.random.rand()*2.0*np.pi
	leadingEdge_sample_theta = np.random.rand()*2.0*np.pi



	### Compute corresponding margin sample centre coordinates for these thetas (take infiltrating margin sample at 75% or tumour radius)
	margin_sample_centreCoordinates = (int(centre_x + ((tumour_radius*0.75) * np.sin(margin_sample_theta))) , int(centre_y + ((tumour_radius*0.75) * np.cos(margin_sample_theta))))
	leadingEdge_sample_centreCoordinates = (int(centre_x + ((tumour_radius*0.90) * np.sin(leadingEdge_sample_theta))) , int(centre_y + ((tumour_radius*0.90) * np.cos(leadingEdge_sample_theta))))



	### Take margin sample 
	SIM_infiltratingMargin_sample = take_sample(margin_sample_centreCoordinates[0] , margin_sample_centreCoordinates[1] , infiltratingMargin_sample_radius)
	SIM_infiltratingMargin_sample_copyNumberOnly_trimmed = [cell[-1] for cell in SIM_infiltratingMargin_sample if (cell[-1] >= 3)]



	### Take leading edge sample 
	SIM_leadingEdge_sample = take_sample(leadingEdge_sample_centreCoordinates[0] , leadingEdge_sample_centreCoordinates[1] , leadingEdge_sample_radius)
	SIM_leadingEdge_sample_copyNumberOnly_trimmed = [cell[-1] for cell in SIM_leadingEdge_sample if (cell[-1] >= 3)]



	# Append a single cell (with the same copy number) to all distributions, to avoid error with Wasserstein metric if distributions are empty		
	np.append(all_patient_measurements[1] , 5)
	SIM_infiltratingMargin_sample_copyNumberOnly_trimmed.append(5)


	# Append a single cell (with the same copy number) to all distributions, to avoid error with Wasserstein metric if distributions are empty
	if (LEADING_EDGE_SAMPLE_DETECTED):
		np.append(all_patient_measurements[2] , 5)
	SIM_leadingEdge_sample_copyNumberOnly_trimmed.append(5)



	# Compute distance between patient and simulated multi-region samples using Wasserstein metric
	core_distance = stats.wasserstein_distance(all_patient_measurements[0] , SIM_core_sample_copyNumberOnly_trimmed)
	infiltratingMargin_distance = stats.wasserstein_distance(all_patient_measurements[1] , SIM_infiltratingMargin_sample_copyNumberOnly_trimmed)
	if (LEADING_EDGE_SAMPLE_DETECTED):
		leadingEdge_distance = stats.wasserstein_distance(all_patient_measurements[2] , SIM_leadingEdge_sample_copyNumberOnly_trimmed)


	# Print data to console
	patient_name = args.patient.split("/")[-2] if (args.patient[-1] == "/") else args.patient.split("/")[-1]
	if (LEADING_EDGE_SAMPLE_DETECTED):
		print(patient_name + "," + s_value + "," + k_value + "," + q_value + "," + seed_value + ",{},{},{},{},{},{},{},{},{},{}".format(round(core_distance , 5) , round(core_sample_radius , 5) , round(infiltratingMargin_distance , 5) , margin_sample_centreCoordinates[0] + minX , margin_sample_centreCoordinates[1] + minY , round(infiltratingMargin_sample_radius , 5) , round(leadingEdge_distance , 5) , leadingEdge_sample_centreCoordinates[0] + minX , leadingEdge_sample_centreCoordinates[1] + minY , round(leadingEdge_sample_radius , 5)))
	else:
		print(patient_name + "," + s_value + "," + k_value + "," + q_value + "," + seed_value + ",{},{},{},{},{},{},nan,nan,nan,nan".format(round(core_distance , 5) , round(core_sample_radius , 5) , round(infiltratingMargin_distance , 5) , margin_sample_centreCoordinates[0] + minX , margin_sample_centreCoordinates[1] + minY , round(infiltratingMargin_sample_radius , 5)))






























