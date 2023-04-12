

import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt 
import sys
import argparse
import seaborn as sns
import pandas as pd
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Rectangle



# Parse command line argument
parser = argparse.ArgumentParser()
parser.add_argument('--path', help='File path to cells.csv file to be plotted', type=str)
parser.add_argument('--plotCopyNumber', help='Colour cells according to their ecDNA copy number', action='store_true')
args = parser.parse_args()



# Import raw data
x , y , mutations = np.loadtxt(args.path , delimiter="," , unpack=True)


# matplotlib.pyplot.pcolor() function requires the input data to be square. So first find if x or y is the larger dimension.
minX = int(np.min(x))
minY = int(np.min(y))
maxX = int(np.max(x))
maxY = int(np.max(y))

largest_dim = np.max([maxX-minX , maxY-minY])







#~~~~~~~~~~~~~~~~~~~~	Plotting    ~~~~~~~~~~~~~~~~~~~~~#
if (args.plotCopyNumber == False):

	# Initiate input array
	toPlot = np.zeros((int(largest_dim)+1 , int(largest_dim)+1))

	# Populate array
	for i in range(len(x)):
		toPlot[int(x[i])-minX][int(y[i])-minY] = np.sign(mutations[i]) + 1



	# Trim array to remove rows and columns of all zeros (maintaining square array)
	while(True):
		if ((np.all(toPlot[0] == 0)) and (np.all(toPlot[:,0] == 0))):
			toPlot = np.delete(toPlot , obj = 0 , axis = 0)
			toPlot = np.delete(toPlot , obj = 0 , axis = 1)
		else:
			break

	while(True):
		if ((np.all(toPlot[-1] == 0)) and (np.all(toPlot[:,-1] == 0))):
			toPlot = np.delete(toPlot , obj = -1 , axis = 0)
			toPlot = np.delete(toPlot , obj = -1 , axis = 1)
		else:
			break



	toPlot_x = np.linspace(0 , toPlot.shape[0] , toPlot.shape[0])
	toPlot_y = np.linspace(0 , toPlot.shape[0] , toPlot.shape[0])
	X , Y = np.meshgrid(toPlot_x , toPlot_y)


	# Set up colour map
	cmap = mpl.colors.ListedColormap(['white' , (10/256, 22/256, 230/256, 0.7) , [1, 0.34117647, 0, 1]])
	norm = mpl.colors.Normalize(vmin = 0, vmax = 3)

	# Plot formatting stuff
	plt.figure(figsize=( 8 , 8 ))
	plt.gca().axis('off')

	# Plot
	plt.pcolor(X , Y , toPlot , cmap = cmap , shading = 'auto' , snap = True)

	# Save plot
	plt.savefig(args.path + ".ecDNA_carryingStatus.png" , dpi = 300 , format = 'png')



else:
	### First, plot empty lattice points and cells with zero ecDNA
	fig, ax = plt.subplots(nrows = 1 , ncols = 2 , gridspec_kw = {'width_ratios': [5, 1]})

	# Initiate input array
	toPlot = np.zeros((int(largest_dim)+1 , int(largest_dim)+1))

	# Populate array
	for i in range(len(x)):
		toPlot[int(x[i])-minX][int(y[i])-minY] = np.sign(mutations[i]) + 1



	toPlot_x = np.linspace(0 , int(largest_dim)+1 , int(largest_dim)+1)
	toPlot_y = np.linspace(0 , int(largest_dim)+1 , int(largest_dim)+1)
	X , Y = np.meshgrid(toPlot_x , toPlot_y)


	# Set up colour map
	cmap = mpl.colors.ListedColormap(['white' , (10/256,22/256,230/256,0.5) , (1,1,1,0)])
	norm = mpl.colors.Normalize(vmin = 0, vmax = 3)

	# Plot formatting stuff
	#plt.figure(figsize=( 8 , 8 ))
	ax[0].axis('off')

	# Plot
	ax[0].pcolor(X , Y , toPlot , cmap = cmap , shading = 'auto' , snap = True)








	### Second, plot cells with ecDNA
	# Initiate input array
	toPlot = np.zeros((int(largest_dim)+1 , int(largest_dim)+1))


	# Populate array
	for i in range(len(x)):
		if (mutations[i] > 0):
			toPlot[int(x[i])-minX][int(y[i])-minY] = mutations[i] + 1



	toPlot_x = np.linspace(0 , int(largest_dim)+1 , int(largest_dim)+1)
	toPlot_y = np.linspace(0 , int(largest_dim)+1 , int(largest_dim)+1)
	X , Y = np.meshgrid(toPlot_x , toPlot_y)


	# Set up colour map
	#plasma = cm.get_cmap('hsv_r', 256)
	plasma = mpl.colormaps['hsv_r']
	#newcolors = plasma(np.linspace(0.5, 1, 256))
	newcolors = plasma([np.sqrt(val) for val in np.linspace(0.65, 1, 256)])
	transparent = np.array([1, 1, 1, 0])
	newcolors[0, :] = transparent
	newcmp = ListedColormap(newcolors)



	# Plot
	if (len([val for val in mutations if val > 0]) > 0):
		coloured_cells = ax[0].pcolor(X , Y , toPlot , cmap = newcmp , shading = 'auto' , vmin = np.min([val for val in mutations if val > 0]) , vmax = np.max(mutations) , snap = True)

		# Colourbar
		cbar = plt.colorbar(coloured_cells , orientation = 'vertical' , aspect = 12 , fraction = 0.05)
		cbar.ax.tick_params(size = 0)
		cbar.set_ticks([])
		#cax = fig.add_axes([0.5, 0.2, 0.2, 0.5])
		#cb = mpl.colorbar.ColorbarBase(cax , orientation = 'vertical')
		#fig.colorbar(coloured_cells, cax=cax, orientation = 'horizontal')








		### Third, plot colourbar and distribution of ecDNA copy numbers
		mutations_dataframe = pd.DataFrame([val for val in mutations if val > 0] , columns = ["#ecDNA"])
		sns.kdeplot(data = mutations_dataframe , y = "#ecDNA" , ax = ax[1] , cut = 0 , zorder = 10 , color = "black" , fill = True , bw_adjust = 0.75)
		sns.despine(bottom = True, left = True)


		# Axes and tick formatting
		ax[1].tick_params(axis = 'y' , labelsize = 10 , width = 0 , pad = 20)
		ax[1].tick_params(axis = 'x' , which = 'both' , bottom = False , labelbottom = False) 
		ax[1].set_xlabel("")

		custom_yticks = [0]
		i = 1
		while(True):
			if (i*100 <= np.max(mutations)):
				custom_yticks.append(i*100)
				i += 1
			else:
				break

		ax[1].set_yticks(custom_yticks)
		ax[1].set_yticklabels([str(val) for val in custom_yticks])

		# Workaround to make plot look nice
		axisRescaleConstant = 0.22*(np.max(mutations) - np.min([val for val in mutations if val > 0]))
		ax[1].set_ylim([ ax[1].get_ylim()[0]-axisRescaleConstant , ax[1].get_ylim()[1]+axisRescaleConstant ])


		# Save plot
		ax[0].set_aspect(1.0)
		#ax[0].margins(x=0, y=-0.25)
		plt.subplots_adjust(wspace = 0.04)
		plt.savefig(args.path + ".ecDNA_copyNumber.png" , dpi = 300 , format = 'png')





	else:
		# Initiate input array
		toPlot = np.zeros((int(largest_dim)+1 , int(largest_dim)+1))

		# Populate array
		for i in range(len(x)):
			toPlot[int(x[i])-minX][int(y[i])-minY] = np.sign(mutations[i]) + 1



		# Trim array to remove rows and columns of all zeros (maintaining square array)
		while(True):
			if ((np.all(toPlot[0] == 0)) and (np.all(toPlot[:,0] == 0))):
				toPlot = np.delete(toPlot , obj = 0 , axis = 0)
				toPlot = np.delete(toPlot , obj = 0 , axis = 1)
			else:
				break

		while(True):
			if ((np.all(toPlot[-1] == 0)) and (np.all(toPlot[:,-1] == 0))):
				toPlot = np.delete(toPlot , obj = -1 , axis = 0)
				toPlot = np.delete(toPlot , obj = -1 , axis = 1)
			else:
				break



		toPlot_x = np.linspace(0 , toPlot.shape[0] , toPlot.shape[0])
		toPlot_y = np.linspace(0 , toPlot.shape[0] , toPlot.shape[0])
		X , Y = np.meshgrid(toPlot_x , toPlot_y)


		# Set up colour map
		cmap = mpl.colors.ListedColormap(['white' , (10/256, 22/256, 230/256, 0.7)])
		norm = mpl.colors.Normalize(vmin = 0, vmax = 2)

		# Plot formatting stuff
		plt.figure(figsize=( 8 , 8 ))
		plt.gca().axis('off')

		# Plot
		plt.pcolor(X , Y , toPlot , cmap = cmap , shading = 'auto' , snap = True)



		# Save plot
		ax[0].set_aspect(1.0)
		#ax[0].margins(x=0, y=-0.25)
		plt.subplots_adjust(wspace = 0.04)
		plt.savefig(args.path + ".ecDNA_copyNumber.png" , dpi = 300 , format = 'png')









