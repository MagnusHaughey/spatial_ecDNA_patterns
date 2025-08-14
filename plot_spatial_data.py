

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
import matplotlib.patches as patches



# Parse command line argument
parser = argparse.ArgumentParser()
parser.add_argument('--path', help = 'File path to cells.csv file to be plotted', type = str)
parser.add_argument('--plotCopyNumber',  help = 'Colour cells according to their ecDNA copy number', action = 'store_true')
args = parser.parse_args()



# Import raw data
x , y , mutations = np.loadtxt(args.path , delimiter = "," , unpack = True , skiprows = 1)


# matplotlib.pyplot.pcolor() function requires the input data to be square. So first find if x or y is the larger dimension.
minX = int(np.min(x))
minY = int(np.min(y))
maxX = int(np.max(x))
maxY = int(np.max(y))

largest_dim = np.max([maxX-minX , maxY-minY])







#~~~~~~~~~~~~~~~~~~~~	Plotting    ~~~~~~~~~~~~~~~~~~~~~#
if (args.plotCopyNumber == False) or (len([val for val in mutations if val > 0]) == 0):	# If plotCopyNumber flag not specified, or specified but there are no ecDNA+ cells, simply plot ecDNA- cells as blue and ecDNA+ cells as red

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


	# Set up colour map
	cmap = mpl.colors.ListedColormap(['white' , (10/256, 22/256, 230/256, 0.7) , [1, 0.34117647, 0, 1]])
	norm = mpl.colors.Normalize(vmin = 0, vmax = 3)


	# Plot formatting
	plt.figure(figsize = ( 8 , 8 ))
	plt.gca().axis('off')


	# Plot
	toPlot_x = np.linspace(0 , toPlot.shape[0] , toPlot.shape[0])
	toPlot_y = np.linspace(0 , toPlot.shape[0] , toPlot.shape[0])
	X , Y = np.meshgrid(toPlot_x , toPlot_y)

	plt.pcolor(X , Y , toPlot , cmap = cmap , shading = 'auto' , snap = True)


	# Save plot
	if (args.plotCopyNumber == False):
		plt.savefig(args.path + ".ecDNA_carryingStatus.png" , dpi = 300 , format = 'png')
	else:
		plt.savefig(args.path + ".ecDNA_copyNumber.png" , dpi = 300 , format = 'png')




else:	## If plotCopyNumber flag is specified, plot ecDNA- cells as blue, and ecDNA+ cells with colour gradient indicating ecDNA copy number


	### First, plot empty lattice points and cells with zero ecDNA
	fig, ax = plt.subplot_mosaic(mosaic = "AABC" , figsize = (8,5) , width_ratios = (3,3,1,3))

	# Initiate input array
	toPlot = np.zeros((int(largest_dim)+1 , int(largest_dim)+1))

	# Populate array
	for i in range(len(x)):
		toPlot[int(x[i])-minX][int(y[i])-minY] = np.sign(mutations[i]) + 1



	# Set up colour map
	cmap = mpl.colors.ListedColormap(['white' , (10/256,22/256,230/256,0.5) , (1,1,1,0)])
	norm = mpl.colors.Normalize(vmin = 0, vmax = 3)

	# Plot formatting stuff
	ax["A"].axis('off')

	# Plot
	toPlot_x = np.linspace(0 , int(largest_dim)+1 , int(largest_dim)+1)
	toPlot_y = np.linspace(0 , int(largest_dim)+1 , int(largest_dim)+1)
	X , Y = np.meshgrid(toPlot_x , toPlot_y)

	ax["A"].pcolor(X , Y , toPlot , cmap = cmap , shading = 'auto' , snap = True)






	### Second, plot cells with ecDNA
	# Initiate input array
	toPlot = np.zeros((int(largest_dim)+1 , int(largest_dim)+1))


	# Populate array
	for i in range(len(x)):
		if (mutations[i] > 0):
			toPlot[int(x[i])-minX][int(y[i])-minY] = mutations[i] + 1




	# Set up colour map
	plasma = mpl.colormaps['hsv_r']
	newcolors = plasma([np.sqrt(val) for val in np.linspace(0.65, 1, 256)])
	transparent = np.array([1, 1, 1, 0])
	newcolors[0, :] = transparent
	newcmp = ListedColormap(newcolors)




	# Plot
	toPlot_x = np.linspace(0 , int(largest_dim)+1 , int(largest_dim)+1)
	toPlot_y = np.linspace(0 , int(largest_dim)+1 , int(largest_dim)+1)
	X , Y = np.meshgrid(toPlot_x , toPlot_y)
	coloured_cells = ax["A"].pcolor(X , Y , toPlot , cmap = newcmp , shading = 'auto' , vmin = np.min([val for val in mutations if val > 0]) , vmax = np.max(mutations) , snap = True)



	# Colourbar
	#cbar = plt.colorbar(coloured_cells , orientation = 'vertical' , aspect = 12 , fraction = 0.05)
	#cbar.ax.tick_params(size = 0)
	#cbar.set_ticks([])



	### Third, plot colourbar and distribution of ecDNA copy numbers
	mutations_dataframe = pd.DataFrame([val for val in mutations if val > 0] , columns = ["#ecDNA"])
	sns.kdeplot(data = mutations_dataframe , y = "#ecDNA" , ax = ax["C"] , cut = 0 , zorder = 10 , color = "black" , fill = True , bw_adjust = 0.75)
	sns.despine(bottom = True, left = True)


	# Axes and tick formatting
	ax["C"].tick_params(axis = 'y' , labelsize = 10 , width = 0 , pad = 20)
	ax["C"].tick_params(axis = 'x' , which = 'both' , bottom = False , labelbottom = False) 
	ax["C"].set_xlabel("")

	
	if (np.max(mutations) > 100):
		custom_yticks = [0]
		i = 1
		while(True):
			if (i*100 <= np.max(mutations)):
				custom_yticks.append(i*100)
				i += 1
			else:
				break


	else:
		custom_yticks = [0]
		i = 1
		while(True):
			if (i*10 <= np.max(mutations)):
				custom_yticks.append(i*10)
				i += 1
			else:
				break

	ax["C"].set_yticks(custom_yticks)
	ax["C"].set_yticklabels([str(val) for val in custom_yticks])



	# Workaround to make plot look nice
	axisRescaleConstant = 0.22*(np.max(mutations) - np.min([val for val in mutations if val > 0]))
	ax["C"].set_ylim([ ax["C"].get_ylim()[0]-axisRescaleConstant , ax["C"].get_ylim()[1]+axisRescaleConstant ])



	# Manually plot color bar legend in RHS plot
	manual_cbar_x = -0.175
	manual_cbar_y = 0.19
	manual_cbar_width = 0.1
	manual_cbar_height = 0.63

	# Plot colorbar outline
	rect = patches.Rectangle((manual_cbar_x , manual_cbar_y), manual_cbar_width, manual_cbar_height, linewidth = 1 , edgecolor = 'black' , facecolor = 'none' , clip_on = False , transform = ax["C"].transAxes)
	ax["C"].add_patch(rect)

	# Manually plot color bar interior
	for i, color in enumerate(newcolors):
		ax["C"].add_patch(plt.Rectangle((manual_cbar_x, manual_cbar_y + manual_cbar_height/len(newcolors) * i) , manual_cbar_width , manual_cbar_height/len(newcolors) , color = color, linewidth = 0, zorder = 0 , clip_on = False , transform = ax["C"].transAxes))



	# Add extra box below color bar to indicate color of cells with 0 ecDNA
	rect = patches.Rectangle((manual_cbar_x , manual_cbar_y-0.02), manual_cbar_width, 0.013, linewidth = 1 , edgecolor = 'black' , facecolor = (10/256,22/256,230/256,0.5) , clip_on = False , transform = ax["C"].transAxes)
	ax["C"].add_patch(rect)


	# Remove middle axes
	ax["B"].set_axis_off()


	# Save plot
	ax["A"].set_aspect(1.0)
	plt.subplots_adjust(wspace = 0.2)
	plt.savefig(args.path + ".ecDNA_copyNumber.png" , dpi = 300 , format = 'png')






