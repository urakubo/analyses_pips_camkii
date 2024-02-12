
import os, glob, pickle, pprint, copy

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1

import utils
import parameters as p
import colormap as c
from scipy import ndimage

plt.rcParams.update(p.rc_param)
	
	
def arrange_graph_bar(ax, panel_dx, y0, panel_size_x, panel_size_y):
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax_h, ax_w = ax.bbox.height, ax.bbox.width
	loc = ax.get_position()
	if y0 is None:
		y0 = loc.y0
	ax.set_position([loc.x0+panel_dx, y0, panel_size_x, panel_size_y])
	
	
	
def make_a_figure( d, target_molecules ):
	
	# Parameters
	vmax       = 0.7

	# Figure specification
	transps = [(0,1,2),(1,2,0),(2,0,1)]
	titles  = [ 'yz', 'xy', 'zx']
	num_columns = 10
	num_rows    = 3
	fig = plt.figure(figsize=(20, 8)) # , tight_layout=False
	
	left, right, bottom, top = 0.0, 0.95, 0.10, 0.99
	wspace, hspace = 0.2, 0.1
	plt.subplots_adjust(left, bottom, right, top, wspace, hspace)
	fig.suptitle( "Time step: {}, MC step: {:e}".format( d['time_frame'], d['mc_step'] ) )
	
	panel_size = 0.15
	
	# Plot profiles
	yloc = []
	for row,(title, transp) in enumerate(zip(titles, transps)):
		axes = []
		column = 1
		ax    = utils.plot_regions_condenstate_from_a_direction_(fig, num_rows, num_columns, row, column, d, transp )
		loc0  = ax.get_position()
		axes.append(ax)
		columns = target_molecules
		axes1 = utils.plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, \
				d, transp, pre_rotated=False, colorbar=False )
		axes.extend(axes1)
		
		yloc.append(loc0.y0)
		for a in axes:
			loc = a.get_position()
			#print('loc ', loc)
			a.set_position([loc.x0, yloc[-1], panel_size, panel_size])
	
	panel_dx = 0.03
	
	column = 7
	
	# Plot RDF
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*2 )
	errorbar= 'shaded' # errorbar='shaded', 'line', or 'final_frame_alone'
	utils.plot_a_rdf( ax, d, errorbar=errorbar, legend=True, target_molecules = target_molecules.keys() )
	arrange_graph_bar(ax, panel_dx, yloc[2], panel_size/2, panel_size)
	
	
	return fig
	
	
if __name__ == '__main__':
	
	# Dataset 1
	# Input files
	dir_edited_data 		= os.path.join('data', 'mix_two_three_components')
	filenames_edited_data = ['CG',\
				'SP',\
				'SPG1',\
				'SPG2']
	targets_molecules 		= [ {'CaMKII':2, 'GluN2B':3} ,\
						{'STG':4,'PSD95':5},\
						{'GluN2B':3, 'STG':4,'PSD95':5},\
						{'GluN2B':3, 'STG':4,'PSD95':5}]
	# Output files
	dir_imgs = os.path.join('imgs', 'mix_two_three_components','profiles')
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	'''
	# Dataset 2
	# Input files
	dir_edited_data 		= os.path.join('data', 'self_affinity')
	filenames_edited_data   = ['CaMKIIalone',\
							'GluN2Balone',\
							'STGalone']
	targets_molecules 		= [ {'CaMKII':2, 'GluN2B':3,'STG':4,'PSD95':5} ] * 3
	# Output files
	dir_imgs = os.path.join('imgs', 'self_affinity','profiles')
	os.makedirs(dir_imgs, exist_ok=True)
	'''
	
	
	for filename, target_molecules in zip(filenames_edited_data, targets_molecules):

		# Load data
		prefix = filename
		sigma = 2
		suffix = 'sigma_{}'.format(sigma)
		print('Target: {}, sigma: {}'.format(filename, sigma))
		d   = utils.load(dir_edited_data, prefix, suffix)
		
		# Make figure
		fig = make_a_figure(d, target_molecules)
		
		# Save figure
		fig.savefig( os.path.join(dir_imgs, '{}_sigma_{}.svg'.format( filename, sigma ) ) )
		fig.savefig( os.path.join(dir_imgs, '{}_sigma_{}.png'.format( filename, sigma ) ) , dpi=150)
		#plt.show()
		plt.clf()
		plt.close(fig=fig)

