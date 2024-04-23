
##
## Plot profile - CaMKII and GluN2B only
## Need d['rdf']['CaMKII_bead']
## 


import os, sys, glob, pickle, pprint, copy

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1

import utils
import parameters as p
import colormap as c
from scipy import ndimage


plt.rcParams.update(p.rc_param)
	
	
def peak(x, c):
    return np.exp(-np.power(x - c, 2) / 16.0)
	
def lin_interp(x, y, i, half):
    return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))
	
def half_max_x(x, y):
    half = max(y)/2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    return lin_interp(x, y, zero_crossings_i[0], half)
	
	
def plot_a_rdf( ax, r, y, legend=True,  ylim = (-0.006,0.66) ):

	color = c.cmap_universal_ratio['CaMKII']
	ax.step( r, y, color=color, label='CaMKII bead')
	if legend==True:
		ax.legend(frameon=False)
	
	ax.set_xlabel('Distance from \n center-of-mass (l.u.)')
	ax.set_ylabel('(beads / voxel)')
	ax.set_xlim(0,40)
	ax.set_ylim(*ylim)
	ax.set_xticks(range(0,50,10))
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	
	return
	
	
def arrange_graph_bar(ax, panel_size_x, panel_size_y):
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax_h, ax_w = ax.bbox.height, ax.bbox.width
	loc = ax.get_position()
	y0 = loc.y0
	ax.set_position([loc.x0, y0, panel_size_x, panel_size_y])
	
	
def plot_bar(ax, data, color, ylim, ylegend, width=0.5):
	ax.bar(*zip(*data.items()), width=width, color=color)
	ax.set_ylabel(ylegend)
	ax.set_ylim(ylim)
	
	
def make_a_figure( d ):
	
	# Figure specification
	transps = [(0,1,2),(1,2,0),(2,0,1)]
	titles  = [ 'yz', 'xy', 'zx']
	num_columns = 10
	num_rows    = 3
	panel_size = 0.15
	left, right, bottom, top = 0.0, 0.95, 0.10, 0.99
	wspace, hspace = 0.2, 0.1
	
	fig = plt.figure(figsize=(20, 8)) # , tight_layout=False
	plt.subplots_adjust(left, bottom, right, top, wspace, hspace)
	fig.suptitle( "Time step: {}, MC step: {:e}".format( d['time_frame'], d['mc_step'] ) )
	
	
	# Plot RDF
	r = d['rdf_bins'][2:-1]
	y = d['rdf']['CaMKII_bead'][2:]
	
	hmx = half_max_x(r, y)
	print('hmx ', hmx)
	
	column = 3
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*1 )
	plot_a_rdf( ax, r, y, legend=True ) # , ylim = (-0.006,0.46)
	arrange_graph_bar(ax, panel_size/2, panel_size)
	ax.set_title( 'r = {:.3f}'.format(hmx) )
	
	return fig, hmx
	
	
if __name__ == '__main__':
	
	
	
	# CG Valency length
	#'''
	subdirs    = ['val{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
	filenames_input  = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_output = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
	dir_input        = 'CG_valency_length'
	dir_target       = 'CG_valency_length'
	#'''
	suffix = 'radius_CaMKII_bead'
	
	# Shared part of initialization
	dir_lammpstrj    = os.path.join('..', 'lammpstrj3', dir_input)
	dir_edited_data  = os.path.join('data3', dir_target)
	dir_imgs = os.path.join('imgs3', dir_target, suffix)
	
	os.makedirs(dir_edited_data, exist_ok=True)
	os.makedirs(dir_imgs, exist_ok=True)
	
	hmxs = {}
	for filename_input, filename_output in zip(filenames_input, filenames_output):
		# Load data
		print("\n"+filename_input)
		sampling_frame = utils.get_num_frames(dir_lammpstrj, filename_input)
		print("The last timeframe was loaded." )
		
		print("sampling_frame ", sampling_frame )
		types, positions_grid_coord,ids_molecule, mc_step = \
			utils.load_lammpstrj( dir_lammpstrj, filename_input, sampling_frame )
		
		
		d = {}
		# RDF, CaMKII binding beads
		rdf, rdf_bins = \
			utils.get_rdf_CaMKII_bead( dir_lammpstrj, filename_input, sampling_frame )
		d['rdf_bins'] = rdf_bins
		d['rdf'] = rdf
		
		# Time info
		d['time_frame'] = sampling_frame
		d['mc_step']    = mc_step
		
		# Make figure
		fig, hmx = make_a_figure(d)
		hmxs[filename_output] = hmx
		# Save figure
		fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename_output ) ) )
		fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename_output ) ) , dpi=150)
		#plt.show()
		plt.clf()
		plt.close(fig=fig)

	# Save the edited data
	prefix = 'Radiuses'
	suffix = 'CaMKII_bead'
	print('{}_{}'.format(prefix, suffix))
	pprint.pprint(hmxs)
	utils.save(dir_edited_data, prefix, suffix, hmxs)

