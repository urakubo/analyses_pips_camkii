##
##
##


import os, sys, glob, pickle, pprint, copy

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1

import utils
import parameters as p
import colormap as c
from scipy import ndimage

import networkx as nx


plt.rcParams.update(p.rc_param)
		
	
	
def arrange_graph_bar(ax, panel_dx, y0, panel_size_x, panel_size_y):
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax_h, ax_w = ax.bbox.height, ax.bbox.width
	loc = ax.get_position()
	if y0 is None:
		y0 = loc.y0
	ax.set_position([loc.x0+panel_dx, y0, panel_size_x, panel_size_y])
	
	
def make_a_figure( d ):
	
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
	#fig.suptitle( "Time step: {}, MC step: {:e}".format( d['time_frame'], d['mc_step'] ) )
	
	panel_size = 0.15
	
	# Plot profiles
	yloc = []
	for row,(title, transp) in enumerate(zip(titles, transps)):
		columns = {'CaMKII':1, 'GluN2B':2,'STG':3,'PSD95':4, 'Shared PSD95':5, 'Unshared PSD95':6}
		axes1 = utils.plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp, pre_rotated=False )
		axes = axes1
		loc0  = axes[0].get_position()
		
		yloc.append(loc0.y0)
		for a in axes:
			loc = a.get_position()
			#print('loc ', loc)
			a.set_position([loc.x0, yloc[-1], panel_size, panel_size])
	
	panel_dx = 0.03
	
	
	#'''
	# Plot RDF
	column = 8
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*2 )
	utils.plot_a_rdf_PSD95( ax, d, legend=True ) # , ylim = (-0.006,0.46)
	arrange_graph_bar(ax, panel_dx, yloc[2], panel_size/2, panel_size)
	#'''
	
	return fig
	
	
	
if __name__ == '__main__':
	
	
	
	# Conc dependnece
	#'''
	filenames    = [str(i).zfill(3) for i in range(48) ]
	dir_target     = 'conc_dependence'
	#'''
	
	# Shared init
	dir_edited_data	= os.path.join('data3', dir_target)
	dir_imgs = os.path.join('imgs3', dir_target,'profiles_shared_PSD95')
	sigma = 2
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	for filename in filenames:
		
		'''
		# Load graph
		d = utils.load(dir_edited_data, filename_edited, suffix)
		multi_graph = d['multi_graph']
		# Find ids_bead for the PSD95 shared by STG and GluN2B
		ids_molecule_shared_PSD, ids_bead_shared_PSD = get_ids_PSD95_shared_by_STG_GluN2B(multi_graph)
		ids_molecule_unshared_PSD, ids_bead_unshared_PSD = get_ids_PSD95_unshared_by_STG_GluN2B(multi_graph)
		
		
		# Load lammpstrj data
		print("\n"+filename_lammpstrj)
		sampling_frame = utils.get_num_frames(dir_lammpstrj, filename_lammpstrj)
		
		types, positions_grid_coord,ids_molecule, mc_step = \
			utils.load_lammpstrj( dir_lammpstrj, filename_lammpstrj, sampling_frame )
		
		# Centering
		center    = utils.get_center_of_mass(types, positions_grid_coord)
		positions_real_coord = utils.centering(positions_grid_coord, center)
		
		# Get concs and condensates
		d = utils.get_concs_and_condensates(types, positions_real_coord, ids_molecule, sigma = 2)
		
		
		# Get the condenate regions of targs_molecule. 'Shared PSD95'
		for m, ids_bead in zip( ['Shared PSD95', 'Unshared PSD95'], [ids_bead_shared_PSD, ids_bead_unshared_PSD] ):
			loc_in_real_coord = positions_real_coord[ids_bead,:]
			loc_in_grid_mesh  = utils.get_hist(loc_in_real_coord)
			conc_in_grid_mesh = ndimage.gaussian_filter(loc_in_grid_mesh, sigma = 2)
			d['locs_in_grid_mesh'][m]  = loc_in_grid_mesh
			d['concs_in_grid_mesh'][m] = conc_in_grid_mesh
		'''
		
		# RDF
		#rdf_grid_points = utils.get_lattice_grids()
		#current_rdfs    = utils.get_a_rdf(types, positions, rdf_grid_points, center )
		#'''
		
		# Load data
		prefix = filename
		suffix = 'sigma_{}'.format(sigma)
		print('Target: {}, sigma: {}'.format(filename, sigma))
		d   = utils.load(dir_edited_data, prefix, suffix)
		
		# Make figure
		fig = make_a_figure(d)
		
		# Save figure
		print(dir_imgs)
		fig.savefig( os.path.join(dir_imgs, filename+'.svg' ) )
		fig.savefig( os.path.join(dir_imgs, filename+'.png' ) , dpi=150)
		#plt.show()
		plt.clf()
		plt.close(fig=fig)

	
	