##
##
##


import os, sys, glob, pickle, pprint, copy

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1


from scipy import ndimage
import networkx as nx

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
from specification_datasets import SpecDatasets



plt.rcParams.update(p.rc_param)
	
	
def arrange_graph_bar(ax, panel_dx, y0, panel_size_x, panel_size_y):
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax_h, ax_w = ax.bbox.height, ax.bbox.width
	loc = ax.get_position()
	if y0 is None:
		y0 = loc.y0
	ax.set_position([loc.x0+panel_dx, y0, panel_size_x, panel_size_y])
	
	
def plot_a_figure( d ):
	
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
		columns = {'CaMKII':1, 'GluN2B':2,'Shared PSD95':3, 'PSD95':4, 'STG':5}
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
	column = 9
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*2 )
	utils.plot_a_rdf_PSD95( ax, d, legend=True ) # , ylim = (-0.006,0.46)
	arrange_graph_bar(ax, panel_dx, yloc[2], panel_size/2, panel_size)
	#'''
	return fig
	
	
	
	
def calc_concs_PSD95_shared_by_STG_GluN2B(dir_lammpstrj, filename_lammpstrj, ids_bead_shared_PSD):
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
	for m, ids_bead in zip( ['Shared PSD95'], [ids_bead_shared_PSD] ):
		loc_in_real_coord = positions_real_coord[ids_bead,:]
		loc_in_grid_mesh  = utils.get_hist(loc_in_real_coord)
		conc_in_grid_mesh = ndimage.gaussian_filter(loc_in_grid_mesh, sigma = 2)
		d['locs_in_grid_mesh'][m]  = loc_in_grid_mesh
		d['concs_in_grid_mesh'][m] = conc_in_grid_mesh
	
	return d
	
	
	
class PlotProfilesSharedPSD95(SpecDatasets):
	def __init__( self ):
		
		pass
		
		
	def run( self ):
		
		# Shared init
		self.dir_imgs = os.path.join(self.dir_imgs_root, 'profiles_shared_PSD95')
		os.makedirs(self.dir_imgs, exist_ok=True)
		print('Img dir: ', self.dir_imgs)
		
		for filename in self.filenames_edited:
			fig = self.make_a_figure( filename )
			self.save_a_figure( fig, filename )
		
		
	def make_a_figure( self, filename ):
		
		# Load graph
		filename_suffix_edited = 'sigma_2'
		filename_suffix_graph  = 'connectivity_graph'
		d_graph = utils.load(self.dir_edited_data, filename, filename_suffix_graph)
		
		multi_graph         = d_graph['multi_graph']
		dir_lammpstrj       = d_graph['dir_lammpstrj']
		filename_lammpstrj  = d_graph['filename']
		ids_bead_shared_PSD = d_graph['ids_bead_shared_PSD']
		sampling_frame      = d_graph['sampling_frame']
		
		# Calc the conc of shared PSD95
		d = calc_concs_PSD95_shared_by_STG_GluN2B(\
			dir_lammpstrj, \
			filename_lammpstrj, \
			ids_bead_shared_PSD)
		
		
		# RDF
		rdf, rdf_bins, rdf_sampling_frames = \
			utils.get_rdfs( dir_lammpstrj, filename_lammpstrj, sampling_frame, multi_graph = multi_graph )
		d['rdf_PSD95_bins'] = rdf_bins
		d['rdf_PSD95_sampling_frames'] = rdf_sampling_frames
		d['rdf_PSD95'] = rdf
		
		
		# Make figure
		return plot_a_figure(d)
		
		
	def save_a_figure( self, fig, filename ):
		# Save figure
		fig.savefig( os.path.join(self.dir_imgs, filename+'.svg' ) )
		fig.savefig( os.path.join(self.dir_imgs, filename+'.png' ) , dpi=150)
		#plt.show()
		plt.clf()
		plt.close(fig=fig)
		
		
if __name__ == '__main__':
	
	obj = PlotProfilesSharedPSD95()
	obj.valency_length() #  conc_dependence(), valency_length(), valency_length_CG(), boundary_conditions2()
	obj.run()
	
	