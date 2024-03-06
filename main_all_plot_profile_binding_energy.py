
import os, glob, pickle, pprint, copy

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1

import utils
import colormap as c
import parameters as p
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
	
	
def plot_vols_condensate_bar(ax, vols):
	cols = [ c.cmap_universal_ratio[t] for t in vols.keys() ]
	ax.bar(*zip(*vols.items()), width=0.5, color=cols )
	ax.set_title('Condensate volume')
	ax.set_ylabel('(lattice units)')
	# ax.set_ylim(0,0.6)
	ax.tick_params(axis='x', rotation=45)
	
	
	
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
	fig.suptitle( "Time step: {}, MC step: {:e}".format( d['time_frame'], d['mc_step'] ) )
	
	panel_size = 0.15
	
	# Plot profiles
	yloc = []
	for row,(title, transp) in enumerate(zip(titles, transps)):
		columns = {'CaMKII':1, 'GluN2B':2, 'STG':3,'PSD95':4}
		axes = utils.plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp, \
				pre_rotated=False, colorbar=False )
		loc0  = axes[0].get_position()
		
		yloc.append(loc0.y0)
		for a in axes:
			loc = a.get_position()
			#print('loc ', loc)
			a.set_position([loc.x0, yloc[-1], panel_size, panel_size])
	
	panel_dx = 0.06
	
	
	# Average concs in CaMKII/ condensate
	column = 5
	targets_condensate = ['CaMKII','STG']
	for i, t in enumerate(targets_condensate):
		ax = fig.add_subplot( num_rows, num_columns, column+num_columns*i )
		utils.plot_concs_condensate_bar(ax, t, d)
		arrange_graph_bar(ax, panel_dx, yloc[i], panel_size/2, panel_size)

	
	# Plot RDF
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*2 )
	errorbar= 'shaded' # errorbar='shaded', 'line', or 'final_frame_alone'
	utils.plot_a_rdf( ax, d, errorbar=errorbar, legend=True )
	arrange_graph_bar(ax, panel_dx, yloc[2], panel_size/2, panel_size)
	
	
	
	# Average anisotropic energy in CaMKII/STG condensate
	
	panel_dx = 0.0
	
	column = 8
	targets_condensate = ['CaMKII','STG']
	ymax = 1.5
	
	type_energy = 'energy_anisotropic'
	for i, target_condensate in enumerate(targets_condensate):
		ax = fig.add_subplot( num_rows, num_columns, column+num_columns*i )
		utils.plot_binding_energy_bar(ax,d, type_energy, target_condensate, ymax)
		arrange_graph_bar(ax, panel_dx, yloc[i], panel_size/2, panel_size)
		
	# Volume of  CaMKII/STG condensate
	ax   = fig.add_subplot( num_rows, num_columns, column+num_columns*2 )
	vols = { t: np.sum( d['region_condensate_in_grid_mesh'][t] ) for t in targets_condensate }
	plot_vols_condensate_bar(ax, vols)
	arrange_graph_bar(ax, panel_dx, yloc[2], panel_size/3, panel_size)
	
	
	column = 10
	targets_condensate = ['CaMKII','STG']
	type_energy = 'energy_isotropic'
	for i, target_condensate in enumerate(targets_condensate):
		ax = fig.add_subplot( num_rows, num_columns, column+num_columns*i )
		utils.plot_binding_energy_bar(ax,d, type_energy, target_condensate, ymax)
		arrange_graph_bar(ax, panel_dx, yloc[i], panel_size/2, panel_size)
	
	return fig
	
	
if __name__ == '__main__':
	
	# Dataset 1
	
	# Linker length
	filenames_edited_data = [str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(5) ]
	dir_target  = 'valency_linker_length'
	
	
	# Shared init
	dir_edited_data	= os.path.join('data2', dir_target)
	dir_imgs = os.path.join('imgs2', dir_target,'binding_energy')
	sigma = 2
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	for filename in filenames_edited_data:
		# Load data
		prefix = filename
		suffix = 'sigma_{}'.format(sigma)
		print('Target: {}, sigma: {}'.format(filename, sigma))
		d   = utils.load(dir_edited_data, prefix, suffix)
		
		# Make figure
		fig = make_a_figure(d)
		
		# Save figure
		fig.savefig( os.path.join(dir_imgs, filename + '.svg' ) )
		fig.savefig( os.path.join(dir_imgs, filename + '.png' ) , dpi=150)
		#plt.show()
		plt.clf()
		plt.close(fig=fig)

