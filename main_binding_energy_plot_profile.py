
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
	# ax.set_aspect("equal")
	#ax.set_aspect(1.0 / ax.get_data_ratio())
	ax_h, ax_w = ax.bbox.height, ax.bbox.width
	#ax.bbox.width = ax_w / 2
	#print('ax_h, ax_w ', ax_h, ax_w)
	# ax.set_aspect('auto')
	#x, y, w, h = ax.get_position().bounds
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
		axes = []
		column = 1
		ax    = utils.plot_regions_condenstate_from_a_direction_(fig, num_rows, num_columns, row, column, d, transp )
		loc0  = ax.get_position()
		axes.append(ax)
		columns = {'CaMKII':2, 'GluN2B':3, 'STG':4,'PSD95':5}
		axes1 = utils.plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp, \
				pre_rotated=False, colorbar=False )
		axes.extend(axes1)
		
		yloc.append(loc0.y0)
		for a in axes:
			loc = a.get_position()
			#print('loc ', loc)
			a.set_position([loc.x0, yloc[-1], panel_size, panel_size])
	
	
	panel_dx = 0.06
	
	
	# Average concs in CaMKII/STG condensate
	column = 6
	targets_condensate = ['CaMKII','STG']
	for i, t in enumerate(targets_condensate):
		data   = { k: d['conc_condensate'][t][k]  for k in utils.molecules_without_all.keys()}
		colormap_conc_bargraph =[c.cmap_universal_ratio[k] for k in data.keys()]
		ax = fig.add_subplot( num_rows, num_columns, column+num_columns*i )
		ax.set_title('Concentration \n in {} condensate'.format( t ))
		ax.bar(*zip(*data.items()), width=0.6, color=colormap_conc_bargraph )
		ax.set_ylim(0,0.5)
		ax.set_ylabel('(beads / volume)')
		arrange_graph_bar(ax, panel_dx, yloc[i], panel_size/2, panel_size)
		ax.tick_params(axis='x', rotation=45)
	
	
	# Plot RDF
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*2 )
	errorbar= 'shaded' # errorbar='shaded', 'line', or 'final_frame_alone'
	utils.plot_a_rdf( ax, d, errorbar=errorbar, legend=True )
	arrange_graph_bar(ax, panel_dx, yloc[2], panel_size/2, panel_size)
	
	
	# Summed energy in CaMKII/STG condensate
	column = 8
	targets_condensate = ['CaMKII','STG']
	for i, t in enumerate(targets_condensate):
		data   = { k: -d['energy_condensate'][t][k]  for k in utils.molecules_with_all.keys()}
		colormap_conc_bargraph =[c.cmap_universal_ratio[k] for k in data.keys()]
		ax = fig.add_subplot( num_rows, num_columns, column+num_columns*i )
		ax.set_title('Total binding energy \n in {} condensate'.format( t ))
		ax.bar(*zip(*data.items()), width=0.6, color=colormap_conc_bargraph )
		ax.set_ylim(0,80000)
		ax.set_ylabel('Summed (1)')
		arrange_graph_bar(ax, panel_dx, yloc[i], panel_size/2, panel_size)
		ax.tick_params(axis='x', rotation=45)	

	# Volume of  CaMKII/STG condensate
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*2 )
	
	targs = ['CaMKII', 'STG']
	vols  = { t: np.sum( d['region_condensate_in_grid_mesh'][t] ) for t in targs }
	plot_vols_condensate_bar(ax, vols)
	arrange_graph_bar(ax, panel_dx, yloc[2], panel_size/3, panel_size)
	
	
	panel_dx = 0.0
	
	
	# Average energy in CaMKII/STG condensate
	column = 10
	targets_condensate = ['CaMKII','STG']
	for i, t in enumerate(targets_condensate):
		data   = { k: -d['energy_condensate'][t][k] / vols[t]  for k in utils.molecules_with_all.keys()}
		colormap_conc_bargraph =[c.cmap_universal_ratio[k] for k in data.keys()]
		ax = fig.add_subplot( num_rows, num_columns, column+num_columns*i )
		ax.set_title('Binding energy per volume \n in {} condensate (All: {:.3g})'.format( t, data['All'] ))
		ax.bar(*zip(*data.items()), width=0.6, color=colormap_conc_bargraph )
		ax.set_ylim(0,3.5)
		ax.set_ylabel('(1 / volume)')
		arrange_graph_bar(ax, panel_dx, yloc[i], panel_size/2, panel_size)
		ax.tick_params(axis='x', rotation=45)	

	
	return fig
	
	
if __name__ == '__main__':
	
	# Dataset 1
	
	# Input files
	dir_edited_data 		= os.path.join('data', 'binding_energy')
	filenames_edited_data = ['PIPS',\
		'partial', \
		'iPIPS', \
		'homo_valence4', \
		'homo_linear']
	
	# Output files
	dir_imgs = os.path.join('imgs', 'binding_energy','profiles')
	os.makedirs(dir_imgs, exist_ok=True)
	
	sigma = 2
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

