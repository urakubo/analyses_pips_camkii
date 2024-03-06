
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
		axes1 = utils.plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp, pre_rotated=False )
		axes.extend(axes1)
		'''
		columns = {'CaMKII':6, 'STG':7}
		axes2 = utils.plot_watershed_region_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp)
		axes.extend(axes2)
		'''
		
		yloc.append(loc0.y0)
		for a in axes:
			loc = a.get_position()
			#print('loc ', loc)
			a.set_position([loc.x0, yloc[-1], panel_size, panel_size])
	
	panel_dx = 0.03
	
	
	# Plot the volume of watershed basin
	'''
	ratio_volumes_watershed = {k: v*100 for k, v in d['ratio_volumes_watershed'].items()}
	cols = [c.cmap_universal_ratio[k] for k in ratio_volumes_watershed.keys()]
	ax = fig.add_subplot( num_rows, num_columns, 8 )
	ax.bar(ratio_volumes_watershed.keys(), ratio_volumes_watershed.values(), width=0.5, color=cols) # , color=cols
	ax.set_title('Volume of \n watershed basin')
	ax.set_ylabel('/ total system volume (%)')
	ax.set_ylim(0,3.0)
	arrange_graph_bar(ax, panel_dx, yloc[0], panel_size/4, panel_size )
	'''
	
	
	# Plot concs in condensates
	column = 9
	targ = 'STG'
	ref  = 'CaMKII'
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*1 )
	utils.plot_conc_condensate_bar(ax, targ, ref, d)
	arrange_graph_bar(ax, panel_dx, yloc[1], panel_size/4, panel_size)
	
	targ = 'CaMKII'
	ref  = 'STG'
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*2 )
	utils.plot_conc_condensate_bar(ax, targ, ref, d)
	arrange_graph_bar(ax, panel_dx, yloc[2], panel_size/4, panel_size)
	
	targs        = ['STG','CaMKII']
	counterparts = ['CaMKII','STG']
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*0 )
	utils.plot_conc_ratio_condensate_bar(ax, targs, counterparts, d)
	arrange_graph_bar(ax, panel_dx, yloc[0], panel_size/4, panel_size )
	
	
	# Average concs in CaMKII/ condensate
	column = 10
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
	
	return fig
	
	
if __name__ == '__main__':
	
	# Datasets
	
	# Conc dependence
	#'''
	dir_target  = 'conc_dependence'
	filenames_edited_data 	= [str(i).zfill(3) for i in range(30) ] # 70
	#'''
	
	# Linker length
	'''
	filenames_edited_data = [str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(5) ]
	dir_target  = 'valency_linker_length'
	'''
	# Figure 1
	#'''
	filenames_edited_data = ['CG','CPG','SP']
	dir_target  = 'types_mixture'
	#'''
	
	# Shared init
	dir_edited_data	= os.path.join('data2', dir_target)
	dir_imgs = os.path.join('imgs2', dir_target,'profiles')
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
		fig.savefig( os.path.join(dir_imgs, '{}_sigma_{}.svg'.format( filename, sigma ) ) )
		fig.savefig( os.path.join(dir_imgs, '{}_sigma_{}.png'.format( filename, sigma ) ) , dpi=150)
		#plt.show()
		plt.clf()
		plt.close(fig=fig)

