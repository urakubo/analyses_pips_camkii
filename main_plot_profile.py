
import os, glob, pickle, pprint, copy

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1

import utils
import colormap as c
from scipy import ndimage


plt.rcParams.update({
                    'pdf.fonttype' : 'truetype',
                    'svg.fonttype' : 'none',
                    'font.family' : 'sans-serif',
                    'font.sans-serif' : 'Arial',
                    'font.style' : 'normal'})
	
	
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
	
	
def plot_concs_condensate_bar(ax, targ, ref, d):
	cc = d['conc_condensate']
	col = c.cmap_universal_ratio[targ]
	ax.bar(['In '+targ+'\n condensates', 'In '+ref+'\n condensates'], [cc[targ][targ], cc[ref][targ]], width=0.5, color=col)
	ax.set_title('Conc of {}'.format(targ))
	ax.set_ylabel('(beads / volume)')
	ax.set_ylim(0,0.6)
	ax.tick_params(axis='x', rotation=45)
	
	
def plot_conc_ratio_condensate_bar(ax, targs, counterparts, d):
	cc = d['conc_condensate']
	conc_ratio = [ cc[t][t] / cc[c][t] for t, c in zip(targs, counterparts)]
	cols = [ c.cmap_universal_ratio[targ] for targ in targs ]
	ax.bar(targs, conc_ratio, width=0.5, color=cols)
	ax.set_title('Partition index')
	ax.set_ylabel('(target / counterpart)')
	ax.set_ylim(0,40)
	
	
def make_a_figure( d ):
	
	# Preparation
	transps = [(0,1,2),(1,2,0),(2,0,1)]
	titles  = [ 'yz', 'xy', 'zx']
	num_columns = 10
	num_rows    = 3
	fig = plt.figure(figsize=(20, 8)) # , tight_layout=False
	
	left, right, bottom, top = 0.0, 0.95, 0.10, 0.99
	wspace, hspace = 0.2, 0.1
	plt.subplots_adjust(left, bottom, right, top, wspace, hspace)
	
	panel_size= 0.15
	
	# Plot profiles
	yloc = []
	for row,(title, transp) in enumerate(zip(titles, transps)):
		axes = []
		column = 1
		ax    = utils.plot_regions_condenstate_from_a_direction(fig, num_rows, num_columns, row, column, d, transp )
		loc0  = ax.get_position()
		axes.append(ax)
		columns = {'CaMKII':2, 'GluN2B':3, 'STG':4,'PSD95':5}
		axes1 = utils.plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp )
		axes.extend(axes1)
		columns = {'CaMKII':6, 'STG':7}
		axes2 = utils.plot_watershed_region_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp)
		axes.extend(axes2)
		
		yloc.append(loc0.y0)
		for a in axes:
			loc = a.get_position()
			#print('loc ', loc)
			a.set_position([loc.x0, yloc[-1], panel_size, panel_size])
	
	
	panel_dx = 0.03
	
	# Plot the volume of watershed basin
	ratio_volumes_watershed = {k: v*100 for k, v in d['ratio_volumes_watershed'].items()}
	cols = [c.cmap_universal_ratio[k] for k in ratio_volumes_watershed.keys()]
	ax = fig.add_subplot( num_rows, num_columns, 8 )
	ax.bar(ratio_volumes_watershed.keys(), ratio_volumes_watershed.values(), width=0.5, color=cols) # , color=cols
	ax.set_title('Volume of \n watershed basin')
	ax.set_ylabel('/ total system volume (%)')
	ax.set_ylim(0,3.0)
	arrange_graph_bar(ax, panel_dx, yloc[0], panel_size/4, panel_size )
	
	
	# Plot concs in condensates
	column = 9
	targ = 'STG'
	ref  = 'CaMKII'
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*1 )
	plot_concs_condensate_bar(ax, targ, ref, d)
	arrange_graph_bar(ax, panel_dx, yloc[1], panel_size/4, panel_size)
	
	targ = 'CaMKII'
	ref  = 'STG'
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*2 )
	plot_concs_condensate_bar(ax, targ, ref, d)
	arrange_graph_bar(ax, panel_dx, yloc[2], panel_size/4, panel_size)
	
	targs        = ['STG','CaMKII']
	counterparts = ['CaMKII','STG']
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*0 )
	plot_conc_ratio_condensate_bar(ax, targs, counterparts, d)
	arrange_graph_bar(ax, panel_dx, yloc[0], panel_size/4, panel_size )
	
	
	# Average concs in CaMKII/ condensate
	column = 10
	targets_condensate = ['CaMKII','STG']
	for i, t in enumerate(targets_condensate):
		data   = { k: d['conc_condensate'][t][k]  for k in utils.molecules_without_all.keys()}
		colormap_conc_bargraph =[c.cmap_universal_ratio[k] for k in data.keys()]
		ax = fig.add_subplot( num_rows, num_columns, column+num_columns*i )
		ax.set_title('{} condensate'.format( t ))
		ax.bar(*zip(*data.items()), width=0.6, color=colormap_conc_bargraph )
		ax.set_ylim(0,0.5)
		ax.set_ylabel('(beads / volume)')
		arrange_graph_bar(ax, panel_dx, yloc[i], panel_size/2, panel_size)
		ax.tick_params(axis='x', rotation=45)
	
	
	# Plot RDF
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*2 )
	r = d['rdf_bins'][1:-1]
	for k in utils.molecules_without_all.keys():
		rdf_mean  = np.mean( d['rdf'][k][1:], axis = 1)
		rdf_error = np.std(  d['rdf'][k][1:], axis = 1)
		ax.fill_between(r, rdf_mean-rdf_error, rdf_mean+rdf_error,color=c.cmap_universal_ratio_light[k])
		ax.plot(r, rdf_mean, color=c.cmap_universal_ratio[k])
	ax.set_xlabel('Distance from \n center-of-mass (l.u.)')
	ax.set_ylabel('(beads / volume)')
	ax.set_xlim(0,30)
	ax.set_ylim(-0.006,0.66)
	ax.set_xticks(range(0,40,10))
	arrange_graph_bar(ax, panel_dx, yloc[2], panel_size/2, panel_size)
	
	
	return fig
	
	
if __name__ == '__main__':
	
	#'''
	# Dataset 1
	filenames	= [	'PIPS',\
					'iPIPS',\
					'PartialE',\
					'Homo']	
	dir_data	= 'data'
	dir_imgs	= 'imgs/profiles_example'
	#'''

	'''
	# Dataset 2
	filenames   = [str(i).zfill(3) for i in range(30) ]
	dir_data	= 'data'
	dir_imgs	= 'imgs/profiles_conc'
	'''
	
	os.makedirs(dir_imgs, exist_ok=True)
	sigmas = [2,3,4]
	for sigma in sigmas:
		for filename in filenames:
			# Load data
			prefix = filename
			suffix = 'sigma_{}'.format(sigma)
			print('Target: {}, sigma: {}'.format(filename, sigma))
			d   = utils.load(dir_data, prefix, suffix)
			
			# Process data
			# labels_watershed, volume_ratios_watershed = watershed_segmentation( d )
			fig = make_a_figure(d)
			
			# Save figure
			fig.savefig( os.path.join(dir_imgs, '{}_sigma_{}_2.svg'.format( filename, sigma ) ) )
			fig.savefig( os.path.join(dir_imgs, '{}_sigma_{}_2.png'.format( filename, sigma ) ) , dpi=150)
			#plt.show()
			plt.clf()
			plt.close(fig=fig)
