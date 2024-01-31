
import os, glob, pickle, pprint, copy

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1

import utils
import colormap as c
from scipy import ndimage


plt.rcParams.update({
                    'pdf.fonttype' : 42,
                    'font.family' : 'sans-serif',
                    'font.sans-serif' : 'Arial',
                    'font.style' : 'normal'})
	
	
def arrange_graph_bar(ax):
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	# ax.set_aspect("equal")
	ax.set_aspect(1.0 / ax.get_data_ratio())
	
	
def plot_concs_condensate_bar(ax, targ, ref, d):
	cc = d['conc_condensate']
	col = c.cmap_universal_ratio[targ]
	ax.bar([targ+'\n condensates', ref+'\n condensates'], [cc[targ][targ], cc[ref][targ]], width=0.5, color=col)
	ax.set_title('Conc of {}'.format(targ))
	ax.set_ylabel('Beads / volume')
	ax.set_ylim(0,0.6)
	ax.tick_params(axis='x', rotation=45)
	arrange_graph_bar(ax)
	
	
def plot_conc_ratio_condensate_bar(ax, targs, counterparts, d):
	cc = d['conc_condensate']
	conc_ratio = [ cc[t][t] / cc[c][t] for t, c in zip(targs, counterparts)]
	cols = [ c.cmap_universal_ratio[targ] for targ in targs ]
	ax.bar(targs, conc_ratio, width=0.5, color=cols)
	ax.set_title('Partition index \n between condensates')
	ax.set_ylabel('Target / counterpart')
	ax.set_ylim(0,40)
	# ax.tick_params(axis='x', rotation=45)
	arrange_graph_bar(ax)
	
	
def make_a_figure(d, dir_imgs, filename):
	
	# Preparation
	transps = [(0,1,2),(1,2,0),(2,0,1)]
	titles  = [ 'yz', 'xy', 'zx']
	num_columns = 9
	num_rows    = 3
	fig = plt.figure(figsize=(18, 8), tight_layout=True)
	
	# Plot profiles
	for row,(title, transp) in enumerate(zip(titles, transps)):
		column = 1
		utils.plot_regions_condenstate_from_a_direction(fig, num_rows, num_columns, row, column, d, transp )
		columns = {'CaMKII':2, 'GluN2B':3, 'STG':4,'PSD95':5}
		utils.plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp )
		columns = {'CaMKII':6, 'STG':7}
		utils.plot_watershed_region_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp)
	
	
	# Plot the volume of watershed basin
	ratio_volumes_watershed = {k: v*100 for k, v in d['ratio_volumes_watershed'].items()}
	cols = [c.cmap_universal_ratio[k] for k in ratio_volumes_watershed.keys()]
	ax = fig.add_subplot( num_rows, num_columns, 8 )
	ax.bar(ratio_volumes_watershed.keys(), ratio_volumes_watershed.values(), width=0.5, color=cols) # , color=cols
	ax.set_title('Volume of \n watershed basin')
	ax.set_ylabel('/ total system volume (%)')
	ax.set_ylim(0,3.0)
	arrange_graph_bar(ax)
	
	
	# Plot concs in condensates
	column = 9
	targ = 'STG'
	ref  = 'CaMKII'
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*1 )
	plot_concs_condensate_bar(ax, targ, ref, d)
	
	targ = 'CaMKII'
	ref  = 'STG'
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*2 )
	plot_concs_condensate_bar(ax, targ, ref, d)
	
	targs        = ['STG','CaMKII']
	counterparts = ['CaMKII','STG']
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*0 )
	plot_conc_ratio_condensate_bar(ax, targs, counterparts, d)
	
	return fig
	
	
if __name__ == '__main__':
	
	'''
	# Dataset 1
	filenames	= [	'PIPS',\
					'iPIPS',\
					'PartialE',\
					'Homo']	
	dir_data	= 'data'
	dir_imgs	= 'imgs/profiles_example'
	'''
	
	#'''
	# Dataset 2
	filenames   = [str(i).zfill(3) for i in range(30) ]
	dir_data	= 'data'
	dir_imgs	= 'imgs/profiles_conc'
	#'''
	
	os.makedirs(dir_imgs, exist_ok=True)
	sigmas = [2,3,4]
	for sigma in sigmas:
		for filename in filenames:
			# Load data
			prefix = filename
			suffix = 'sigma_{}'.format(sigma)
			print('Target: {}, sigma: {}'.format(filename, sigma))
			d      = utils.load(dir_data, prefix, suffix)
			
			# Process data
			# labels_watershed, volume_ratios_watershed = watershed_segmentation( d )
			fig = make_a_figure(d, dir_imgs, filename)
			
			# Save figure
			fig.savefig( os.path.join(dir_imgs, '{}_sigma_{}.pdf'.format( filename, sigma ) ) )
			fig.savefig( os.path.join(dir_imgs, '{}_sigma_{}.png'.format( filename, sigma ) ) , dpi=150)
			#plt.show()
			plt.clf()
			plt.close(fig=fig)
