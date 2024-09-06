
import os, sys, glob, pickle, pprint, copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import mpl_toolkits.axes_grid1

from scipy.interpolate import griddata, RegularGridInterpolator
import utils
import colormap as c
import parameters as p


plt.rcParams.update(p.rc_param)

	
	
if __name__ == '__main__':
	
	
	'''	
	data  = ave_num_binding_GluN2B
	title = 'Number of GluN2B bound to one PSD95'
	filename_output = 'num_GluN2B_bound_to_one_PSD95'
	colormap   =  c.cmap_white_purple_universal
	levels  = np.linspace(0,3,10)
	
	data  = ave_num_binding_GluN2B_CaMKII
	title = 'Number of GluN2B bound to one CaMKII'
	filename_output = 'num_GluN2B_bound_to_one_CaMKII'
	colormap   =  c.cmap_white_green_universal
	levels  = np.linspace(0,12,7)
	
	data  = ave_num_binding_PSD95_STG
	title = 'Number of PSD95 bound to one STG'
	filename_output = 'num_PSD95_bound_to_one_STG'
	colormap   = c.cmap_white_red_universal
	levels  = np.linspace(0,4,6)
	
	
	data  = ave_num_binding_STG
	title = 'Number of STG bound to one PSD95'
	filename_output = 'num_STG_bound_to_one_PSD95'
	colormap   = c.cmap_white_red_universal
	levels  = np.linspace(0,3,6)
	'''
	#'''
	data  = ratio_binding_Both
	title = 'Ratio of PSD95 bound to both GluN2B and STG'
	filename_output = 'PSD95_bound_to_both_GluN2B_STG'
	colormap   =  plt.colormaps['Greys']
	#levels  = np.linspace(0,0.8,9)
	levels  = np.linspace(0,1.0,11)
	#'''
	
	
	# Dataset 1
	species = 'PSD95_connection' # 'STG','GluN2B', 'PSD95','CaMKII', 'PSD95_connection'
	
	dir_target = 'conc_dependence'
	
	dir_edited_data  = os.path.join('data3',dir_target)
	
	dir_imgs = os.path.join('imgs3', dir_target)
	os.makedirs(dir_imgs, exist_ok=True)
	suffix = 'connection'
	
	
	STG    = [540, 1620, 2160, 2700, 3240, 4520] 
	GluN2B = [570, 1080, 4320, 6480, 8640, 10800, 12960, 17280]
	
	volume = np.prod(p.space_np)
	STG    = [ s / volume for s in STG    ]
	GluN2B = [ n / volume for n in GluN2B ]
	
	STG    = np.array( STG ) * 1000
	GluN2B = np.array( GluN2B ) * 100
	
	num_rows		= len( GluN2B )
	num_columns		= len( STG )
	
	
	#
	data = np.zeros([num_columns, num_rows], dtype = 'float')
	#
	for i, stg in enumerate(STG):
		for j, glun in enumerate(GluN2B):
			# Load data
			id = i + j * len(STG)
			prefix = str(id).zfill(3)
			d      = utils.load(dir_edited_data, prefix, suffix)
			print('Target: {}'.format(prefix))
			
			# sys.exit(0)
			title = prefix # prefix, None
			row    = num_rows-j-1
			column = i+1
			print('column: ', column, ', row: ', row)

			species = 'PSD95'
			ave_num_binding_GluN2B[i,j] = np.average(d[species]['nums']['GluN2B_PSD95'])
			ave_num_binding_STG[i,j] = np.average(d[species]['nums']['STG_PSD95'])

			species = 'CaMKII'
			ave_num_binding_GluN2B_CaMKII[i,j] = np.average(d[species]['nums']['GluN2B_CaMKII'])

			species = 'STG'
			ave_num_binding_PSD95_STG[i,j] = np.average(d[species]['nums']['STG_PSD95'])

			species = 'PSD95_connection'
			total_number = sum(d[species].values())
			ratio_binding_None[i,j]   = d[species]['None'] / total_number
			ratio_binding_STG[i,j]    = d[species]['STG only'] / total_number
			ratio_binding_GluN2B[i,j] = d[species]['GluN2B only'] / total_number
			ratio_binding_Both[i,j]   = d[species]['Both'] / total_number
	
	
	
	#colormap   =  plt.colormaps['Oranges']
	
	
	fig  = plt.figure(figsize=(5, 5))
	fig.subplots_adjust(wspace=0.4,  hspace=0.6)
	
	
	# cmap_gray_cr_pk_gray # c.cmap_white_green_universal, plt.colormaps['jet']# 'winter', etc
	
	
	ax = fig.add_subplot( 1, 1, 1 )
	cs, cb = utils.plot_a_panel(ax, data, STG, GluN2B, colormap, levels)
	ax.set_title( title )
	ax.set_xlabel('STG (beads / voxel) x 10-3')
	ax.set_ylabel('GluN2B (beads / voxel) x 10-2')
		
	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename_output ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename_output ) ) , dpi=150)
	
	plt.show()
	plt.clf()
	plt.close(fig=fig)
	
	
	