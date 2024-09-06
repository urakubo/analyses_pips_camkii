
import os, sys, glob, pickle, pprint
import numpy as np
import math

import matplotlib
import matplotlib.pyplot as plt

import utils
import parameters as p
import colormap as c

plt.rcParams.update(p.rc_param)
plt.rcParams.update( {'font.size': 5} )
	
if __name__ == '__main__':

	
	# Dataset 1
	species       = 'CaMKII' # 'STG','GluN2B', 'PSD95','CaMKII'
	type_analysis = 'distribution'
	# 'average' and 'distribution' for all,
	# 'CaMKII' and 'PSD95' for GluN2B
	# 'ratio' for PSD95
	
	ymax = {'STG'   : {'average': 4, 'distribution':4000},\
			'GluN2B': {'average': 2, 'distribution':8000}, \
			'PSD95' : {'average': 3, 'distribution':3000},
			'CaMKII': {'average':12, 'distribution':1000}}
	ymax['GluN2B']['CaMKII'] =  10000
	ymax['GluN2B']['PSD95']  =  10000
	ymax['PSD95']['ratio']   =  10000
	
	dir_target  = 'conc_dependence'
	dir_target2 = 'connectivity'
	
	dir_edited_data  = os.path.join('data3',dir_target)
	
	dir_imgs = os.path.join('imgs3', dir_target, dir_target2)
	os.makedirs(dir_imgs, exist_ok=True)
	suffix = 'connectivity_graph'
	
	
	STG    = [540, 1620, 2160, 2700, 3240, 4520] 
	GluN2B = [570, 1080, 4320, 6480, 8640, 10800, 12960, 17280]
	
	volume = np.prod(p.space_np)
	STG    = [ s / volume for s in STG    ]
	GluN2B = [ n / volume for n in GluN2B ]
	
	num_rows		= len( GluN2B )
	num_columns		= len( STG )
	
	vals = np.zeros([num_rows, num_columns])
	fig  = plt.figure(figsize=(10, 10), tight_layout=True)
	plt.suptitle('Species: {}, {}'.format( species, type_analysis ), )
	#
	for i, stg in enumerate(STG):
		for j, glun in enumerate(GluN2B):
			# Load data
			id = i + j * len(STG)
			prefix = str(id).zfill(3)
			d      = utils.load(dir_edited_data, prefix, suffix)
			print('Target: {}'.format(prefix))
			
			dist = d[species][type_analysis]
			
			title = prefix # prefix, None
			row    = num_rows-j-1
			column = i+1
			print('column: ', column, ', row: ', row)
			ax = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
			ax.set_title(prefix)
			ax.bar(*zip(*dist.items()), width=0.5)
			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)
			ax.set_ylim(0,ymax[species][type_analysis])
			ax.set_ylabel('(Number)')
			# ax.set_aspect(1.0 / ax.get_data_ratio())

			if  type_analysis in ['distribution', 'CaMKII', 'PSD95'] :
				ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

	
	# Save figure
	filename = '{}_{}'.format(species, type_analysis)
	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename ) ) , dpi=150)
	plt.show()
	plt.clf()
	plt.close(fig=fig)
	
	
	