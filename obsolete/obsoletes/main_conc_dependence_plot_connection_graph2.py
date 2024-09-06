
import os, sys, glob, pickle, pprint
import numpy as np
import math

import matplotlib
import matplotlib.pyplot as plt

import networkx as nx

sys.path.append('../')

import utils
import parameters as p
import colormap as c

plt.rcParams.update(p.rc_param)
plt.rcParams.update( {'font.size': 6} )

def plot_a_number_cluster(ax, centrality, prefix):
	
	'''
	nums_connected = [len(c) for c in nx.connected_components(g)]
	u, counts = np.unique(nums_connected, return_counts=True)
	print('numb  ',u)
	print('freq  ',counts)
	print('\n')
	
	'''
	
	
	ax.grid()
	ax.set_title(prefix)
	#ax.bar(u, counts, width=0.5)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	#ax.set_ylim(0,ymax)
	#ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
	# ax.plot(u, counts, 'k-')
	
	ax.hist(centrality.values(), range = (0.0001, 0.04), bins= 40, log=True )
	ax.set_ylabel('(Number)')
	ax.set_ylabel('(Num of molecules)')
	# ax.set_aspect(1.0 / ax.get_data_ratio())


if __name__ == '__main__':

	
	# Dataset 1
	
	
	dir_target = 'conc_dependence'
	dir_edited_data  = 'conc_dependence'
	dir_edited_data  = os.path.join('..','data2',dir_edited_data)
	
	dir_imgs = os.path.join('imgs2', dir_target)
	os.makedirs(dir_imgs, exist_ok=True)
	suffix = 'betweenness'
	
	
	
		
	
	STG    = [540 , 1620, 2160,  2700,  3240, 4520] 
	GluN2B = [1080, 4320, 8640, 10800, 12960]
	
	volume = np.prod(p.space_np)
	STG    = [ s / volume for s in STG    ]
	GluN2B = [ n / volume for n in GluN2B ]
	
	num_rows		= len( GluN2B )
	num_columns		= len( STG )
	
	vals = np.zeros([num_rows, num_columns])
	fig  = plt.figure(figsize=(10, 10), tight_layout=True)
	#
	for i, stg in enumerate(STG):
		for j, glun in enumerate(GluN2B):
			# Load data
			id = i + j * len(STG)
			prefix = str(id).zfill(3)
			centrality = utils.load(dir_edited_data, prefix, suffix)
			print('Target: {}'.format(prefix))
			
			# sys.exit(0)
			title = prefix # prefix, None
			row    = num_rows-j-1
			column = i+1
			print('column: ', column, ', row: ', row)
			ax = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
			plot_a_number_cluster(ax, centrality, prefix)
	
	# Save figure
	filename = suffix
	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename ) ) , dpi=150)
	plt.show()
	plt.clf()
	plt.close(fig=fig)
	
	
	
	
	
	
	
	
	