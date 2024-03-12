
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
	

	dir_target  = 'valency_linker_length'
	dir_edited_data	= os.path.join('..','data2', dir_target)
	dir_imgs = os.path.join('imgs2', dir_target, 'binding')
	os.makedirs(dir_imgs, exist_ok=True)
	suffix = 'betweenness'
	
	valency = list(range(2,14,2)) 
	linker_length  = [1, 3, 5, 7, 9]
	
	fnames_valency       = { v: str(v).zfill(2) for v in valency }
	fnames_linker_length = {ll: str(i).zfill(3) for i,ll in enumerate(linker_length) }
	
	num_rows		= len( valency )
	num_columns		= len( linker_length )
	
	vals = np.zeros([num_rows, num_columns])
	fig  = plt.figure(figsize=(10, 10), tight_layout=True)
	#
	
	for i, v in enumerate(valency):
		for j, ll in enumerate(linker_length):
			# Load data
			prefix = fnames_valency[v]+'_'+fnames_linker_length[ll]

			centrality = utils.load(dir_edited_data, prefix, suffix)
			print('Target: {}'.format(prefix))
			# sys.exit(0)
			title = prefix # prefix, None
			row    = num_rows-i-1
			column = j+1
			# print('column: ', column, ', row: ', row)
			ax = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
			plot_a_number_cluster(ax, centrality, prefix)
	
			#if prefix == '12_000':
			#	sys.exit(0)
	
	# Save figure
	filename = suffix
	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename ) ) , dpi=150)
	plt.show()
	plt.clf()
	plt.close(fig=fig)
	
	
	