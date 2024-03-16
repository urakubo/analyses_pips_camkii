
import os, sys, glob, pickle, pprint
import numpy as np
import math

import matplotlib
import matplotlib.pyplot as plt

sys.path.append('../')

import utils
import parameters as p
import colormap as c

plt.rcParams.update(p.rc_param)
plt.rcParams.update( {'font.size': 5} )



def plot_a_binding(ax, multi_graph, prefix, species, distribution_or_average = 'distribution'):
	
	
	ids = [ i for i, attr in multi_graph.nodes('species') if attr == species ]
	
	if species == 'CaMKII' and distribution_or_average == 'distribution':
		numbers_of_connections = [ multi_graph.degree[id] for id in ids ]
		dist = {'{:d} GluN2B'.format(i): numbers_of_connections.count(i) for i in range(13)}
		ymax = 1000
		#nx.get_node_attributes(multi_graph, 'species').keys()
	elif species == 'CaMKII' and distribution_or_average == 'average':
		numbers_of_connections = [ multi_graph.degree[id] for id in ids ]
		dist = {'GluN2B \n {:.2f}'.format(ave_GluN2B): np.average(numbers_of_connections)}
		ymax = 12
	if species == 'STG' and distribution_or_average == 'distribution':
		numbers_of_connections = [ multi_graph.degree[id] for id in ids ]
		dist = {'{:d} PSD95'.format(i): numbers_of_connections.count(i) for i in range(5)}
		ymax = 4000
		#nx.get_node_attributes(multi_graph, 'species').keys()
	elif species == 'GluN2B' and distribution_or_average == 'distribution':
		ymax = 8000
		reference_types_connection = [\
			[], \
			['GluN2B_CaMKII'], \
			['GluN2B_CaMKII', 'GluN2B_CaMKII'], \
			['GluN2B_PSD95'], \
			['GluN2B_PSD95', 'GluN2B_PSD95'], \
			['GluN2B_CaMKII', 'GluN2B_PSD95']]
		titles = [\
			'None', \
			'CaMKII, None',\
			'CaMKII, CaMKII',\
			'PSD95, None', \
			'PSD95, PSD95', \
			'CaMKII, PSD95' ]
		dist = {t: 0 for t in titles}
		edges_from_molecules = [ multi_graph.edges(id, keys=True, data=True) for id in ids ]
		connections_from_one_molecules = [ [e[3]['type_connection'] for e in es] for es in edges_from_molecules ]
		for connections_from_a_molecule in connections_from_one_molecules:
			for j, ref in enumerate( reference_types_connection ):
				if utils.equal_list(connections_from_a_molecule, ref):
					dist[titles[j]] += 1
	elif species == 'PSD95' and distribution_or_average == 'distribution':
		ymax = 3000
		reference_types_connection = [\
			[], \
			['STG_PSD95'],['STG_PSD95']*2,['STG_PSD95']*3, \
			['GluN2B_PSD95'], ['GluN2B_PSD95']*2, ['GluN2B_PSD95']*3,\
			['STG_PSD95', 'GluN2B_PSD95'],\
			['STG_PSD95']*2 + ['GluN2B_PSD95'],\
			['STG_PSD95'] + ['GluN2B_PSD95']*2]
		titles = [\
			'None', \
			'1 STG', '2 STG', '3 STG', \
			'1 GluN2B', '2 GluN2B', '3 GluN2B', \
			'1 STG, 1 GluN2B', '2 STG, 1 GluN2B', '1 STG, 2 GluN2B'
			 ]
		dist = {t: 0 for t in titles}
		edges_from_molecules = [ multi_graph.edges(id, keys=True, data=True) for id in ids ]
		connections_from_one_molecules = [ [e[3]['type_connection'] for e in es] for es in edges_from_molecules ]
		for connections_from_a_molecule in connections_from_one_molecules:
			for j, ref in enumerate( reference_types_connection ):
				if utils.equal_list(connections_from_a_molecule, ref):
					dist[titles[j]] += 1

	elif species == 'PSD95' and distribution_or_average == 'average':

		numbers_of_connections = [ multi_graph.degree[id] for id in ids ]
		dist = {'GluN2B \n {:.2f}'.format(ave_GluN2B): np.average(numbers_of_connections)}
		ymax = 12

		edges_from_molecules = [ multi_graph.edges(id, keys=True, data=True) for id in ids ]
		connections_from_one_molecules = [ [e[3]['type_connection'] for e in es] for es in edges_from_molecules ]
		
		
		ymax = 3000
		reference_types_connection = [\
			[], \
			['STG_PSD95'],['STG_PSD95']*2,['STG_PSD95']*3, \
			['GluN2B_PSD95'], ['GluN2B_PSD95']*2, ['GluN2B_PSD95']*3,\
			['STG_PSD95', 'GluN2B_PSD95'],\
			['STG_PSD95']*2 + ['GluN2B_PSD95'],\
			['STG_PSD95'] + ['GluN2B_PSD95']*2]
		titles = [\
			'None', \
			'1 STG', '2 STG', '3 STG', \
			'1 GluN2B', '2 GluN2B', '3 GluN2B', \
			'1 STG, 1 GluN2B', '2 STG, 1 GluN2B', '1 STG, 2 GluN2B'
			 ]
		dist = {t: 0 for t in titles}
		edges_from_molecules = [ multi_graph.edges(id, keys=True, data=True) for id in ids ]
		connections_from_one_molecules = [ [e[3]['type_connection'] for e in es] for es in edges_from_molecules ]
		for connections_from_a_molecule in connections_from_one_molecules:
			for j, ref in enumerate( reference_types_connection ):
				if utils.equal_list(connections_from_a_molecule, ref):
					dist[titles[j]] += 1


	ax.set_title(prefix)
	ax.bar(*zip(*dist.items()), width=0.5)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_ylim(0,ymax)
	ax.set_ylabel('(Number)')
	# ax.set_aspect(1.0 / ax.get_data_ratio())

	if distribution_or_average == 'distribution':
		ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

if __name__ == '__main__':

	
	# Dataset 1
	species = 'STG' # 'STG','GluN2B', 'PSD95','CaMKII', 'PSD95_connection'
	dir_target  = 'conc_dependence'
	dir_target2 = 'nums_connection'
	
	dir_edited_data  = os.path.join('data3',dir_target)
	
	dir_imgs = os.path.join('imgs3', dir_target, dir_target2)
	os.makedirs(dir_imgs, exist_ok=True)
	suffix = 'graph2'
	
	
	STG    = [540, 1620, 2160, 2700, 3240, 4520] 
	GluN2B = [570, 1080, 4320, 6480, 8640, 10800, 12960, 17280]
	
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
			g      = utils.load(dir_edited_data, prefix, suffix)
			print('Target: {}'.format(prefix))
			graph = g[0]
			multi_graph = g[1]
			#sys.exit(0)
			title = prefix # prefix, None
			row    = num_rows-j-1
			column = i+1
			print('column: ', column, ', row: ', row)
			ax = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
			plot_a_binding(ax, multi_graph, prefix, species)
	
	# Save figure
	filename = 'ave_num_binding_{}'.format(species)
	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename ) ) , dpi=150)
	plt.show()
	plt.clf()
	plt.close(fig=fig)
	
	
	