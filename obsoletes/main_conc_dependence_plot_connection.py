
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



def get_connection_statistics(multi_graph, species, type_analysis):
	
	ids = [ i for i, attr in multi_graph.nodes('species') if attr == species ]
	numbers_of_connections = [ multi_graph.degree[id] for id in ids ]
	
	edges_from_molecules = [ multi_graph.edges(id, keys=True, data=True) for id in ids ]
	connections_from_one_molecules = [ [e[3]['type_connection'] for e in es] for es in edges_from_molecules ]
	
	if species == 'GluN2B' and type_analysis in ['CaMKII','PSD95']:
		ymax = 10000
		types_connection = [[len(v) for k, v in multi_graph[id].items() if multi_graph.nodes[k]['species'] == type_analysis] for id in ids]
		reference_types_connection = {\
			'No {}'.format(type_analysis) :[],\
			'One {} (single)'.format(type_analysis) :[1],\
			'One {} (double)'.format(type_analysis) :[2],\
			'Two {}s'.format(type_analysis) :[1, 1]\
			}
		dist = {k: types_connection.count(v) for k, v in reference_types_connection.items()}
		
	elif species == 'CaMKII' and type_analysis == 'distribution':
		ymax = 1000
		dist = {'{:d} GluN2B'.format(i): numbers_of_connections.count(i) for i in range(13)}
		
	elif species == 'CaMKII' and type_analysis == 'average':
		ymax = 12
		ave_numbers_of_connections = np.average( numbers_of_connections )
		dist = {'GluN2B \n {:.2f}'.format(ave_numbers_of_connections): ave_numbers_of_connections}
		
	elif species == 'STG' and type_analysis == 'distribution':
		ymax = 4000
		dist = {'{:d} PSD95'.format(i): numbers_of_connections.count(i) for i in range(5)}
		
	elif species == 'STG' and type_analysis == 'average':
		ymax = 4
		dist = {'PSD95 \n {:.2f}'.format(numbers_of_connections): numbers_of_connections}
		
	elif species == 'GluN2B' and type_analysis == 'average':
		ymax = 4
		num_GluN2B_CaMKII = np.average( [ c.count('GluN2B_CaMKII') for c in connections_from_one_molecules ] )
		num_GluN2B_PSD95  = np.average( [ c.count('GluN2B_PSD95') for c in connections_from_one_molecules ] )
		dist = {'GluN2B_CaMKII \n {:.2f}'.format(num_GluN2B_CaMKII): num_GluN2B_CaMKII,
				'GluN2B_PSD95  \n {:.2f}'.format(num_GluN2B_PSD95) : num_GluN2B_PSD95 }
		
	elif species == 'GluN2B' and type_analysis == 'distribution':
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
		for connections_from_a_molecule in connections_from_one_molecules:
			for j, ref in enumerate( reference_types_connection ):
				if utils.equal_list(connections_from_a_molecule, ref):
					dist[titles[j]] += 1

	elif species == 'PSD95' and type_analysis == 'average':
		ymax = 4
		num_STG_PSD95    = np.average( [c.count('STG_PSD95')  for c in connections_from_one_molecules ]   )
		num_GluN2B_PSD95 = np.average( [c.count('GluN2B_PSD95') for c in connections_from_one_molecules ] )
		dist = {'STG_PSD95 \n {:.2f}'.format(num_STG_PSD95): num_STG_PSD95, \
				'GluN2B_PSD95 \n {:.2f}'.format(num_GluN2B_PSD95): num_GluN2B_PSD95 }

	elif species == 'PSD95' and type_analysis == 'distribution':
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
		for connections_from_a_molecule in connections_from_one_molecules:
			for j, ref in enumerate( reference_types_connection ):
				if utils.equal_list(connections_from_a_molecule, ref):
					dist[titles[j]] += 1
					
	elif species == 'PSD95' and type_analysis == 'ratio':
		ymax = 4
		types_connection = []
		for c in connections_from_one_molecules:
			if (c not in ['STG_PSD95']) and (c not in ['GluN2B_PSD95']):
				types_connection.append('None')
			elif (c in ['STG_PSD95']) and (c not in ['GluN2B_PSD95']):
				types_connection.append('STG only')
			elif (c not in ['STG_PSD95']) and (c in ['GluN2B_PSD95']):
				types_connection.append('PSD95 only')
			elif (c in ['STG_PSD95']) and (c in ['GluN2B_PSD95']):
				types_connection.append('Both')
		dist = {t: types_connection.count(t) for t in [ 'None', 'STG only', 'PSD95 only', 'Both' ]}
	
	return dist, ymax
	
	
if __name__ == '__main__':

	
	# Dataset 1
	species = 'GluN2B' # 'STG','GluN2B', 'PSD95','CaMKII', 'PSD95_connection'
	type_analysis = 'CaMKII'
	# 'average' and 'distribution' for all,
	# 'CaMKII' and 'PSD95' for GluN2B
	# 'ratio' for PSD95
	
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
	plt.suptitle('Species: {}, {}'.format( species, type_analysis ), )
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
			# sys.exit(0)
			
			dist, ymax = get_connection_statistics(multi_graph, species, type_analysis)
			
			title = prefix # prefix, None
			row    = num_rows-j-1
			column = i+1
			print('column: ', column, ', row: ', row)
			ax = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
			ax.set_title(prefix)
			ax.bar(*zip(*dist.items()), width=0.5)
			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)
			ax.set_ylim(0,ymax)
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
	
	
	