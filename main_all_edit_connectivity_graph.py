
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
	
	
def get_graphs(ids_molecule, types, bp):
	
	unique_ids_molecule = np.unique( ids_molecule )
	graph = nx.Graph()
	multi_graph = nx.MultiGraph()
	CaMKII_or_GluN2B = np.zeros( len(ids_molecule), dtype = 'int')
	
	# Create nodes
	for id in unique_ids_molecule:
		flags = (id == ids_molecule)
		ids_bead         = np.nonzero(flags)
		types_bead       = types[flags]
		ids_partner_bead = bp[flags]
		species = [k for k, i in p.molecules_without_all.items() if types_bead[0] in i['id']][0]
		
		graph.add_node(id, \
			species    = species, \
			ids_bead   = ids_bead, \
			types_bead = types_bead, \
			num_beads  = np.sum(flags), \
			ids_partner_bead =ids_partner_bead)
		
		multi_graph.add_node(id,\
			species    = species, \
			ids_bead   = ids_bead, \
			types_bead = types_bead, \
			num_beads  = np.sum(flags), \
			ids_partner_bead =ids_partner_bead)
		
	# Make connection
	list_ids = list(graph.nodes.keys())
	ids_bead_already_connected = np.zeros_like( bp, dtype='int' )
	for i in range(ids_molecule.shape[0]):
		if 	(bp[i] >= 0) and \
			(ids_bead_already_connected[i] == 0) and \
			(ids_bead_already_connected[bp[i]] == 0):
			
			ids_bead_already_connected[i]     = 1
			ids_bead_already_connected[bp[i]] = 1
			
			id_molecule = ids_molecule[i]
			id_molecule_partner = ids_molecule[bp[i]]
			
			species = graph.nodes[id_molecule]['species']
			species_partner    = graph.nodes[id_molecule_partner]['species']
			connecting_species = [species, species_partner]
			if ('GluN2B' in connecting_species) and ('CaMKII' in connecting_species):
				type_connection = 'GluN2B_CaMKII'
			elif ('GluN2B' in connecting_species) and ('PSD95' in connecting_species):
				type_connection = 'GluN2B_PSD95'
			elif ('STG' in connecting_species) and ('PSD95' in connecting_species):
				type_connection = 'STG_PSD95'
			else:
				raise ValueError("Erronous connection: {}", connecting_species)
			
			graph.add_edge(id_molecule, id_molecule_partner, type_connection = type_connection)
			multi_graph.add_edge(id_molecule, id_molecule_partner, type_connection = type_connection)
	
	return graph, multi_graph
	
	
def get_connection_statistics(multi_graph, species, type_analysis):
	
	ids = [ i for i, attr in multi_graph.nodes('species') if attr == species ]
	numbers_of_connections = [ multi_graph.degree[id] for id in ids ]
	
	edges_from_molecules = [ multi_graph.edges(id, keys=True, data=True) for id in ids ]
	connections_from_one_molecules = [ [e[3]['type_connection'] for e in es] for es in edges_from_molecules ]
	
	
	if species == 'GluN2B' and type_analysis in ['CaMKII','PSD95']:
		types_connection = [[len(v) for k, v in multi_graph[id].items() if multi_graph.nodes[k]['species'] == type_analysis] for id in ids]
		reference_types_connection = {\
			'No {}'.format(type_analysis) :[],\
			'One {} (single)'.format(type_analysis) :[1],\
			'One {} (double)'.format(type_analysis) :[2],\
			'Two {}s'.format(type_analysis) :[1, 1]\
			}
		dist = {k: types_connection.count(v) for k, v in reference_types_connection.items()}
		
	elif species in ['CaMKII','STG'] and type_analysis == 'distribution':
		partner      = {'CaMKII':'GluN2B', 'STG': 'PSD95'}
		binding_num  = {'CaMKII':13, 'STG':5}
		dist = {'{:d} {}'.format(i, partner[species]): numbers_of_connections.count(i) for i in range(binding_num[species])}
		
	elif species in ['CaMKII','STG'] and type_analysis == 'average':
		partner = {'CaMKII':'GluN2B', 'STG': 'PSD95'}
		ave_numbers_of_connections = np.average( numbers_of_connections )
		dist = {partner[species]: ave_numbers_of_connections}
		
	elif species == 'GluN2B' and type_analysis == 'average':
		num_GluN2B_CaMKII = np.average( [ c.count('GluN2B_CaMKII') for c in connections_from_one_molecules ] )
		num_GluN2B_PSD95  = np.average( [ c.count('GluN2B_PSD95') for c in connections_from_one_molecules ] )
		dist = {'GluN2B_CaMKII \n {:.2f}'.format(num_GluN2B_CaMKII): num_GluN2B_CaMKII,
				'GluN2B_PSD95  \n {:.2f}'.format(num_GluN2B_PSD95) : num_GluN2B_PSD95 }
		
	elif species == 'GluN2B' and type_analysis == 'distribution':
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
		num_STG_PSD95    = np.average( [c.count('STG_PSD95')  for c in connections_from_one_molecules ]   )
		num_GluN2B_PSD95 = np.average( [c.count('GluN2B_PSD95') for c in connections_from_one_molecules ] )
		dist = {'STG_PSD95 \n {:.2f}'.format(num_STG_PSD95): num_STG_PSD95, \
				'GluN2B_PSD95 \n {:.2f}'.format(num_GluN2B_PSD95): num_GluN2B_PSD95 }

	elif species == 'PSD95' and type_analysis == 'distribution':
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
	
	return dist
	
	
if __name__ == '__main__':
	
	
	# Valency length
	'''
	subdirs    = ['val_{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
	filenames_input  = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_output = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
	dir_input        = 'valency_length'
	dir_edited_data  = 'valency_length'
	'''
	
	# Conc dependnece
	#'''
	filenames_output = [str(i).zfill(3) for i in range(48) ]
	filenames_input  = ['R2_{}.lammpstrj'.format(f) for f in filenames_output ] #70
	dir_input        = 'conc_dependence'
	dir_edited_data  = 'conc_dependence'
	#'''
	
	# Shared part of initialization
	dir_lammpstrj    = os.path.join('..', 'lammpstrj3', dir_input)
	dir_edited_data  = os.path.join('data3',dir_edited_data)
	os.makedirs(dir_edited_data, exist_ok=True)
	suffix = 'connectivity_graph'
	
	#
	for filename_input, filename_output in zip(filenames_input, filenames_output):
		
		# Check sampling point
		print("\n"+filename_input)
		sampling_frame = utils.get_num_frames(dir_lammpstrj, filename_input)
		# print("The last timeframe was sampled: ", sampling_frame )
		
		
		# Load data
		types, positions_grid_coord,ids_molecule, mc_step = \
			utils.load_lammpstrj( dir_lammpstrj, filename_input, sampling_frame )
		bp = utils.load_lammpstrj_binding_partners( dir_lammpstrj, filename_input, sampling_frame )
		
		
		# Generate graph
		graph, multi_graph = get_graphs(ids_molecule, types, bp)
		d = {}
		d['graph'] = graph
		d['multi_graph'] = multi_graph
		
		for species in ['STG','GluN2B', 'PSD95','CaMKII']:
			for type_analysis in ['average', 'distribution']:
				d[species] = {}
				d[species][type_analysis] = get_connection_statistics(multi_graph, species, type_analysis)
		
		species = 'GluN2B'
		for type_analysis in ['CaMKII', 'PSD95']:
			d[species][type_analysis] = get_connection_statistics(multi_graph, species, type_analysis)
		
		species = 'PSD95'
		type_analysis = 'ratio'
		d[species][type_analysis] = get_connection_statistics(multi_graph, species, type_analysis)
		
		
		# Save the edited data
		prefix = filename_output
		utils.save(dir_edited_data, prefix, suffix, d)
		
