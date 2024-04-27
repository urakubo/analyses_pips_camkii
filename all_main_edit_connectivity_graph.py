
import os, sys, glob, pickle, pprint
import numpy as np
import math


import skimage
from scipy import ndimage

import matplotlib
import matplotlib.pyplot as plt

import networkx as nx


import utils
import parameters as p
import colormap as c


	
def get_surface_condensate_CaMKII(types, positions, ids_molecule, sigma=2):
	
	# Parameters
	
	positions_in_grid_coord = np.floor(positions + p.center_np).astype('int')
	locs_in_grid_mesh = \
		{v['id']:  np.ones([len(p.edge0),len(p.edge1),len(p.edge2)], dtype='int')*(-1) for v in p.subunits.values()}
	
	for i, (t, pos) in enumerate( zip(types, positions_in_grid_coord) ):
		locs_in_grid_mesh[t][pos[0],pos[1],pos[2]] = i
	
	# Define the region periphery.
	region_periphery = utils.get_periphery_in_grid_mesh()
	
	# Get the CaMKII condenate regions.
	conc_CaMKII_in_grid_mesh = [ locs_in_grid_mesh[id] >= 0 for id in p.molecules_without_all['CaMKII']['id'] ]
	conc_CaMKII_in_grid_mesh = sum(conc_CaMKII_in_grid_mesh).astype('float')
	conc_CaMKII_in_grid_mesh = ndimage.gaussian_filter(conc_CaMKII_in_grid_mesh, sigma = sigma)
	condensate_CaMKII_in_grid_mesh = utils.get_high(conc_CaMKII_in_grid_mesh)


	# Get the condensate of all.
	conc_all_in_grid_mesh = [ v >= 0 for v in locs_in_grid_mesh.values()]
	conc_all_in_grid_mesh = sum( conc_all_in_grid_mesh ).astype('float')
	
	# print('conc_all_in_grid_mesh.shape : ', conc_all_in_grid_mesh.shape ) 
	conc_all_in_grid_mesh = ndimage.gaussian_filter(conc_all_in_grid_mesh, sigma = sigma)
	conc_all_periphery    = np.sum(conc_all_in_grid_mesh * region_periphery ) / np.sum( region_periphery )
	condensate_all_in_grid_mesh = utils.get_high(conc_all_in_grid_mesh-conc_all_periphery)


	# Get the interface region
	footprint_inside  = skimage.morphology.ball(1) #### "2": r = 3, "1": r = 2
	footprint_outside = skimage.morphology.ball(15) #### "2": r = 3, "1": r = 2
	dil1 = skimage.morphology.binary_dilation( condensate_CaMKII_in_grid_mesh, footprint=footprint_outside )
	ero1 = skimage.morphology.binary_erosion( condensate_CaMKII_in_grid_mesh, footprint=footprint_inside )
	region_interface_CaMKII = np.logical_xor(dil1, ero1)
	
	dil2 = skimage.morphology.binary_dilation( condensate_all_in_grid_mesh, footprint=footprint_outside  )
	ero2 = skimage.morphology.binary_erosion( condensate_all_in_grid_mesh , footprint=footprint_inside )
	region_interface_all = np.logical_xor(dil2, ero2)
	
	beads_interface_CaMKII = {v['id']: locs_in_grid_mesh[v['id']][region_interface_CaMKII] for v in p.subunits.values() } 
	beads_interface_CaMKII = {k: v[v != -1] for k, v in beads_interface_CaMKII.items() } 
	beads_interface_all    = {v['id']: locs_in_grid_mesh[v['id']][region_interface_all] for v in p.subunits.values() } 
	beads_interface_all    = {k: v[v != -1] for k, v in beads_interface_all.items() } 
	
	# Summary
	d = {
			'locs_in_grid_mesh'				: locs_in_grid_mesh,
			'conc_CaMKII_in_grid_mesh'		: conc_CaMKII_in_grid_mesh,
			'condensate_CaMKII_in_grid_mesh': condensate_CaMKII_in_grid_mesh,
			'conc_all_in_grid_mesh'			: conc_all_in_grid_mesh,
			'condensate_all_in_grid_mesh'	: condensate_all_in_grid_mesh,
			'region_interface_CaMKII'		: region_interface_CaMKII,
			'region_interface_all'			: region_interface_all,
			'beads_interface_CaMKII'		: beads_interface_CaMKII,
			'beads_interface_all'			: beads_interface_all,
		}
	return d
	
	
def get_graphs(ids_molecule, types, bp, positions_grid_coord):
	
	unique_ids_molecule = np.unique( ids_molecule )
	multi_graph = nx.MultiGraph()
	simple_graph_CaMKII_GluN2B = nx.Graph()
	# CaMKII_or_GluN2B = np.zeros( len(ids_molecule), dtype = 'int')
	
	# Create nodes
	for id in unique_ids_molecule:
		flags = (id == ids_molecule)
		ids_bead         = np.nonzero(flags)
		types_bead       = types[flags]
		positions_bead   = positions_grid_coord[flags,:]
		
		ids_partner_bead = bp[flags]
		species = [k for k, i in p.molecules_without_all.items() if types_bead[0] in i['id']][0]
		
		
		multi_graph.add_node(id,\
			species    = species, \
			ids_bead   = ids_bead, \
			types_bead = types_bead, \
			num_beads  = np.sum(flags), \
			positions_grid_coord = positions_bead, \
			ids_partner_bead =ids_partner_bead)
		
		if species in ['CaMKII','GluN2B']:
			simple_graph_CaMKII_GluN2B.add_node(id, \
				species    = species, \
				ids_bead   = ids_bead, \
				types_bead = types_bead, \
				num_beads  = np.sum(flags), \
				positions_grid_coord = positions_bead, \
				ids_partner_bead =ids_partner_bead)
		
	# Make connection
	list_ids = list(multi_graph.nodes.keys())
	ids_bead_already_connected = np.zeros_like( bp, dtype='int' )
	for i in range(ids_molecule.shape[0]):
		if 	(bp[i] >= 0) and \
			(ids_bead_already_connected[i] == 0) and \
			(ids_bead_already_connected[bp[i]] == 0):
			
			ids_bead_already_connected[i]     = 1
			ids_bead_already_connected[bp[i]] = 1
			
			id_molecule = ids_molecule[i]
			id_molecule_partner = ids_molecule[bp[i]]
			
			species = multi_graph.nodes[id_molecule]['species']
			species_partner = multi_graph.nodes[id_molecule_partner]['species']
			connecting_species = [species, species_partner]
			if ('GluN2B' in connecting_species) and ('CaMKII' in connecting_species):
				type_connection = 'GluN2B_CaMKII'
			elif ('GluN2B' in connecting_species) and ('PSD95' in connecting_species):
				type_connection = 'GluN2B_PSD95'
			elif ('STG' in connecting_species) and ('PSD95' in connecting_species):
				type_connection = 'STG_PSD95'
			else:
				raise ValueError("Erronous connection: {}", connecting_species)
			multi_graph.add_edge(id_molecule, id_molecule_partner, type_connection = type_connection, id_bead1 = i ,id_bead2 = bp[i])
			
			if type_connection == 'GluN2B_CaMKII':
				simple_graph_CaMKII_GluN2B.add_edge(id_molecule, id_molecule_partner, type_connection = type_connection, id_bead1 = i ,id_bead2 = bp[i])
	
	return multi_graph, simple_graph_CaMKII_GluN2B
	
	
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
		dist = {'GluN2B_CaMKII': num_GluN2B_CaMKII,
				'GluN2B_PSD95' : num_GluN2B_PSD95 }
		
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
		dist = {'STG_PSD95': num_STG_PSD95, \
				'GluN2B_PSD95': num_GluN2B_PSD95 }

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
		#print('connections_from_one_molecules')
		#print(connections_from_one_molecules)
		for c in connections_from_one_molecules:
			if ('STG_PSD95' not in c) and ('GluN2B_PSD95' not in c):
				types_connection.append('None')
			elif ('STG_PSD95' in c) and ('GluN2B_PSD95' not in c):
				types_connection.append('STG only')
			elif ('STG_PSD95' not in c) and ('GluN2B_PSD95' in c):
				types_connection.append('PSD95 only')
			elif ('STG_PSD95' in c) and ('GluN2B_PSD95' in c):
				types_connection.append('Both')
		dist = {t: types_connection.count(t) for t in [ 'None', 'STG only', 'PSD95 only', 'Both' ]}
	
	return dist
	
	
if __name__ == '__main__':
	
	
	# CG Valency length
	'''
	subdirs    = ['val{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
	filenames_lammpstrj = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_edited    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
	dir_target          = 'CG_valency_length'
	'''
	
	
	# Small colony 
	subdirs    = ['CaMKII_432_GluN2Bc_8640', 'CaMKII_864_GluN2Bc_8640']
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
	filenames_lammpstrj = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_edited    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(3) for id_f in range(7)]
	dir_target  = 'small_colony'
	
	
	# Shared part of initialization
	dir_lammpstrj    = os.path.join('..', 'lammpstrj4', dir_target)
	dir_edited_data  = os.path.join('data4', dir_target)
	os.makedirs(dir_edited_data, exist_ok=True)
	suffix = 'connectivity_graph'
	
	#
	for filename_lammpstrj, filename_edited in zip(filenames_lammpstrj, filenames_edited):
		
		# Check sampling point
		print("\n"+filename_lammpstrj)
		sampling_frame = utils.get_num_frames(dir_lammpstrj, filename_lammpstrj)
		# print("The last timeframe was sampled: ", sampling_frame )
		
		
		# Load data
		types, positions_,ids_molecule, mc_step = \
			utils.load_lammpstrj( dir_lammpstrj, filename_lammpstrj, sampling_frame )
		bp = utils.load_lammpstrj_binding_partners( dir_lammpstrj, filename_lammpstrj, sampling_frame )
		
		
		# Centering
		center_of_mass = utils.get_center_of_mass(types, positions_) # reference_molecule_for_centering = 'CaMKII'
		positions = utils.centering(positions_, center_of_mass)
		
		
		# Generate graph
		multi_graph, simple_graph_CaMKII_GluN2B = get_graphs(ids_molecule, types, bp, positions)
		d = {}
		d['multi_graph'] = multi_graph
		d['simple_graph_CaMKII_GluN2B'] = simple_graph_CaMKII_GluN2B
		
		for species in ['STG','GluN2B', 'PSD95','CaMKII']:
			d[species] = {}
			for type_analysis in ['average', 'distribution']:
				d[species][type_analysis] = get_connection_statistics(multi_graph, species, type_analysis)
		
		species = 'GluN2B'
		for type_analysis in ['CaMKII', 'PSD95']:
			d[species][type_analysis] = get_connection_statistics(multi_graph, species, type_analysis)
		
		species = 'PSD95'
		type_analysis = 'ratio'
		d[species][type_analysis] = get_connection_statistics(multi_graph, species, type_analysis)
		
		
		d['condensate_CaMKII'] = get_surface_condensate_CaMKII(types, positions, ids_molecule)
		
		# Save the edited data
		prefix = filename_edited
		utils.save(dir_edited_data, prefix, suffix, d)
		
