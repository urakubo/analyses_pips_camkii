
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
	
	
def get_clusters(ids_molecule, types, bp):
	
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
	
	
if __name__ == '__main__':
	
	
	# Valency length
	#'''
	subdirs    = ['val_{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
	filenames_input  = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_output = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
	dir_input        = 'valency_length'
	dir_edited_data  = 'valency_length'
	#'''
	
	# Conc dependnece
	'''
	filenames_output = [str(i).zfill(3) for i in range(48) ]
	filenames_input  = ['R2_{}.lammpstrj'.format(f) for f in filenames_output ] #70
	dir_input        = 'conc_dependence'
	dir_edited_data  = 'conc_dependence'
	'''
	
	# Shared part of initialization
	dir_lammpstrj    = os.path.join('..', 'lammpstrj3', dir_input)
	dir_edited_data  = os.path.join('data3',dir_edited_data)
	os.makedirs(dir_edited_data, exist_ok=True)
	suffix = 'graph2'
	
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
		g = get_clusters(ids_molecule, types, bp)
		
		
		# Save the edited data
		prefix = filename_output
		utils.save(dir_edited_data, prefix, suffix, g)
		
