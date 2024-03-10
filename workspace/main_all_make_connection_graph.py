
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
	
	
	#
	# grid_coord: Molecular coordinate in the grid space (0, 1,..., 119) 
	# (numpy uint) [(x0,y0,z0), (x1,y1,z1), ..., (xn,yn,zn)]
	#
	# real_coord: Molecular coordinate in the real space [-60, 60) 
	# (numpy float) [(x0,y0,z0), (x1,y1,z1), ..., (xn,yn,zn)]
	#
	# grid_mesh : molecular existance in the 3D space
	# (numpy bool in 3D space) (p.space[0], p.space[1], p.space[2])
	#
	
	
def get_clusters(ids_molecule, types, bp):
	
	unique_ids_molecule = np.unique( ids_molecule )
	
	# Create a graph and nodes
	graph = nx.Graph()
	for id in unique_ids_molecule:
		flags = (id == ids_molecule)
		ids_bead = np.nonzero(flags)
		types_bead = types[flags]
		ids_partner_bead = bp[flags]
		molecular_species = [k for k, i in p.molecules_without_all.items() if types_bead[0] in i['id']]
		graph.add_node(id, molecular_species = molecular_species[0], \
			ids_bead   = ids_bead, \
			types_bead = types_bead, \
			num_beads  = np.sum(flags), \
			ids_partner_bead =ids_partner_bead)
	
	# Make connection
	for i in range(ids_molecule.shape[0]):
		if bp[i] >= 0:
			id_molecule = ids_molecule[i]
			# print('bp[i] ', bp[i])
			id_molecule_partner = ids_molecule[bp[i]]
			graph.add_edge(id_molecule, id_molecule_partner)
	return graph
	
	
if __name__ == '__main__':
	
	#'''
	# Dataset 2:  Valency length
	subdirs    = ['val{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(5)]
	filenames_input  = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_output = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(5) ]
	dir_input        = 'Feb_Figure3'
	dir_edited_data  = 'valency_linker_length'
	#'''
	
	# Shared part of initialization
	dir_lammpstrj    = os.path.join('..','..', 'lammpstrj2', dir_input)
	dir_edited_data  = os.path.join('..','data2',dir_edited_data)
	os.makedirs(dir_edited_data, exist_ok=True)
	suffix = 'graph'
	
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
		graph = get_clusters(ids_molecule, types, bp)
		
		# Save the edited data
		prefix = filename_output
		utils.save(dir_edited_data, prefix, suffix, graph)
		
