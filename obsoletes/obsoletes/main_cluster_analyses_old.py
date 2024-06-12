
import os, sys, glob, pickle, pprint
import numpy as np
import math

import matplotlib
import matplotlib.pyplot as plt

import networkx as nx
from networkx.utils import pairwise

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
	
	
def get_clusters(positions_grid_coord, ids_molecule, types, energy_aniso, ids_grid_space):
	
	
	unique_ids_molecule = np.unique( ids_molecule )
	possible_types_connection = (\
		{ p.subunits['GluN2B binding site']['id'], p.subunits['CaMKII binding site']['id'] },\
		{ p.subunits['STG binding site']['id'], p.subunits['PSD']['id'] },\
		{ p.subunits['STG binding site']['id'], p.subunits['GluN2B binding site']['id']})
	print('possible_types_connection')
	print( possible_types_connection )
	# 
	
	# Create a graph and nodes
	graph = nx.Graph()
	for id in unique_ids_molecule:
		types_bead = types[id == ids_molecule]
		molecular_species = [k for k, i in p.molecules_without_all.items() if types_bead[0] in i['id']]
		graph.add_node(id, molecular_species = molecular_species[0], types_bead = types_bead)
	
	# Make connection
	for id in unique_ids_molecule:
		is_target_molecule = (id == ids_molecule)
		targ_types     = types[is_target_molecule]
		targ_energies  = energy_aniso[is_target_molecule]
		targ_positions = positions_grid_coord[is_target_molecule,:]
		
		has_energy     = (targ_energies < 0)
		targ_types     = targ_types[has_energy]
		targ_energies  = targ_energies[has_energy]
		targ_positions = targ_positions[has_energy,:]
		
		for ii in range( targ_positions.shape[0] ):
			counter_jj = 0
			for jj in range( p.neighbors26.shape[0] ):
				loc = targ_positions[ii,:] + p.neighbors26[jj,:]
				has_same_energy = (targ_energies[ii] == ids_grid_space['energy'][loc[0],loc[1],loc[2]])
				
				snd_type_bead = targ_types[ii]
				rec_type_bead = ids_grid_space['type_bead'][loc[0],loc[1],loc[2]]
				types_bead = {snd_type_bead, rec_type_bead}
				is_possible_connection = types_bead in possible_types_connection
				if has_same_energy and is_possible_connection:
					graph.add_edge(id, ids_grid_space['id_molecule'][loc[0],loc[1],loc[2]])
					print(counter_jj, types_bead, end=' , ')
					counter_jj += 1
	
	return graph
	
	
if __name__ == '__main__':

	i = 2
	
	# Dataset 1
	dir_input       = 'Feb_Figure2'
	fnames = [27,21,15,9]
	filenames_input  = [ 'R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in fnames ]
	filenames_output = [ '{}'.format(str(i).zfill(3)) for i in fnames ]

	# Shared part of initialization
	dir_lammpstrj    = os.path.join('..','..', 'lammpstrj2', dir_input)
	dir_edited_data  = os.path.join('data2','cluster_analyses')
	os.makedirs(dir_edited_data, exist_ok=True)
	
	#
	for filename_input, filename_output in zip(filenames_input, filenames_output):
		
		
		# Check sampling point
		print("\n"+filename_input)
		sampling_frame = utils.get_num_frames(dir_lammpstrj, filename_input)
		print("The last timeframe was sampled: ", sampling_frame )
		
		
		# Load data
		types, positions_grid_coord,ids_molecule, mc_step = \
			utils.load_lammpstrj( dir_lammpstrj, filename_input, sampling_frame )
		energy = utils.load_lammpstrj_binding_energy( dir_lammpstrj, filename_input, sampling_frame )
		energy_aniso = energy['energy_anisotropic']
		
		
		# Centering
		center = utils.get_center_of_mass(types, positions_grid_coord)
		center = np.floor(center).astype('int')
		positions_grid_coord = utils.centering_lattice_space(positions_grid_coord, center).astype('int')
		
		
		# Generate the grid space for ids.
		# Assign ids into the grid space
		'''
		ids_grid_space = {}
		ids_grid_space['id_bead'] = np.zeros(p.space, dtype='int')
		ids_grid_space['id_bead'][:,:] = np.iinfo('int').min
		ids_grid_space['type_bead'] = np.zeros(p.space, dtype='int')
		ids_grid_space['type_bead'][:,:] = np.iinfo('int').min
		ids_grid_space['id_molecule'] = np.zeros(p.space, dtype='int')
		ids_grid_space['id_molecule'][:,:] = np.iinfo('int').min
		ids_grid_space['energy']   = np.zeros(p.space, dtype='float')
		ids_grid_space['energy'][:,:] = np.nan
		for i in range(positions_grid_coord.shape[0]):
			ids_grid_space['id_bead'][\
				positions_grid_coord[i,0],
				positions_grid_coord[i,1],
				positions_grid_coord[i,2]] = i
			ids_grid_space['type_bead'][\
				positions_grid_coord[i,0],
				positions_grid_coord[i,1],
				positions_grid_coord[i,2]] = types[i]
			ids_grid_space['id_molecule'][\
				positions_grid_coord[i,0],
				positions_grid_coord[i,1],
				positions_grid_coord[i,2]] = ids_molecule[i]
			ids_grid_space['energy'][\
				positions_grid_coord[i,0],
				positions_grid_coord[i,1],
				positions_grid_coord[i,2]] = energy_aniso[i]
		graph = get_clusters(positions_grid_coord, ids_molecule, types, energy_aniso, ids_grid_space)
		'''
		
		
		
		# Save the edited data
		prefix = filename_output
		utils.save(dir_edited_data, prefix, suffix, graph)
		
		
