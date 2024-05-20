
import os, sys, glob, pickle, pprint
import numpy as np
import math



import matplotlib
import matplotlib.pyplot as plt

import networkx as nx

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph
	
	
if __name__ == '__main__':
	

	# Conc dependnece
	'''
	filenames_output = [str(i).zfill(3) for i in range(81) ]
	filenames_input  = ['R2_{}.lammpstrj'.format(f) for f in filenames_output ]
	dir_input        = 'conc_dependence'
	dir_edited_data  = 'conc_dependence'
	'''
	'''
	filenames_output = [str(i).zfill(3) for i in range(9) ]
	filenames_input  = ['R2_{}.lammpstrj'.format(f) for f in filenames_output ]
	dir_input        = os.path.join('conc_dependence', '0.33')
	dir_edited_data  = 'conc_dependence_0.33'
	'''
	
	# Valency length
	#'''
	subdirs    = ['val_{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
	filenames_input  = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_output = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
	dir_input        = 'valency_length'
	dir_edited_data  = 'valency_length'
	#'''
	
	
	# CG Valency length
	'''
	subdirs    = ['val{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
	filenames_input  = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_output = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
	dir_input        = 'CG_valency_length'
	dir_edited_data  = 'CG_valency_length'
	'''
	
	
	# Shared part of initialization
	dir_lammpstrj    = os.path.join('..', 'lammpstrj4', dir_input)
	dir_edited_data  = os.path.join('data4',dir_edited_data)
	os.makedirs(dir_edited_data, exist_ok=True)
	suffix = 'connectivity_graph'
	
	#
	for filename_input, filename_output in zip(filenames_input, filenames_output):
		
		# Check sampling point
		print("\n"+filename_input)
		sampling_frame = utils.get_num_frames(dir_lammpstrj, filename_input)
		# print("The last timeframe was sampled: ", sampling_frame )
		
		
		# Load data
		types, positions_,ids_molecule, mc_step = \
			utils.load_lammpstrj( dir_lammpstrj, filename_input, sampling_frame )
		bp = utils.load_lammpstrj_binding_partners( dir_lammpstrj, filename_input, sampling_frame )
		
		
		# Centering
		center_of_mass = utils.get_center_of_mass(types, positions_) # reference_molecule_for_centering = 'CaMKII'
		positions_real_coord = utils.centering(positions_, center_of_mass)
		
		
		# Generate graph
		multi_graph = utils_graph.get_multi_graph(ids_molecule, types, bp, positions_real_coord)
		d = {}
		d['multi_graph']   = multi_graph
		d['dir_lammpstrj'] = dir_lammpstrj
		d['filename']      = filename_input
		d['sampling_frame']= sampling_frame
		
		for species in ['STG','GluN2B', 'PSD95','CaMKII']:
			d[species] = {}
			for type_analysis in ['average', 'distribution']:
				d[species][type_analysis] = utils_graph.get_connection_statistics(multi_graph, species, type_analysis)
		
		species = 'GluN2B'
		for type_analysis in ['CaMKII', 'PSD95']:
			d[species][type_analysis] = utils_graph.get_connection_statistics(multi_graph, species, type_analysis)
		
		species = 'PSD95'
		type_analysis = 'ratio'
		d[species][type_analysis] = utils_graph.get_connection_statistics(multi_graph, species, type_analysis)
		
		
		_, ids_bead_shared_PSD = utils.get_ids_PSD95_shared_by_STG_GluN2B(multi_graph, shared_or_unshared = 'shared')
		d['ids_bead_shared_PSD'] = ids_bead_shared_PSD
		
		
		d['condensate_CaMKII'] = utils_graph.get_surface_condensate_CaMKII(types, positions_real_coord, ids_molecule)
		
		# Save the edited data
		prefix = filename_output
		utils.save(dir_edited_data, prefix, suffix, d)
		
