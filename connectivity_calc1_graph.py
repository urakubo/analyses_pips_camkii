
import os, sys, glob, pickle, pprint
import numpy as np

import skimage
from scipy import ndimage

import networkx as nx

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph
	
	
	
def process_and_save2(types, positions_grid_coord, ids_molecule, bp, filename_edited, suffix, mc_step, sampling_frame ):

	# Centering
	# reference_molecule_for_centering = 'CaMKII'
	center_of_mass = utils.get_center_of_mass(types, positions_grid_coord)
	positions_real_coord = utils.centering(positions_grid_coord, center_of_mass)
	
	# Generate graph
	multi_graph = utils_graph.get_multi_graph(ids_molecule, types, bp, positions_real_coord)
	d = {}
	d['mc_step']        = mc_step
	d['sampling_frame'] = sampling_frame
	d['multi_graph']    = multi_graph
	#d['simple_graph_CaMKII_GluN2B'] = simple_graph_CaMKII_GluN2B
	
	# Save the edited data
	prefix = filename_edited
	utils.save(dir_edited_data, prefix, suffix, d)
	
	
def repeat_for_valency_length(filenames_lammpstrj, filenames_edited):
	#
	suffix = 'connectivity_graph'
	for filename_lammpstrj, filename_edited in zip(filenames_lammpstrj, filenames_edited):
		
		# Check sampling point
		sampling_frame = utils.get_num_frames(dir_lammpstrj, filename_lammpstrj)
		
		# Load data
		types, positions,ids_molecule, mc_step = \
			utils.load_lammpstrj( dir_lammpstrj, filename_lammpstrj, sampling_frame )
		bp = utils.load_lammpstrj_binding_partners( dir_lammpstrj, filename_lammpstrj, sampling_frame )
		
		print("\n"+filename_lammpstrj)
		# print("The last timeframe was sampled: ", sampling_frame )
		
		process_and_save2(types, positions, ids_molecule, bp, filename_edited, suffix, mc_step, sampling_frame )
		

def repeat_for_time_development(filename_lammpstrj, filename_edited, max_backward_frames_for_sampling, num_skip_frames_for_sampling = 10):
	
	num_frames = utils.get_num_frames(dir_lammpstrj, filename_lammpstrj)
	for i in range(0, max_backward_frames_for_sampling):
		# Check sampling point
		sampling_frame = num_frames + (i - max_backward_frames_for_sampling)*num_skip_frames_for_sampling
		suffix = 'connectivity_graph_{}'.format(i)
		
		# Load data
		types, positions,ids_molecule, mc_step = \
			utils.load_lammpstrj( dir_lammpstrj, filename_lammpstrj, sampling_frame )
		bp = utils.load_lammpstrj_binding_partners( dir_lammpstrj, filename_lammpstrj, sampling_frame )
		
		process_and_save2(types, positions, ids_molecule, bp, filename_edited, suffix, mc_step, sampling_frame )
		
		print("\n"+filename_lammpstrj)
		print("{} / {} ".format(i, max_backward_frames_for_sampling) )
		print("The sampling timeframe was: ", sampling_frame )
		print("The sampling MC step   was: ", mc_step )
		
	
	
if __name__ == '__main__':
	
	
	
	# Small colony 2
	#'''
	subdirs    = ['val_{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
	filenames_lammpstrj = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_edited    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7)]
	dir_target  = 'small_colony2'
	#'''
	
	# repeat_for_valency_length(filenames_lammpstrj, filenames_edited)
	
	i = 7*5+1 # val_12\R2_002
	i = 7*4+0 # val_10\R2_000
	i = 7*5+0 # val_12\R2_000
	i = 7*4+1 # val_10\R2_001
	i = 7*3+0 # val_08\R2_000
	
	
	# Small colony
	'''
	subdirs    = ['CaMKII_432_GluN2Bc_8640', 'CaMKII_864_GluN2Bc_8640']
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
	filenames_lammpstrj = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_edited    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(3) for id_f in range(7)]
	dir_target  = 'small_colony'
	i = 6 # 5
	num_skip_frames_for_sampling = 1
	i = 5 #
	num_skip_frames_for_sampling = 3
	'''
	
	# Shared part of initialization
	dir_lammpstrj    = os.path.join('..', 'lammpstrj4', dir_target)
	dir_edited_data  = os.path.join('data4', dir_target)
	os.makedirs(dir_edited_data, exist_ok=True)
	

	#'''
	filename_lammpstrj = filenames_lammpstrj[i]
	filename_edited    = filenames_edited[i]
	
	print('filename_lammpstrj ', filename_lammpstrj)
	print('filename_edited    ', filename_edited   )
	max_backward_frames_for_sampling = 80
	repeat_for_time_development(filename_lammpstrj, \
		filename_edited, \
		max_backward_frames_for_sampling, \
		num_skip_frames_for_sampling = num_skip_frames_for_sampling)	
	#'''
	
	'''
	for filename_lammpstrj, filename_edited in zip(filenames_lammpstrj, filenames_edited):
		print('filename_lammpstrj ', filename_lammpstrj)
		print('filename_edited    ', filename_edited   )
		max_backward_frames_for_sampling = 80
		num_skip_frames_for_sampling = 10
		repeat_for_time_development(filename_lammpstrj, \
			filename_edited, \
			max_backward_frames_for_sampling, \
			num_skip_frames_for_sampling = num_skip_frames_for_sampling)
	'''
	
	# 

