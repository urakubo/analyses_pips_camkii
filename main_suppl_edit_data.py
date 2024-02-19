
import os, sys, glob, pickle, pprint
import numpy as np

import utils
	
	
if __name__ == '__main__':
	
	'''
	# Dataset 1
	filenames_lammpstrj = ['CG.lammpstrj',\
				'SP.lammpstrj',\
				'SPG1.lammpstrj',\
				'SPG2.lammpstrj']
	filenames_edited_data = ['CG',\
				'SP',\
				'SPG1',\
				'SPG2']
	target_dir = 'mix_two_three_components'
	
	# Dataset 2
	filenames_lammpstrj =  ['OnlyCaMKIICaMKIIAffinity.lammpstrj',\
							'OnlyGluN2BGluN2B.lammpstrj',\
							'onlySTGSTG.lammpstrj']
	filenames_edited_data = ['CaMKIIalone',\
							'GluN2Balone',\
							'STGalone']
	target_dir = 'self_affinity'
	'''
	
	# Dataset 2
	filenames_lammpstrj =  ['Activated CaMKII_2.lammpstrj',\
							'Unactivated CamKII_0.lammpstrj',\
							'Unactivated CamKII_1.lammpstrj']
	filenames_edited_data = ['activated_CaMKII_2',\
							'unactivated_CamKII_0',\
							'unactivated_CamKII_1']
	target_dir = 'figure_1c'

	#
	dir_lammpstrj    = os.path.join('..', 'lammpstrj', target_dir)
	dir_edited_data  = os.path.join('data',target_dir)
	
	
	for filename_input, filename_output in zip(filenames_lammpstrj, filenames_edited_data):
		
		
		
		# Load data
		print("\n"+filename_input)
		sampling_frame = utils.get_num_frames(dir_lammpstrj, filename_input)
		
		print("sampling_frame ", sampling_frame )
		types, positions_grid_coord,ids_molecule, mc_step = \
			utils.load_lammpstrj( dir_lammpstrj, filename_input, sampling_frame )
		print("The last timeframe was loaded." )
		
		# Centering
		center    = utils.get_center_of_mass(types, positions_grid_coord)
		positions_real_coord = utils.centering(positions_grid_coord, center)
		
		# Blur
		sigma = 2
		print('sigma ', sigma)
		# Get concs and condensates
		d = utils.get_concs_and_condensates(types, positions_real_coord, ids_molecule, sigma)
		
		# Watershed
		'''
		labels_watershed_in_grid_mesh, ratio_volumes_watershed = utils.watershed_segmentation( d )
		d['labels_watershed_in_grid_mesh'] = labels_watershed_in_grid_mesh
		d['ratio_volumes_watershed']       = ratio_volumes_watershed
		'''
		
		
		# RDF
		rdf, rdf_bins, rdf_sampling_frames = \
			utils.get_rdfs( dir_lammpstrj, filename_input, sampling_frame )
		d['rdf_bins'] = rdf_bins
		d['rdf_sampling_frames'] = rdf_sampling_frames
		d['rdf'] = rdf
		
		
		# Time info
		d['time_frame'] = sampling_frame
		d['mc_step']    = mc_step
		
		
		# Save the edited data
		prefix = filename_output
		suffix = 'sigma_{}'.format(sigma)
		utils.save(dir_edited_data, prefix, suffix, d)
		
		
		
