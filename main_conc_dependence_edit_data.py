
import os, sys, glob, pickle, pprint
import numpy as np

import utils
	
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
	
	
if __name__ == '__main__':
	
	
	# Dataset
	dir_lammpstrj    = os.path.join('..', 'lammpstrj','figure2_data')
	dir_edited_data  = os.path.join('data', 'conc_dependence')
	filenames_output = [str(i).zfill(3) for i in range(70) ]
	filenames_input  = ['R1_{}.lammpstrj'.format(f) for f in filenames_output ] #70
	
	
	# Init
	os.makedirs(dir_edited_data, exist_ok=True)
	
	for filename_input, filename_output in zip(filenames_input, filenames_output):
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
		
