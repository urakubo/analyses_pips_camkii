
import os, sys, glob, pickle, pprint
import numpy as np

import utils
	
	
	#
	# grid_coord: Molecular locations in the grid space (0, 1,..., 119) 
	# (numpy uint) [(x0,y0,z0), (x1,y1,z1), ..., (xn,yn,zn)]
	#
	# real_coord: Molecular locations in the real space [-60, 60) 
	# (numpy float) [(x0,y0,z0), (x1,y1,z1), ..., (xn,yn,zn)]
	#
	# grid_mesh : molecular existance in a 3D numpy variable
	# (True/False in the 3D space (util.space[0], util.space[1], util.space[2]))
	#
	
	
	

if __name__ == '__main__':
	
	
	# Dataset
	dir_lammpstrj    = './../lammpstrj/figure2_data'
	dir_edited_data  = 'data/conc_dependence'
	filenames_output = [str(i).zfill(3) for i in range(70) ]
	filenames_input  = ['R1_{}.lammpstrj'.format(f) for f in filenames_output ] #70
	
	# RDF
	num_sampling_frames = 5
	sampling_interval   = 2
	
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
		
		# RDF
		rdf, rdf_bins, rdf_sampling_frames = \
			utils.get_rdfs( dir_lammpstrj, filename_input, sampling_frame )
		
		sigmas = [2]
		for sigma in sigmas:
			print('sigma ', sigma)
			# Get concs
			d = utils.get_concs_and_condensates(types, positions_real_coord, ids_molecule, sigma)
			# Watershed
			labels_watershed_in_grid_mesh, ratio_volumes_watershed = utils.watershed_segmentation( d )
			d['labels_watershed_in_grid_mesh'] = labels_watershed_in_grid_mesh
			d['ratio_volumes_watershed']       = ratio_volumes_watershed
			
			# RDF
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

