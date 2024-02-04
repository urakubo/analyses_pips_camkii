
import os, sys, glob, pickle, pprint
import numpy as np
from scipy import ndimage

from skimage.segmentation import watershed
import utils
	
	
	
def get_rdfs( dir_input, filename_input, target_time_frame ):
	
	# Parameters
	rdf_bins          = np.arange(0, 35) # You can set "np.arange(0, 35, 2)"
	rdf_grid_points   = utils.get_lattice_grids()
	num_frames        = 5
	sampling_interval = 2
	
	# Target frames
	if target_time_frame == 0:
		rdf_target_frames = [0]
	else:
		rdf_target_frames = list(range( target_time_frame - (num_frames - 1)*sampling_interval, \
					target_time_frame + sampling_interval,\
					sampling_interval))
	print('rdf_target_frames ', rdf_target_frames)
	rdf = utils.get_rdfs_from_multiple_frames(dir_input, filename_input, \
			rdf_target_frames, rdf_bins, rdf_grid_points)
	
	return rdf, rdf_bins, rdf_target_frames
	#
def get_average_center_of_mass( dir_lammpstrj, filename_input, sampling_time_frames ):
	centers = np.zeros( (len(sampling_time_frames), 3) )
	for i, target_frame in enumerate( sampling_time_frames ):
		types, positions, _ = utils.load_data( dir_lammpstrj, filename_input, target_frame )
		centers[i,:] = utils.get_center_of_mass(types, positions)
	center = np.mean(centers, axis = 0)
	return center

if __name__ == '__main__':
	
	# Input files
	dir_lammpstrj   = './../lammpstrj/mix_two_three_components'
	filenames_lammpstrj   = [		'CG.lammpstrj',\
							'SP.lammpstrj',\
							'SPG1.lammpstrj',\
							'SPG2.lammpstrj']
	
	# Output files
	dir_edited_data = 'data'
	filenames_edited_data = ['CG',\
							'SP',\
							'SPG1',\
							'SPG2']
	
	# Parameters
	sigma = 2
	sampling_time_frames = list(range(0,200,20))
	print('Sampling time frame: ', sampling_time_frames)
	
	
	for filename_input, filename_output in zip(filenames_lammpstrj, filenames_edited_data):
		
		print('filename_input ', filename_input)
		# Get centers-of-mass
		# center = get_average_center_of_mass( dir_lammpstrj, filename_input, sampling_time_frames )
		# The above is problematic if the centers are located at the edges of ends.
		center = utils.center # Or the above
		print(center)
		
		
		for target_frame in sampling_time_frames:
			
			print('target_time_frame ', target_frame)
			
			# Load data
			types, positions,ids_molecule = utils.load_data( dir_lammpstrj, filename_input, target_frame )
			
			# Centering
			positions = utils.centering(positions, center)
			
			# Get concs and condensates
			d = utils.get_concs_and_condensates(types, positions, ids_molecule, sigma)
			
			# RDF
			rdf, rdf_bins, rdf_target_frames = \
				get_rdfs( dir_lammpstrj, filename_input, target_frame )
			d['rdf_bins'] = rdf_bins
			d['rdf_target_frames'] = rdf_target_frames
			d['rdf'] = rdf
			
			# Save the edited data
			prefix = filename_output
			suffix = 'frame_{}'.format(target_frame)
			utils.save(dir_edited_data, prefix, suffix, d)

