
import os, sys, glob, pickle, pprint
import numpy as np
from scipy import ndimage

from skimage.segmentation import watershed
import utils
	
	#
	# grid_coord: Molecular locations in the grid space [0, 1,..., 119]
	# (numpy uint) [(x0,y0,z0), (x1,y1,z1), ..., (xn,yn,zn)]
	#
	# real_coord: Molecular locations in the real space [-60, 60) 
	# (numpy float) [(x0,y0,z0), (x1,y1,z1), ..., (xn,yn,zn)]
	#
	# grid_mesh : molecular existance in a 3D numpy variable
	# (True/False in the 3D space (p.space[0], p.space[1], p.space[2]))
	#
	
def get_rdfs( dir_input, filename_input, target_frame, center ):
	
	# Parameters
	rdf_bins          = np.arange(0, 35) # You can set "np.arange(0, 35, 2)"
	rdf_grid_points   = utils.get_lattice_grids()
	num_sampling_frames = 5
	sampling_interval   = 2
	
	# Target frames
	rdf_sampling_frames = list(range( target_frame - (num_sampling_frames - 1)*sampling_interval, \
			target_frame + sampling_interval,\
			sampling_interval))
	
	if np.any(np.array(rdf_sampling_frames) < 0):
		rdf_sampling_frames = [0]

	print('rdf_sampling_frames ', rdf_sampling_frames)
	rdf = utils.get_rdfs_from_multiple_frames(dir_input, filename_input, \
			rdf_sampling_frames, rdf_bins, rdf_grid_points, center = center)
	
	return rdf, rdf_bins, rdf_sampling_frames
	#

	
if __name__ == '__main__':
	
	# Input files
	dir_lammpstrj   = './../lammpstrj/mix_two_three_components'
	filenames_lammpstrj = ['CG.lammpstrj',\
				'SP.lammpstrj',\
				'SPG1.lammpstrj',\
				'SPG2.lammpstrj']
	sampling_interval = 20
	sampling_time_frames_ = list(range(0,200+sampling_interval,sampling_interval))
	sampling_time_framess = 4*[sampling_time_frames_]
	
	# Output files
	dir_edited_data = 'data'
	filenames_edited_data = ['CG',\
				'SP',\
				'SPG1',\
				'SPG2']
	calc_center = 'each step' # 'none', 'each step', 'sampling time average'
	
	
	dir_lammpstrj   = './../lammpstrj/self_affinity'
	filenames_lammpstrj =  ['OnlyCaMKIICaMKIIAffinity.lammpstrj',\
							'OnlyGluN2BGluN2B.lammpstrj',\
							'onlySTGSTG.lammpstrj']
	dir_edited_data = 'data'
	filenames_edited_data = ['CaMKIIalone',\
							'GluN2Balone',\
							'STGalone']
	
	flag_calc_center = False
	sampling_time_framess = [ [0, 100, 200], [0, 100, 187] , [0, 100, 195] ]
	
	
	# Parameters
	sigma = 2
	
	
	for filename_input, filename_output, sampling_time_frames in zip(filenames_lammpstrj, filenames_edited_data, sampling_time_framess):
		
		times_frames = utils.get_num_frames( dir_lammpstrj, filename_input )
		print('\nfilename_input      : ', filename_input)
		print('Time frames           : ', times_frames)
		print('Sampling time frames  : ', sampling_time_frames)
		
		
	for filename_input, filename_output, sampling_time_frames in zip(filenames_lammpstrj, filenames_edited_data, sampling_time_framess):
		
		print('\nfilename_input      : ', filename_input)
		
		# Centering
		if calc_center == 'sampling time average':
			center = utils.get_average_center_of_mass( dir_lammpstrj, filename_input, sampling_time_frames )
		elif calc_center == 'none':
			center = utils.center_np
		print('center : ', center)
		
		
		for target_frame in sampling_time_frames:
			
			# Load data
			types, positions_grid_coord, ids_molecule, mc_step = utils.load_lammpstrj(dir_lammpstrj, filename_input, target_frame)
			print('target_frame : ', target_frame)
			print('mc_step      : ', mc_step)
			
			# Centering
			if calc_center == 'each step':
				center = utils.get_center_of_mass(types, positions_grid_coord)
			positions_real_coord = utils.centering(positions_grid_coord, center)
			
			# Get concs and condensates
			d = utils.get_concs_and_condensates(types, positions_real_coord, ids_molecule, sigma)
			
			# RDF
			rdf, rdf_bins, rdf_sampling_frames = \
				get_rdfs( dir_lammpstrj, filename_input, target_frame, center )
			d['rdf_bins'] = rdf_bins
			d['rdf_sampling_frames'] = rdf_sampling_frames
			d['rdf'] = rdf
			
			# Time info
			d['time_frame'] = target_frame
			d['mc_step'] = mc_step
			
			# Save the edited data
			prefix = filename_output
			suffix = 'frame_{}'.format(target_frame)
			utils.save(dir_edited_data, prefix, suffix, d)

