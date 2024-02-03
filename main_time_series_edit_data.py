
import os, sys, glob, pickle, pprint
import numpy as np
from scipy import ndimage

from skimage.segmentation import watershed
import utils
	
	
	#
	# real_coord: coordinate of molecules in the real space 
	# (float) [(x0,y0,z0), (x1,y1,z1), ..., (xn,yn,zn)]
	#
	# grid_coord: coordinate of molecules in the grid space (0, 1,..., 119) 
	# (uint) [(x0,y0,z0), (x1,y1,z1), ..., (xn,yn,zn)]
	#
	# grid_mesh : molecules are embedded in a 3D numpy variable
	# (True/False in the 3D space (util.space[0], util.space[1], util.space[2]))
	#
def get_concs_and_condensates(types, positions, ids_molecule, sigma=2):
	
	# Parameters
	targs_molecule  = utils.molecules_with_all.keys() # ['GluN2B', 'CaMKII', 'STG', 'PSD95', 'All']
	
	# Get the locations of target molecules separately.
	flag_type          = {t: [True if k in utils.molecules_with_all[t]['id'] else False for k in types] for t in targs_molecule}
	locs_in_real_coord = {k: positions[flag_type[k],:] for k in targs_molecule } 
	locs_in_grid_mesh  = {k: utils.get_hist(locs_in_real_coord[k]) for k in targs_molecule}
	
	# Get a peripheral region (a region outside of a sphere) 
	# and obtain the concentrations in this area (diluted region concentration for a baseline).
	region_periphery = utils.get_periphery_in_grid_mesh()
	concs_periphery  = {t: np.sum(locs_in_grid_mesh[t] * region_periphery ) / np.sum( region_periphery ) for t in targs_molecule}
	
	# Get the condenate regions of targs_molecule.
	concs_in_grid_mesh = {t: ndimage.gaussian_filter(locs_in_grid_mesh[t], sigma = sigma) for t in targs_molecule}
	regions_condensate = {t: utils.get_high(concs_in_grid_mesh[t]-concs_periphery[t]) for t in targs_molecule}
	
	
	def get_concs_condensate(ref_molecule):
		return {t: np.sum(locs_in_grid_mesh[t] * regions_condensate[ref_molecule])/ np.sum( regions_condensate[ref_molecule] ) \
				for t in targs_molecule}
	concs_condensate = {t: get_concs_condensate(t) for t in targs_molecule}
	
	
	# Summary
	d = {
			'locs_in_grid_mesh':	locs_in_grid_mesh,
			'concs_in_grid_mesh':	concs_in_grid_mesh,
			'region_condensate':	regions_condensate,
			'conc_periphery'	:	concs_periphery,
			'conc_condensate'	:	concs_condensate,
		}
	return d
	
	
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
			
			# Get concs
			d = get_concs_and_condensates(types, positions, ids_molecule, sigma)
			
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

