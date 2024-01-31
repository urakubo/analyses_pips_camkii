
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
def get_concs_core_shell_region(types, positions, ids_molecule, sigma=2):
	
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
	
	
	# The below was prepared for visualization; however, it did not contribute to the improvement.
	regions_condensate0_25 = {t: utils.get_high(concs_in_grid_mesh[t]-concs_periphery[t], th=0.25) for t in targs_molecule}
	
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
			'region_condensate0_25': regions_condensate0_25
		}
	return d
	
	
	
def watershed_segmentation( d ):
	markers = np.zeros(utils.space, dtype = int)
	loc = int(utils.space[0]/2), int(utils.space[1]/2), int(utils.space[2]/2)
	markers[ loc ]    = 1
	markers[10,10,10] = 2
	
	
	tot_volume = utils.space[0] * utils.space[1] * utils.space[2]
	targets_watershed = ['CaMKII', 'STG']
	
	labels_watershed_in_grid_mesh = {}
	ratio_volumes_watershed = {}
	for j, t_watershed in enumerate( targets_watershed ):
		concs_in_grid_mesh = d['concs_in_grid_mesh'][t_watershed]
		label_watershed = watershed(concs_in_grid_mesh, markers=markers) == 1
		labels_watershed_in_grid_mesh[t_watershed] = label_watershed
		ratio_volumes_watershed[t_watershed] = np.sum( label_watershed ) / tot_volume
	return labels_watershed_in_grid_mesh, ratio_volumes_watershed
	

def get_rdf( num_frames, dir_input, filename_input, reference_molecule_for_centering ):
	
	# Parameters
	rdf_bins        = np.arange(0, 35) # You can set "np.arange(0, 35, 2)"
	rdf_grid_points = utils.get_lattice_grids()
	num_target_frames = 10
	sampling_interval = 2
	
	# Target frames
	rdf_target_frames = list(range( num_frames - num_target_frames*sampling_interval, \
					num_frames,\
					sampling_interval))
	rdf = utils.get_rdf_from_multiple_frame(dir_input, filename_input, \
			rdf_target_frames, rdf_bins, rdf_grid_points, reference_molecule_for_centering)
	
	return rdf, rdf_bins, rdf_target_frames
	#


if __name__ == '__main__':
	# Parameters
	
	'''
	# Dataset 1
	dir_input       = './../231215TravelingSpeed/data'
	filenames_input = [		'length3.lammpstrj',\
							'length9.lammpstrj',\
							'Partial_engulfment.lammpstrj',\
							'Homogeneous.lammpstrj']
	
	dir_output       = 'data'
	filenames_output = [	'PIPS',\
							'iPIPS',\
							'PartialE',\
							'Homo']
	'''
	
	
	#'''
	# Dataset 2
	dir_input        = './../231018Graphs/data'
	filenames_input  = ['R3_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(30) ]
	dir_output       = 'data'
	filenames_output = [str(i).zfill(3) for i in range(30) ]
	#'''
	
	
	reference_molecule_for_centering = 'All'
	for filename_input, filename_output in zip(filenames_input, filenames_output):
		# Load data
		print("\n"+filename_input)
		num_frames     = utils.get_num_frames(dir_input, filename_input)
		print("num_frames ", num_frames )
		types, positions,ids_molecule = utils.load_data( dir_input, filename_input, num_frames )
		print("The last timeframe was loaded." )
		
		# Centering
		center    = utils.get_center(types, positions, reference_molecule_for_centering)
		positions = utils.centering(positions, center)
		
		# RDF
		rdf, rdf_bins, rdf_target_frames = \
			get_rdf(num_frames, dir_input, filename_input, reference_molecule_for_centering)
		
		sigmas = [2,3,4]
		for sigma in sigmas:
			print('sigma ', sigma)
			# Get concs
			d = get_concs_core_shell_region(types, positions, ids_molecule, sigma)
			# Watershed
			labels_watershed_in_grid_mesh, ratio_volumes_watershed = watershed_segmentation( d )
			d['labels_watershed_in_grid_mesh'] = labels_watershed_in_grid_mesh
			d['ratio_volumes_watershed']       = ratio_volumes_watershed
			# RDF
			d['rdf_bins'] = rdf_bins
			d['rdf_target_frames'] = rdf_target_frames
			d['rdf'] = rdf
			
			# Save the edited data
			prefix = filename_output
			suffix = 'sigma_{}'.format(sigma)
			utils.save(dir_output, prefix, suffix, d)

