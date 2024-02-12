
import os, sys, glob, pickle, pprint
import numpy as np
import utils

from scipy import ndimage

	
	
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
	
	
def get_region_condensate_conc_energy(types, positions, ids_molecule, energy, sigma=2):
	# 
	energy /= 2
	
	# Parameters
	targs_molecule  = utils.molecules_with_all.keys() # ['GluN2B', 'CaMKII', 'STG', 'PSD95', 'All']
	
	# Get the locations of target molecules separately.
	flag_type          = {t: [True if k in utils.molecules_with_all[t]['id'] else False for k in types] for t in targs_molecule}
	locs_in_real_coord = {k: positions[flag_type[k],:] for k in targs_molecule }
	engy_in_real_coord = {k: energy[flag_type[k]] for k in targs_molecule }
	
	locs_in_grid_mesh  = {k: utils.get_hist(locs_in_real_coord[k]) for k in targs_molecule}
	engy_in_grid_mesh  = {k: utils.get_sum_energy(locs_in_real_coord[k], engy_in_real_coord[k]) for k in targs_molecule}
	
	# Get a peripheral region (a region outside of a sphere) 
	# and obtain the concentrations in this area (diluted region concentration for a baseline).
	region_periphery = utils.get_periphery_in_grid_mesh()
	concs_periphery  = {t: np.sum( locs_in_grid_mesh[t] * region_periphery ) / np.sum( region_periphery ) for t in targs_molecule}
	
	# Get the condenate regions of targs_molecule.
	concs_in_grid_mesh = {t: ndimage.gaussian_filter(locs_in_grid_mesh[t], sigma = sigma) for t in targs_molecule}
	regions_condensate_in_grid_mesh = {t: utils.get_high(concs_in_grid_mesh[t]-concs_periphery[t]) for t in targs_molecule}
	
	def get_concs_condensate(ref_molecule):
		return {t: np.sum(locs_in_grid_mesh[t] * regions_condensate_in_grid_mesh[ref_molecule])/ \
				np.sum( regions_condensate_in_grid_mesh[ref_molecule] ) for t in targs_molecule }
	concs_condensate = {t: get_concs_condensate(t) for t in targs_molecule}
	
	#
	def get_sum_energy_condensate(ref_molecule):
		return {t: np.sum(engy_in_grid_mesh[t] * regions_condensate_in_grid_mesh[ref_molecule]) for t in targs_molecule }
	energy_condensate = {t: get_sum_energy_condensate(t) for t in targs_molecule}
	
	# Rotate particles and obtain their concentrations and 
	rotated_in_real_coord = utils.rotate_particles_in_CaMKII_PSD95_direction( locs_in_real_coord )
	rotated_in_grid_mesh  = {k: utils.get_hist(rotated_in_real_coord[k]) for k in targs_molecule}
	rotated_concs_in_grid_mesh = \
		{t: ndimage.gaussian_filter(rotated_in_grid_mesh[t], sigma = sigma) for t in targs_molecule}
	rotated_regions_condensate_in_grid_mesh = \
		{t: utils.get_high(rotated_concs_in_grid_mesh[t]-concs_periphery[t]) for t in targs_molecule}
	
	# Summary
	d = {
			'locs_in_grid_mesh':	locs_in_grid_mesh,
			'concs_in_grid_mesh':	concs_in_grid_mesh,
			'region_condensate_in_grid_mesh': regions_condensate_in_grid_mesh,
			#
			'rotated_concs_in_grid_mesh':	rotated_concs_in_grid_mesh,
			'rotated_region_condensate_in_grid_mesh': rotated_regions_condensate_in_grid_mesh,
			#
			'conc_periphery'	:	concs_periphery,
			'conc_condensate'	:	concs_condensate,
			'energy_condensate'	:	energy_condensate
		}
	return d
	
	
if __name__ == '__main__':
	
	
	# Dataset
	dir_lammpstrj    = os.path.join('..', 'lammpstrj', 'binding_energy')
	dir_edited_data  = os.path.join('data', 'binding_energy')
	filenames_input  = ['File_PIPS_ctl.lammpstrj',\
		'File_partial_ctl.lammpstrj',\
		'File_iPIPS_length12.lammpstrj',\
		'File_homogeneous_valence4.lammpstrj',\
		'File_homogeneous_linear.lammpstrj']
	
	filenames_output = ['PIPS',\
		'partial', \
		'iPIPS', \
		'homo_valence4', \
		'homo_linear']
	
	
	
	# Init
	os.makedirs(dir_edited_data, exist_ok=True)
	
	for filename_input, filename_output in zip(filenames_input, filenames_output):
		
		# Load data
		sampling_frame = utils.get_num_frames(dir_lammpstrj, filename_input)
		
		types, positions_grid_coord,ids_molecule, mc_step, energy = \
			utils.load_lammpstrj_binding_energy( dir_lammpstrj, filename_input, sampling_frame )
		
		print("\n"+filename_input)
		print("sampling_frame ", sampling_frame )
		print("The last timeframe was loaded." )
		
		
		# Centering
		center    = utils.get_center_of_mass(types, positions_grid_coord)
		positions_real_coord = utils.centering(positions_grid_coord, center)
		
		# Get concs and energy
		d = get_region_condensate_conc_energy(types, positions_real_coord, ids_molecule, energy)
		
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
		suffix = 'sigma_{}'.format(2)
		utils.save(dir_edited_data, prefix, suffix, d)

