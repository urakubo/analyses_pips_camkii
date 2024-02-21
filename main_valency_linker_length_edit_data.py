
import os, sys, glob, pickle, pprint
import numpy as np
import utils
import parameters as p

from scipy import ndimage
	
	
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
	
	
#def get_region_condensate_conc_energy(types, positions, ids_molecule, energy, sigma=2):
def get_region_condensate_conc_energy(types, positions, ids_molecule, sigma=2):
	# 
	#energy /= 2
	
	# Parameters
	targs_molecule  = p.molecules_with_all.keys() # ['GluN2B', 'CaMKII', 'STG', 'PSD95', 'All']
	
	# Get the locations of target molecules separately.
	flag_type          = {t: [True if k in p.molecules_with_all[t]['id'] else False for k in types] for t in targs_molecule}
	locs_in_real_coord = {k: positions[flag_type[k],:] for k in targs_molecule }
	
	locs_in_grid_mesh  = {k: utils.get_hist(locs_in_real_coord[k]) for k in targs_molecule}
	
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
	
	
	# Summary
	d = {	
			'locs_in_grid_mesh':	locs_in_grid_mesh,
			'concs_in_grid_mesh':	concs_in_grid_mesh,
			'region_condensate_in_grid_mesh': regions_condensate_in_grid_mesh,
			#
			'conc_periphery'	:	concs_periphery,
			'conc_condensate'	:	concs_condensate
		}	
	
	
	# Rotate particles and obtain the concentration.
	
	
	rotated_in_real_coord = utils.rotate_particles_in_CaMKII_PSD95_direction( locs_in_real_coord )
	rotated_in_grid_mesh  = {k: utils.get_hist(rotated_in_real_coord[k]) for k in targs_molecule}
	rotated_concs_in_grid_mesh = \
		{t: ndimage.gaussian_filter(rotated_in_grid_mesh[t], sigma = sigma) for t in targs_molecule}
	rotated_regions_condensate_in_grid_mesh = \
		{t: utils.get_high(rotated_concs_in_grid_mesh[t]-concs_periphery[t]) for t in targs_molecule}
	
	d['rotated_concs_in_grid_mesh']             = rotated_concs_in_grid_mesh
	d['rotated_region_condensate_in_grid_mesh'] = rotated_regions_condensate_in_grid_mesh
	
	
	'''
	engy_in_real_coord = {k: energy[flag_type[k]] for k in targs_molecule }
	engy_in_grid_mesh  = {k: utils.get_sum_energy(locs_in_real_coord[k], engy_in_real_coord[k]) for k in targs_molecule}
	def get_sum_energy_condensate(ref_molecule):
		return {t: np.sum(engy_in_grid_mesh[t] * regions_condensate_in_grid_mesh[ref_molecule]) for t in targs_molecule }
	energy_condensate      = {t: get_sum_energy_condensate(t) for t in targs_molecule}
	d['energy_condensate'] = energy_condensate
	'''
	
	return d
	
	
if __name__ == '__main__':
	
	
	# Dataset
	filenames_output = [str(i).zfill(3) for i in range(24) ] # ['024']
	filenames_input  = ['R1_{}.lammpstrj'.format(f) for f in filenames_output ]
	
	target_dir = 'valency_linker_length'
	
	#
	dir_lammpstrj    = os.path.join('..', 'lammpstrj', target_dir)
	dir_edited_data  = os.path.join('data',target_dir)
	
	
	for filename_input, filename_output in zip(filenames_input, filenames_output):
		
		# Load data
		sampling_frame = utils.get_num_frames(dir_lammpstrj, filename_input)
		
		#types, positions_grid_coord,ids_molecule, mc_step, energy = \
		#	utils.load_lammpstrj_binding_energy( dir_lammpstrj, filename_input, sampling_frame )
		
		types, positions_grid_coord,ids_molecule, mc_step = \
			utils.load_lammpstrj( dir_lammpstrj, filename_input, sampling_frame )
		
		print("\n"+filename_input)
		print("sampling_frame ", sampling_frame )
		print("The last timeframe was loaded." )
		
		
		# Centering
		center    = utils.get_center_of_mass(types, positions_grid_coord)
		positions_real_coord = utils.centering(positions_grid_coord, center)
		
		# Get concs and energy
		#d = get_region_condensate_conc_energy(types, positions_real_coord, ids_molecule, energy)
		d = get_region_condensate_conc_energy( types, positions_real_coord, ids_molecule )
		
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
		
