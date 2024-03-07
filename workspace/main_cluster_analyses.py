
import os, sys, glob, pickle, pprint
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import networkx as nx
from networkx.utils import pairwise

sys.path.append('../')

import utils
import parameters as p
import colormap as c

	
	
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
	
def get_concs_and_condensates(types, positions, ids_molecule, energy = None, sigma=2):
	
	# Parameters
	targs_molecule  = p.molecules_with_all.keys() # ['GluN2B', 'CaMKII', 'STG', 'PSD95', 'All']
	
	# Get the locations of target molecules separately.
	flag_type          = {t: [True if k in p.molecules_with_all[t]['id'] else False for k in types] for t in targs_molecule}
	locs_in_real_coord = {k: positions[flag_type[k],:] for k in targs_molecule } 
	locs_in_grid_mesh  = {k: get_hist(locs_in_real_coord[k]) for k in targs_molecule}
	
	# Get a peripheral region (a region outside of a sphere) 
	# and obtain the concentrations in this area (diluted region concentration for a baseline).
	region_periphery = get_periphery_in_grid_mesh()
	concs_periphery  = {t: np.sum(locs_in_grid_mesh[t] * region_periphery ) / np.sum( region_periphery ) for t in targs_molecule}
	
	# Get the condenate regions of targs_molecule.
	concs_in_grid_mesh = {t: ndimage.gaussian_filter(locs_in_grid_mesh[t], sigma = sigma) for t in targs_molecule}
	regions_condensate_in_grid_mesh = {t: get_high(concs_in_grid_mesh[t]-concs_periphery[t]) for t in targs_molecule}
	regions_condensate_in_grid_mesh['dilute'] = region_periphery
	
	def get_concs_condensate(ref_molecule):
		return {t: np.sum(locs_in_grid_mesh[t] * regions_condensate_in_grid_mesh[ref_molecule])/ \
				np.sum( regions_condensate_in_grid_mesh[ref_molecule] ) for t in targs_molecule }
	concs_condensate = {t: get_concs_condensate(t) for t in targs_molecule}
	
	
	#iso        = energy['energy_isotropic']
	#aniso      = energy['energy_anisotropic']
	#aniso_self = energy['energy_anisotropic_self']
	
	for k, e in energy.items():
		e /= 2
		engy_in_real_coord = {m: e[flag_type[m]] for m in targs_molecule }
		engy_in_grid_mesh  = {m: get_sum_energy(locs_in_real_coord[m], engy_in_real_coord[m]) for m in targs_molecule}
		def get_sum_energy_condensate(ref_molecule):
			return {m: np.sum(engy_in_grid_mesh[m] * regions_condensate_in_grid_mesh[ref_molecule]) for m in targs_molecule }
		engy_condensate    = {r: get_sum_energy_condensate(r) for r in targs_molecule}
		engy_condensate['dilute'] = {m: np.sum(engy_in_grid_mesh[m] * region_periphery) for m in targs_molecule }
		
		d[k] = engy_condensate
	
	
	# Summary
	d = {
			'locs_in_grid_mesh':	locs_in_grid_mesh,
			'concs_in_grid_mesh':	concs_in_grid_mesh,
			'region_condensate_in_grid_mesh': regions_condensate_in_grid_mesh,
			#
			'conc_periphery'	:	concs_periphery,
			'conc_condensate'	:	concs_condensate,
		}

	
	return d
	
	
if __name__ == '__main__':

	i = 2
	
	# Dataset 1
	dir_input       = './../lammpstrj2/Feb_Figure2'
	fnames = [27,21,15,9]
	filenames_input  = [ 'R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in fnames ]
	dir_edited_data  = os.path.join('conc_dependence', 'cluster_analyses' )

	# Shared part of initialization
	dir_lammpstrj    = os.path.join('..','..', 'lammpstrj2', dir_lammpstrj)
	dir_edited_data  = os.path.join('data2',dir_edited_data)
	os.makedirs(dir_edited_data, exist_ok=True)
	
	#
	for filename_input, filename_output in zip(filenames_input, filenames_output):
		
		
		# Check sampling point
		print("\n"+filename_input)
		sampling_frame = utils.get_num_frames(dir_lammpstrj, filename_input)
		print("The last timeframe was sampled: ", sampling_frame )
		
		
		# Load data
		types, positions_grid_coord,ids_molecule, mc_step = \
			utils.load_lammpstrj( dir_lammpstrj, filename_input, sampling_frame )
		energy = utils.load_lammpstrj_binding_energy( dir_lammpstrj, filename_input, sampling_frame )
		
		
		# Centering
		center    = utils.get_center_of_mass(types, positions_grid_coord)
		positions_real_coord = utils.centering(positions_grid_coord, center)
		
		
		# Get concs and condensates
		d = utils.get_concs_and_condensates(types, positions_real_coord, ids_molecule, energy=energy, sigma = 2)
		
		
		# Time info
		d['time_frame'] = sampling_frame
		d['mc_step']    = mc_step
		
		
		# Save the edited data
		prefix = filename_output
		suffix = 'sigma_2'
		utils.save(dir_edited_data, prefix, suffix, d)
		
		
