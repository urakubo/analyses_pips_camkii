
import os, sys, glob, pickle, pprint
import numpy as np
import math

import matplotlib
import matplotlib.pyplot as plt

sys.path.append('../')

import utils
import parameters as p
import colormap as c
	
	
def get_molecular_binding(ids_molecule, types, energy_aniso):
	
	unique_ids_molecule = np.unique( ids_molecule )
	'''
	possible_types_connection = (\
		{ p.subunits['GluN2B binding site']['id'], p.subunits['CaMKII binding site']['id'] },\
		{ p.subunits['STG binding site']['id'], p.subunits['PSD']['id'] },\
		{ p.subunits['STG binding site']['id'], p.subunits['GluN2B binding site']['id']})
	print('possible_types_connection')
	print( possible_types_connection )
	'''
	
	GluN2B_PSD95  = -4.75
	GluN2B_CaMKII = -3.5
	STG_PSD95     = -3.0
	binding_energy = {0.0: 0, -4.75: 'GluN2B_PSD95', -3.5:'GluN2B_CaMKII', -3.0:'STG_PSD95'}
	# 
	
	# Create a graph and nodes
	molecular_binding = {}
	for id in unique_ids_molecule:

		
		is_target_molecule = (id == ids_molecule)
		targ_types         = types[is_target_molecule]
		targ_energies      = energy_aniso[is_target_molecule]
		targ_positions     = positions_grid_coord[is_target_molecule,:]
		
		molecular_species = [k for k, i in p.molecules_without_all.items() if targ_types[0] in i['id']]
		molecular_species = molecular_species[0]
		
		
		targ_types_binding = [binding_energy[e] for e in targ_energies]
		nums = {b: targ_types_binding.count(b) for b in ['GluN2B_PSD95', 'GluN2B_CaMKII', 'STG_PSD95']}
		
		molecular_binding[id] = {\
			'id': id,
			'species': molecular_species, \
			'energies': targ_energies, \
			'binding_types':targ_types_binding,\
			'nums': nums
			}
		if molecular_species == 'PSD95':
			if (nums['GluN2B_PSD95'] == 0) and (nums['STG_PSD95'] == 0):
				connection = 'None'
			elif (nums['GluN2B_PSD95'] > 0) and (nums['STG_PSD95'] == 0):
				connection = 'GluN2B only'
			elif (nums['GluN2B_PSD95'] == 0) and (nums['STG_PSD95'] > 0):
				connection = 'STG only'
			elif (nums['GluN2B_PSD95'] > 0) and (nums['STG_PSD95'] > 0):
				connection = 'Both'
			molecular_binding[id]['connection'] = connection
		
	return molecular_binding
	
	
if __name__ == '__main__':
	
	#'''
	# Dataset 1: Conc dependnece
	dir_input       = 'Feb_Figure2'
	filenames_output = [str(i).zfill(3) for i in range(30) ]
	filenames_input  = ['R2_{}.lammpstrj'.format(f) for f in filenames_output ] #70
	dir_edited_data  = 'conc_dependence_cluster'
	#'''
	
	'''
	# Dataset 2:  Valency length
	subdirs    = ['val{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(5)]
	filenames_input  = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_output = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(5) ]
	dir_input        = 'Feb_Figure3'
	dir_edited_data  = 'valency_linker_length'
	'''
	
	# Shared part of initialization
	dir_lammpstrj    = os.path.join('..','..', 'lammpstrj2', dir_input)
	dir_edited_data  = os.path.join('..','data2',dir_edited_data)
	os.makedirs(dir_edited_data, exist_ok=True)
	suffix = 'connection'
	
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
		energy_aniso = energy['energy_anisotropic']
		
		# Molecular binding
		molecular_binding = get_molecular_binding(ids_molecule, types, energy_aniso)
		
		
		# Re-sort
		molecular_binding_species = {k: {} for k in p.molecules_without_all.keys()}
		for k, v in molecular_binding.items():
			molecular_binding_species[v['species']][k] = v
		
		# PSD95 shared connection
		seq  = [ v['connection'] for v in molecular_binding_species['PSD95'].values()]
		molecular_binding_species['PSD95_connection'] = \
			{k: seq.count(k) for k in ['None', 'STG only', 'GluN2B only','Both']}
		
		# Average binding number per molecule
		for species in p.molecules_without_all.keys():
			nums = {k :[ v['nums'][k] for v in molecular_binding_species[species].values() ]
				for k in ['GluN2B_PSD95', 'GluN2B_CaMKII', 'STG_PSD95']}
			molecular_binding_species[species]['nums'] = nums
		
		# Save the edited data
		prefix = filename_output
		utils.save(dir_edited_data, prefix, suffix, molecular_binding_species)
		
		
