
import os, sys, glob, pickle, pprint
import numpy as np
import math

import matplotlib
import matplotlib.pyplot as plt

sys.path.append('../')

import utils
import parameters as p
import colormap as c
	
	
def get_molecular_binding(ids_molecule, types, bp):
	
	unique_ids_molecule = np.unique( ids_molecule )
	
	'''
	GluN2B_PSD95  = -4.75
	GluN2B_CaMKII = -3.5
	STG_PSD95     = -3.0
	binding_energy = {0.0: 0, -4.75: 'GluN2B_PSD95', -3.5:'GluN2B_CaMKII', -3.0:'STG_PSD95'}
	'''
	
	# Create a graph and nodes
	types_molecule_for_bead = np.zeros_like(bp, dtype = 'int') # CaMKII: 1, GluN2B: 2,  PSD95: 3, STG: 4
	id_molecular_type = {'CaMKII': 1, 'GluN2B': 2, 'PSD95': 3, 'STG': 4}
	
	for id in unique_ids_molecule:
		flags_target_molecule_for_bead = (id == ids_molecule)
		targ_types = types[flags_target_molecule_for_bead]
		molecular_species = [k for k, i in p.molecules_without_all.items() if targ_types[0] in i['id']][0]
		flag_target_molecule = (id == ids_molecule)
		types_molecule_for_bead[flag_target_molecule] = id_molecular_type[molecular_species]
	
	
	molecular_binding = {}
	for id in unique_ids_molecule:
		
		flags_target_molecule_for_bead = (id == ids_molecule)
		targ_types             = types[flags_target_molecule_for_bead]
		targ_bp			       = bp[flags_target_molecule_for_bead]
		ids_partner_bead       = targ_bp[targ_bp >=0]
		types_partner_molecule = np.array([types_molecule_for_bead[id] for id in ids_partner_bead])
		
		# print('types_partner_molecule ', types_partner_molecule)
		
		molecular_species   = [k for k, i in p.molecules_without_all.items() if targ_types[0] in i['id']][0]
		
		
		
		if molecular_species == 'CaMKII':
			nums = {	'GluN2B_CaMKII': ids_partner_bead.shape[0],\
						'GluN2B_PSD95' : 0,\
						'STG_PSD95'    : 0 \
						}
		elif molecular_species == 'GluN2B':
			nums = {	'GluN2B_CaMKII': np.sum(types_partner_molecule == id_molecular_type['CaMKII']),\
						'GluN2B_PSD95' : np.sum(types_partner_molecule == id_molecular_type['PSD95' ]),\
						'STG_PSD95'    : 0\
						}
		elif molecular_species == 'PSD95':
			nums = {	'GluN2B_CaMKII': 0,\
						'GluN2B_PSD95' : np.sum(types_partner_molecule == id_molecular_type['GluN2B']), \
						'STG_PSD95'    : np.sum(types_partner_molecule == id_molecular_type['STG'])\
						}
		elif molecular_species == 'STG':
			nums = {	'GluN2B_CaMKII': 0,\
						'GluN2B_PSD95' : 0,\
						'STG_PSD95'    : ids_partner_bead.shape[0]\
						}
		else:
			print('Error: Molecular type selection : ', molecular_species )
			
		
		molecular_binding[id] = {\
			'id': id,
			'species': molecular_species, \
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
			molecular_binding[id]['PSD95_connection'] = connection
		
	return molecular_binding
	
	
if __name__ == '__main__':
	
	'''
	# Dataset 1: Conc dependnece
	dir_input    = 'conc_dependence'
	filenames_output = [str(i).zfill(3) for i in range(48) ]
	filenames_input  = ['R2_{}.lammpstrj'.format(f) for f in filenames_output ]
	dir_edited_data  = 'conc_dependence'
	'''
	
	#'''
	# Dataset 2:  Valency length
	
	subdirs    = ['val_{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
	filenames_input  = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_output = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
	dir_edited_data  = 'valency_length'
	dir_input        = 'valency_length'
	#'''
	
	# Shared part of initialization
	dir_lammpstrj    = os.path.join('..', 'lammpstrj3', dir_input)
	dir_edited_data  = os.path.join('data3',dir_edited_data)
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
		bp = utils.load_lammpstrj_binding_partners( dir_lammpstrj, filename_input, sampling_frame )
		
		
		# Molecular binding
		molecular_binding = get_molecular_binding(ids_molecule, types, bp)
		
		
		# Re-sort
		molecular_binding_species = {k: {} for k in p.molecules_without_all.keys()}
		for k, v in molecular_binding.items():
			molecular_binding_species[v['species']][k] = v
		
		# PSD95 shared connection
		seq  = [ v['PSD95_connection'] for v in molecular_binding_species['PSD95'].values()]
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
		
		
