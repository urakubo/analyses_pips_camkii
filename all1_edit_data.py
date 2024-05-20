
import os, sys, glob, pickle, pprint
import numpy as np

import lib.utils as utils
import lib.parameters as p


	
if __name__ == '__main__':
	
	
	'''
	# Valency length
	valencies = range(2,16,2)
	lengths   = range(7)
	subdirs    = ['val_{}'.format(i) for i in valencies]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in lengths]
	filenames_input  = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_output = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in valencies for id_f in lengths ]
	dir_lammpstrj    = 'valency_length'
	dir_edited_data  = 'valency_length'
	
	
	# Conc dependence
	filenames_output = [str(i).zfill(3) for i in range(81) ]
	filenames_input  = ['R2_{}.lammpstrj'.format(f) for f in filenames_output ]
	dir_lammpstrj    = 'conc_dependence'
	dir_edited_data  = 'conc_dependence'
	
	filenames_output = [str(i).zfill(3) for i in range(9) ]
	filenames_input  = ['R2_{}.lammpstrj'.format(f) for f in filenames_output ]
	dir_lammpstrj    = os.path.join('conc_dependence', '0.33')
	dir_edited_data  = 'conc_dependence_0.33'
	
	
	# other_conditions
	filenames_output = ['CPG', 'SP', 'SPG','CG_000','CG_001','CG_002','CG_003']
	filenames_input  = [os.path.join('CPG','R2_000.lammpstrj'), \
						os.path.join('SP','R2_000.lammpstrj'), \
						os.path.join('SPG','R2_000.lammpstrj'), \
						os.path.join('binary_CG','R2_000.lammpstrj'), \
						os.path.join('binary_CG','R2_001.lammpstrj'), \
						os.path.join('binary_CG','R2_002.lammpstrj'), \
						os.path.join('binary_CG','R2_003.lammpstrj') \
						]
	
	filenames_output = ['SPG','PIPS']
	filenames_input  = [os.path.join('SPG','SPG_unactivated_trj.lammpstrj'), \
						os.path.join('SPG','PIPS_activated_trj.lammpstrj')\
						]
	dir_edited_data  = 'special_conditions'
	dir_lammpstrj    = '.'
	'''
	
	
	# Valency and linker length of CG
	subdirs    = ['val{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
	filenames_input  = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_output = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
	dir_lammpstrj    = 'CG_valency_length'
	dir_edited_data  = 'CG_valency_length'
	
	#filenames_input  = [filenames_input[7*5+2] , filenames_input[7*5+6]]
	#filenames_output = [filenames_output[7*5+2], filenames_output[7*5+6]]
	
	
	# Shared part of initialization
	dir_lammpstrj    = os.path.join('..', 'lammpstrj4', dir_lammpstrj)
	dir_edited_data  = os.path.join('data4',dir_edited_data)
	os.makedirs(dir_edited_data, exist_ok=True)
	
	
	for filename_input, filename_output in zip(filenames_input, filenames_output):
		
		# Load data
		print("\n"+filename_input)
		sampling_frame = utils.get_num_frames(dir_lammpstrj, filename_input)
		print("The last timeframe was loaded: ", sampling_frame )
		types, positions_grid_coord,ids_molecule, mc_step = \
			utils.load_lammpstrj( dir_lammpstrj, filename_input, sampling_frame )
		
		
		# Centering
		center    = utils.get_center_of_mass(types, positions_grid_coord)
		positions_real_coord = utils.centering(positions_grid_coord, center)
		
		
		# Get concs and condensates
		d = utils.get_concs_and_condensates(types, positions_real_coord, ids_molecule, \
			sigma = 2)
		
		
		# RDF
		rdf, rdf_bins, rdf_sampling_frames = \
			utils.get_rdfs( dir_lammpstrj, filename_input, sampling_frame )
		d['rdf_bins'] = rdf_bins
		d['rdf_sampling_frames'] = rdf_sampling_frames
		d['rdf'] = rdf
		
		
		# Time info
		d['time_frame'] = sampling_frame
		d['mc_step']    = mc_step
		
		
		# Other info
		d['dir_lammpstrj'] = dir_lammpstrj
		d['filename_lammpstrj'] = filename_input
		d['sampling_frame']= sampling_frame
		
		
		# Save the edited data
		prefix = filename_output
		suffix = 'sigma_2'
		utils.save(dir_edited_data, prefix, suffix, d)
