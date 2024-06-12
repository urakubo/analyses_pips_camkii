
import os, sys, glob, pickle, pprint
import numpy as np
import utils
import parameters as p

if __name__ == '__main__':
	
	
	# Dataset1
	'''
	filenames_output = [str(i).zfill(3) for i in range(24) ] # ['024']
	filenames_input  = ['R1_{}.lammpstrj'.format(f) for f in filenames_output ]
	target_dir = 'valency_linker_length'
	dir_lammpstrj    = os.path.join('..', 'lammpstrj', target_dir)
	dir_edited_data  = os.path.join('data',target_dir)
	'''
	
	# Dataset2
	subdirs    = ['val{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(5)]
	filenames_input  = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_output = [str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(5) ]

	target_dir = 'valency_linker_length'
	dir_lammpstrj    = os.path.join('..', 'lammpstrj2', 'Feb_Figure3')
	dir_edited_data  = os.path.join('data2',target_dir)
	
	
	for filename_input, filename_output in zip(filenames_input, filenames_output):
		# Load data
		print("\n"+filename_input)
		sampling_frame = utils.get_num_frames(dir_lammpstrj, filename_input)
		
		print("sampling_frame ", sampling_frame )
		types, positions_grid_coord,ids_molecule, mc_step = \
			utils.load_lammpstrj( dir_lammpstrj, filename_input, sampling_frame )
		energy = \
			utils.load_lammpstrj_binding_energy( dir_lammpstrj, filename_input, sampling_frame )
		print("The last timeframe was loaded." )
		
		# Centering
		center    = utils.get_center_of_mass(types, positions_grid_coord)
		positions_real_coord = utils.centering(positions_grid_coord, center)
		
		# Get concs and condensates
		d = utils.get_concs_and_condensates(types, positions_real_coord, ids_molecule, energy=energy, sigma = 2)
		
		
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
		
