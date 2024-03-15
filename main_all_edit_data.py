
import os, sys, glob, pickle, pprint
import numpy as np

import utils
	
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
	
	
if __name__ == '__main__':
	
	# Valency length
	'''
	subdirs    = ['val_{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
	filenames_input  = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_output = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
	dir_lammpstrj    = 'valency_length'
	dir_edited_data  = 'valency_length'
	has_energy = False
	'''
	
	# Conc dependnece
	'''
	filenames_output = [str(i).zfill(3) for i in range(48) ]
	filenames_input  = ['R2_{}.lammpstrj'.format(f) for f in filenames_output ] #70
	dir_lammpstrj    = 'conc_dependence'
	dir_edited_data  = 'conc_dependence'
	has_energy = False
	'''
	
	# Mixtures
	#'''
	filenames_output = ['CG','CPG','PG','SP','CGSP'] # ,'SPG'
	filenames_input  = ['{}.lammpstrj'.format(n) for n in filenames_output]
	dir_lammpstrj    = 'mixtures'
	dir_edited_data  = 'mixtures'
	has_energy = False
	#'''
	# Long GluN2B
	'''
	filenames_output = [str(i).zfill(3) for i in range(10,15) ]
	filenames_input  = ['R2_{}.lammpstrj'.format(f) for f in filenames_output ]
	dir_lammpstrj    = os.path.join('0')
	dir_edited_data  = 'GluN2B_length4'
	has_energy = False
	'''	
	# Short GluN2B
	'''
	filenames_output = [str(i).zfill(3) for i in range(5) ]
	filenames_input  = ['R2_{}.lammpstrj'.format(f) for f in filenames_output ]
	dir_lammpstrj    = os.path.join('0')
	dir_edited_data  = 'GluN2B_length1'
	has_energy = False
	'''
	
	# Shared part of initialization
	dir_lammpstrj    = os.path.join('..', 'lammpstrj3', dir_lammpstrj)
	dir_edited_data  = os.path.join('data3',dir_edited_data)
	os.makedirs(dir_edited_data, exist_ok=True)
	
	
	for filename_input, filename_output in zip(filenames_input, filenames_output):
		# Load data
		print("\n"+filename_input)
		sampling_frame = utils.get_num_frames(dir_lammpstrj, filename_input)
		print("The last timeframe was loaded." )
		
		print("sampling_frame ", sampling_frame )
		types, positions_grid_coord,ids_molecule, mc_step = \
			utils.load_lammpstrj( dir_lammpstrj, filename_input, sampling_frame )
		
		if has_energy == True:
			energy = \
				utils.load_lammpstrj_binding_energy( dir_lammpstrj, filename_input, sampling_frame )
		else:
			energy = None
		
		# Centering
		center    = utils.get_center_of_mass(types, positions_grid_coord)
		positions_real_coord = utils.centering(positions_grid_coord, center)
		
		# Get concs and condensates
		d = utils.get_concs_and_condensates(types, positions_real_coord, ids_molecule, energy=energy, sigma = 2)
		
		# Watershed
		'''
		labels_watershed_in_grid_mesh, ratio_volumes_watershed = utils.watershed_segmentation( d )
		d['labels_watershed_in_grid_mesh'] = labels_watershed_in_grid_mesh
		d['ratio_volumes_watershed']       = ratio_volumes_watershed
		'''
		
		
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
		suffix = 'sigma_2'
		utils.save(dir_edited_data, prefix, suffix, d)
		
