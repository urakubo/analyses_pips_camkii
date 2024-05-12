
import os, sys, glob, pickle, pprint
import numpy as np

import lib.utils as utils
import lib.parameters as p

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
	
	has_multi_graph = False
	# Valency length
	valencies = range(2,16,2)
	lengths   = range(7)
	subdirs    = ['val_{}'.format(i) for i in valencies]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in lengths]
	filenames_input  = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_output = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in valencies for id_f in lengths ]
	dir_lammpstrj    = 'valency_length'
	dir_edited_data  = 'valency_length'
	has_energy = False
	has_multi_graph = False
	
	#filenames_input  = [filenames_input[-1]]
	#filenames_output = [filenames_output[-1]]
	
	# Shared part of initialization
	dir_lammpstrj    = os.path.join('..', 'lammpstrj4', dir_lammpstrj)
	dir_edited_data  = os.path.join('data4',dir_edited_data)
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
		
		if has_multi_graph == True:
			d_graph = utils.load(dir_edited_data, filename_output,  'connectivity_graph')
			multi_graph = d_graph['multi_graph']
		else:
			multi_graph = None
		
		
		# Centering
		center    = utils.get_center_of_mass(types, positions_grid_coord)
		positions_real_coord = utils.centering(positions_grid_coord, center)
		
		# Get concs and condensates
		d = utils.get_concs_and_condensates(types, positions_real_coord, ids_molecule, \
			energy=energy,\
			multi_graph = multi_graph,\
			sigma = 2)
		
		
		
		# RDF
		#'''
		rdf, rdf_bins, rdf_sampling_frames = \
			utils.get_rdfs( dir_lammpstrj, filename_input, sampling_frame )
		d['rdf_bins'] = rdf_bins
		d['rdf_sampling_frames'] = rdf_sampling_frames
		d['rdf'] = rdf
		#'''
		
		# RDF, CaMKII binding beads
		'''
		rdf, rdf_bins = \
			utils.get_rdf_CaMKII_bead( dir_lammpstrj, filename_input, sampling_frame )
		d['rdf_bins'] = rdf_bins
		d['rdf'] = rdf
		'''
		
		# RDF, PSD95
		'''
		rdf, rdf_bins, rdf_sampling_frames = \
			utils.get_rdfs( dir_lammpstrj, filename_input, sampling_frame, multi_graph = multi_graph )
		d['rdf_PSD95_bins'] = rdf_bins
		d['rdf_PSD95_sampling_frames'] = rdf_sampling_frames
		d['rdf_PSD95'] = rdf
		'''
		
		# Time info
		d['time_frame'] = sampling_frame
		d['mc_step']    = mc_step
		
		
		# Save the edited data
		prefix = filename_output
		suffix = 'sigma_2'
		utils.save(dir_edited_data, prefix, suffix, d)
		
