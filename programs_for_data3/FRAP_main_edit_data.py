


import os, sys, glob, pickle, pprint, copy, pickle
import numpy as np

from ovito.io import import_file


import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot, patches

import networkx as nx

import utils
import parameters as p
import colormap as c

rm = 15

def get_molecular_ids_in_unbleached_area(position):
	x, y, z               = position[:,0], position[:,1], position[:,2]
	xc, yc, zc, rm_       = 0, rm, 0, rm
	ids_unbleached_moleules = ( (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc) > rm_ * rm_ )
	return ids_unbleached_moleules



if __name__ == '__main__':
	
	##
	## Define Input
	##
	
	# CG_ Valency length
	#'''
	subdirs    = ['val{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
	filenames_input  = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_output = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
	dir_input    = 'CG_valency_length'
	#'''
	
	#print('Input data')
	#pprint.pprint(filenames_input)
	
	frame_photobleach = 5
	
	# Shared init
	dir_lammpstrj    = os.path.join('..', 'lammpstrj3', dir_input)
	dir_edited_data  = os.path.join('data3',dir_input)
	os.makedirs(dir_edited_data, exist_ok=True)
	
	##
	##
	##
	##
	for filename_input, filename_output in zip(filenames_input, filenames_output):
		
		
		# Load data
		print("\n"+filename_input)
		num_frames = utils.get_num_frames(dir_lammpstrj, filename_input)
		print("The last timeframe was loaded for the centering. Num frames :", num_frames)
		
		types, positions_grid_coord,ids_molecule, mc_step = \
			utils.load_lammpstrj( dir_lammpstrj, filename_input, num_frames )
		center    = utils.get_center_of_mass(types, positions_grid_coord)
		
		print('Load data.')
		data_all   = import_file(os.path.join(dir_lammpstrj, filename_input), input_format= "lammps/dump" )
		time_steps = np.array([data_all.compute(i).attributes['Timestep'] for i in range(num_frames)])
		time_photobleach= time_steps[frame_photobleach]
		#sys.exit(0)
		
		types_     = []
		positions_ = []
		molecular_numbers_in_target_area = []
		for i in range(num_frames):
			
			# Prepare data
			data_frame        = data_all.compute(i)
			types, positions_grid_coord, ids_molecule = utils.decode_data(data_frame)
			position_centered = utils.centering(positions_grid_coord, center)
			
			# Photobleach
			if i == frame_photobleach:
				ids_unbleached_moleules = get_molecular_ids_in_unbleached_area(position_centered)
			
			# Remove the bleached molecules
			if i >= frame_photobleach:
				types             = types[ids_unbleached_moleules]
				position_centered = position_centered[ids_unbleached_moleules,:]
			
			# Store the molecular numbers in the bleached area.
			ids_molecules_bleached_area   = ~get_molecular_ids_in_unbleached_area(position_centered)
			types_molecules_bleached_area = types[ids_molecules_bleached_area]
			max_id = 7
			molecular_numbers_in_target_area.append([np.count_nonzero(types_molecules_bleached_area == i) for i in range(max_id)])	
			
			# Store the location of all unbleached molecules.
			types_.append( types )
			positions_.append( position_centered )
		
		
		print("Save concentration data.")
		molecular_concentration_in_target_area = np.array(molecular_numbers_in_target_area) / (4/3*np.pi*rm*rm*rm)
		
		d = {}
		d['time_steps'] = time_steps - time_photobleach
		d['molecular_concentration_in_target_area'] = molecular_concentration_in_target_area
		d['types'] = types_
		d['positions'] = positions_
		
		
		# Save the edited data
		prefix = filename_output
		suffix = 'FRAP'
		utils.save(dir_edited_data, prefix, suffix, d)
		
	
	
