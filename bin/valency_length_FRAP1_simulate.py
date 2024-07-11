
import os, sys, glob, pickle, pprint, copy, pickle
import numpy as np

# from sklearn import mixture
#https://matsci.org/t/compatibility-issue-between-python-ovito-library-and-matplotlib/50794
os.environ['OVITO_GUI_MODE'] = '1'
from ovito.io import import_file

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph


#from specification_datasets_FRAP import SpecDatasetsFRAP
from specification_datasets import SpecDatasets


rm = 15

def get_molecular_ids_in_unbleached_area(position):
	x, y, z               = position[:,0], position[:,1], position[:,2]
	xc, yc, zc, rm_       = 0, rm, 0, rm
	ids_unbleached_moleules = ( (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc) > rm_ * rm_ )
	return ids_unbleached_moleules


class SimulatePhotobleach(SpecDatasets):

		
	def __init__( self ):
		#SpecDatasetsFRAP.__init__(self, target_dataset )
		self.center = None
		self.suffix = 'FRAP'
		self.num_skip_frames = 1
		
	def repeat_runs( self ):
		
		print('Repeat runs.')
		print('Target dir : ', self.dir_lammpstrj)
		for filename_lammpstrj, filename_edited, (nf_before_pb, nf_after_pb) in \
			zip( self.filenames_lammpstrj, self.filenames_edited, self.set_frames_before_after ):
			
			# Load data
			data_all   = import_file(
				os.path.join(self.dir_lammpstrj, filename_lammpstrj),
				input_format= "lammps/dump" )
			num_frames = data_all.source.num_frames
			print('num_frames     : ', num_frames)
			
			# Time
			# nf_before_pb, nf_after_pb  = self.set_frames_before_after_photobleach(valency, length)
			target_frames     = list(range(num_frames - nf_after_pb - nf_before_pb, num_frames, self.num_skip_frames))
			frame_photobleach = num_frames - nf_after_pb
			time_steps        = np.array([data_all.compute(t).attributes['Timestep'] for t in target_frames])
			time_photobleach  = data_all.compute(frame_photobleach).attributes['Timestep']
			
			# Run simulation
			print('Run simulation.')
			molecular_numbers_in_target_area = \
				self.run_a_photobleach(target_frames, frame_photobleach, data_all)
			
			# Save
			volume = (4/3*np.pi*rm*rm*rm)
			molecular_concentration_in_target_area = np.array(molecular_numbers_in_target_area) / volume
			d = {}
			d['time_steps'] = time_steps - time_photobleach
			d['molecular_concentration_in_target_area'] = molecular_concentration_in_target_area
			
			print('Save file      : ', filename_edited )
			self.save_time_series( d, filename_edited )
			print('Finished.\n')
				
				
	def save_time_series( self, d, filename_edited ):
		# Save the edited data
		prefix = filename_edited
		utils.save(self.dir_edited_data, prefix, self.suffix, d)
		
		
	def run_a_photobleach(self, target_frames, frame_photobleach, data_all):
		
		max_type_bead = 7
		molecular_numbers_in_target_area = []
		for target_frame in target_frames:
			
			# Prepare data
			data_frame        = data_all.compute(target_frame)
			types, positions_grid_coord, ids_molecule = utils.decode_data(data_frame)
			
			
			# Centering
			if self.center == None:
				center  = utils.get_center_of_mass(types, positions_grid_coord)
			else:
				center	= self.center
			position_centered = utils.centering(positions_grid_coord, center)
			
			
			# Photobleach
			if target_frame == frame_photobleach:
				ids_unbleached_moleules = get_molecular_ids_in_unbleached_area(position_centered)
			
			
			# Remove the bleached molecules
			if target_frame >= frame_photobleach:
				types             = types[ids_unbleached_moleules]
				position_centered = position_centered[ids_unbleached_moleules,:]
			
			
			# Store the molecular numbers in the bleached area.
			ids_molecules_bleached_area   = ~get_molecular_ids_in_unbleached_area(position_centered)
			types_molecules_bleached_area = types[ids_molecules_bleached_area]
			
			
			nums_molecules = [np.count_nonzero(types_molecules_bleached_area == i) for i in range(max_type_bead)]
			molecular_numbers_in_target_area.append(nums_molecules)
			
			
		return molecular_numbers_in_target_area 
	
	
