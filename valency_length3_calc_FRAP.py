
import os, sys, glob, pickle, pprint, copy, pickle
import numpy as np

from ovito.io import import_file


import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph

rm = 15

def get_molecular_ids_in_unbleached_area(position):
	x, y, z               = position[:,0], position[:,1], position[:,2]
	xc, yc, zc, rm_       = 0, rm, 0, rm
	ids_unbleached_moleules = ( (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc) > rm_ * rm_ )
	return ids_unbleached_moleules


class SimulatePhotobleach():

	def set_frames_photobleach_colony2( self, valencies,lengths ):
		num_frames_after_photobleach  = {}
		num_frames_before_photobleach = {}
		for v in valencies:
			for l in lengths:
				filename_edited = str(valency).zfill(2)+'_'+str(length).zfill(3)
				if length == 0 and valency in [10, 12]:
					num_frames_after_photobleach[filename_edited]  = 900
					num_frames_before_photobleach[filename_edited] = 90
				elif length == 0 and valency == 8:
					num_frames_after_photobleach[filename_edited]  = 900
					num_frames_before_photobleach[filename_edited] = 90
				elif length >= 4:
					num_frames_after_photobleach[filename_edited]  = 50
					num_frames_before_photobleach[filename_edited] = 20
				else:
					num_frames_after_photobleach[filename_edited]  = 300
					num_frames_before_photobleach[filename_edited] = 100
		
		return num_frames_before_photobleach, num_frames_after_photobleach
		
		
	def set_frames_photobleach_colony3( self, valencies,lengths ):
		num_frames_after_photobleach  = {}
		num_frames_before_photobleach = {}
		for valency in valencies:
			for length in lengths:
				filename_edited = str(valency).zfill(2)+'_'+str(length).zfill(3)
				num_frames_after_photobleach[filename_edited]  = 100
				num_frames_before_photobleach[filename_edited] = 30
		
		return num_frames_before_photobleach, num_frames_after_photobleach
		
		
	def __init__( self ):
		
		# Small colony 2
		'''
		self.valencies = range(4,14,2)
		self.lengths   = range(7) # [1, 2, 3, 4, 5, 6, 9]
		dir_target  = 'small_colony2'
		self.num_frames_before_photobleach, self.num_frames_after_photobleach = \
			self.set_frames_photobleach_colony2(self.valencies, self.lengths)
		'''
		
		# Small colony 3
		self.valencies = range(4,14,2)
		self.lengths   = range(4,7) # [1, 2, 3, 4, 5, 6, 9]
		dir_target  = 'small_colony3'
		self.num_frames_before_photobleach, self.num_frames_after_photobleach = \
			self.set_frames_photobleach_colony3(self.valencies, self.lengths)
		
		
		# Shared init
		self.dir_lammpstrj    = os.path.join('..', 'lammpstrj4', dir_target)
		self.dir_edited_data  = os.path.join('data4',dir_target)
		os.makedirs(self.dir_edited_data, exist_ok=True)
		
		
	def repeat_runs( self ):
		
		print('Repeat runs.')
		for valency in self.valencies:
			for length in self.lengths:
				
				# Filename
				subdir    = 'val_{}'.format(valency)
				filename  = 'R2_{}.lammpstrj'.format( str(length).zfill(3) )
				filename_lammpstrj = os.path.join(subdir, filename)
				filename_edited    = str(valency).zfill(2)+'_'+str(length).zfill(3)
				
				# Load data
				print('Simulation file: ', filename_edited )
				data_all   = import_file(os.path.join(self.dir_lammpstrj, filename_lammpstrj), input_format= "lammps/dump" )
				num_frames = data_all.source.num_frames
				print('num_frames : ', num_frames)
				
				# Time
				nf_after_pb  = self.num_frames_after_photobleach[filename_edited]
				nf_before_pb = self.num_frames_before_photobleach[filename_edited]
				target_frames     = list(range(num_frames - nf_after_pb - nf_before_pb, num_frames))
				frame_photobleach = num_frames - nf_after_pb
				time_steps        = np.array([data_all.compute(t).attributes['Timestep'] for t in target_frames])
				time_photobleach  = data_all.compute(frame_photobleach).attributes['Timestep']
				
				# Run simulation
				molecular_numbers_in_target_area = \
					self.run_a_photobleach(target_frames, frame_photobleach, data_all)
				print('Finished.')
				
				# Save
				volume = (4/3*np.pi*rm*rm*rm)
				molecular_concentration_in_target_area = np.array(molecular_numbers_in_target_area) / volume
				d = {}
				d['time_steps'] = time_steps - time_photobleach
				d['molecular_concentration_in_target_area'] = molecular_concentration_in_target_area
				
				self.save_time_series(d, filename_edited)
				
				
	def save_time_series( self, d, filename_edited):
		# Save the edited data
		print("Save.")
		prefix = filename_edited
		suffix = 'FRAP'
		utils.save(self.dir_edited_data, prefix, suffix, d)
		

	def run_a_photobleach(self, target_frames, frame_photobleach, data_all):
		
		max_type_bead = 7
		molecular_numbers_in_target_area = []
		for target_frame in target_frames:
			
			# Prepare data
			data_frame        = data_all.compute(target_frame)
			types, positions_grid_coord, ids_molecule = utils.decode_data(data_frame)
			
			
			# Centering
			center    = utils.get_center_of_mass(types, positions_grid_coord)
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



if __name__ == '__main__':
	
	sim = SimulatePhotobleach()
	sim.repeat_runs()
	
