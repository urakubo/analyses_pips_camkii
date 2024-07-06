
import os, sys, glob, pickle, pprint
import numpy as np

import skimage
from scipy import ndimage

import networkx as nx

os.environ['OVITO_GUI_MODE'] = '1'
from ovito.io import import_file

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph

from specification_datasets import SpecDatasets
	
	
	
def process_and_save2(types, positions_grid_coord, ids_molecule, bp, mc_step, sampling_frame, \
					dir_edited_data, filename_edited, suffix ):

	# Centering
	# reference_molecule_for_centering = 'CaMKII'
	center_of_mass = utils.get_center_of_mass(types, positions_grid_coord)
	positions_real_coord = utils.centering(positions_grid_coord, center_of_mass)
	
	# Generate graph
	multi_graph = utils_graph.get_multi_graph(ids_molecule, types, bp, positions_real_coord)
	d = {}
	d['mc_step']        = mc_step
	d['sampling_frame'] = sampling_frame
	d['multi_graph']    = multi_graph
	#d['simple_graph_CaMKII_GluN2B'] = simple_graph_CaMKII_GluN2B
	
	# Save the edited data
	prefix = filename_edited
	utils.save(dir_edited_data, prefix, suffix, d)
	
	
	
class MakeConnectivityGraphTimeDev(SpecDatasets):
	def __init__( self ):
		
		self.max_num_backward_frames_for_sampling = 80
		self.num_skip_frames_for_sampling = 10
		
		
	def inspect( self ):
		
		for i, (f_lammpstrj, f_edited) in enumerate( zip(self.filenames_lammpstrj, self.filenames_edited) ):
			print('ID {}: {}, {}'.format(i, f_lammpstrj, f_edited))
		print()
		
		
	def repeat_for_valency_length( self ):
		
		suffix = 'connectivity_graph'
		for filename_lammpstrj, filename_edited in zip(self.filenames_lammpstrj, self.filenames_edited):
			
			# Load data
			sampling_frame = utils.get_num_frames(self.dir_lammpstrj, filename_lammpstrj)
			
			types, positions,ids_molecule, mc_step = \
				utils.load_lammpstrj( self.dir_lammpstrj, filename_lammpstrj, sampling_frame )
			bp = utils.load_lammpstrj_binding_partners( self.dir_lammpstrj, filename_lammpstrj, sampling_frame )
			
			print("\n"+filename_lammpstrj)
			# print("The last timeframe was sampled: ", sampling_frame )
			
			process_and_save2(types, positions, ids_molecule, bp, mc_step, sampling_frame, \
				self.dir_edited_data, filename_edited, suffix )
		
		
	def repeat_for_time_development( self, i ):
		
		filename_lammpstrj = self.filenames_lammpstrj[i]
		filename_edited    = self.filenames_edited[i]
		print('filename_lammpstrj ', filename_lammpstrj)
		print('filename_edited    ', filename_edited   )
		
		# Load data
		data_all   = import_file(os.path.join(self.dir_lammpstrj, filename_lammpstrj), input_format= "lammps/dump" )
		num_frames = data_all.source.num_frames
		
		for i in range(0, self.max_num_backward_frames_for_sampling):
			
			# Get sampling frame
			sampling_frame = num_frames + (i - self.max_num_backward_frames_for_sampling)*self.num_skip_frames_for_sampling
			suffix = 'connectivity_graph_{}'.format(i)
			
			# Decode timeframe data
			data_target_frame   = data_all.compute(sampling_frame)
			types, positions, ids_molecule = utils.decode_data(data_target_frame)
			mc_step 			= data_target_frame.attributes['Timestep']
			bp 	= np.array( data_target_frame.particles['bp'] ).astype('int')
			
			#
			process_and_save2(types, positions, ids_molecule, bp, mc_step, sampling_frame, \
				self.dir_edited_data, filename_edited, suffix )
			
			print("\n"+filename_lammpstrj)
			print("{} / {} ".format(i, self.max_num_backward_frames_for_sampling) )
			print("The sampling timeframe was: ", sampling_frame )
			print("The sampling MC step   was: ", mc_step )
			
			
if __name__ == '__main__':
	
	obj = MakeConnectivityGraphTimeDev()
	obj.valency_length_small_colony2()
	obj.inspect()
	
	# repeat_for_valency_length(filenames_lammpstrj, filenames_edited)
	
	i = 7*5+1 # val_12\R2_002
	i = 7*4+0 # val_10\R2_000
	i = 7*5+0 # val_12\R2_000
	i = 7*4+1 # val_10\R2_001
	#i = 7*3+0 # val_08\R2_000
	obj.repeat_for_time_development( i )
	#obj.repeat_for_valency_length()
	
	
