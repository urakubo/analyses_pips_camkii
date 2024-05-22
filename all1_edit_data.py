
import os, sys, glob, pickle, pprint
import numpy as np

import lib.utils as utils
import lib.parameters as p
from lib.specification_datasets import SpecDatasets



class EditData(SpecDatasets):
	def __init__( self ):
		
		pass
		
		
	def run( self ):
		
		self.dir_lammpstrj    = os.path.join('..', 'lammpstrj4', self.dir_target)
		self.dir_edited_data  = os.path.join('data4',self.dir_target)
		os.makedirs(self.dir_edited_data, exist_ok=True)
		
		for filename_input, filename_output in zip(self.filenames_lammpstrj, self.filenames_edited):
			self.edit_a_dataset( filename_input, filename_output )
		
		
	def edit_a_dataset( self, filename_input, filename_output ):
		
		# Load data
		print("\n"+filename_input)
		sampling_frame = utils.get_num_frames(self.dir_lammpstrj, filename_input)
		print("The last timeframe was loaded: ", sampling_frame )
		types, positions_grid_coord,ids_molecule, mc_step = \
			utils.load_lammpstrj( self.dir_lammpstrj, filename_input, sampling_frame )
		
		
		# Centering
		center    = utils.get_center_of_mass(types, positions_grid_coord)
		positions_real_coord = utils.centering(positions_grid_coord, center)
		
		
		# Get concs and condensates
		d = utils.get_concs_and_condensates(types, positions_real_coord, ids_molecule, \
			sigma = 2)
		
		
		# RDF
		rdf, rdf_bins, rdf_sampling_frames = \
			utils.get_rdfs( self.dir_lammpstrj, filename_input, sampling_frame )
		d['rdf_bins'] = rdf_bins
		d['rdf_sampling_frames'] = rdf_sampling_frames
		d['rdf'] = rdf
		
		
		# Time info
		d['time_frame'] = sampling_frame
		d['mc_step']    = mc_step
		
		
		# Other info
		d['dir_lammpstrj'] = self.dir_lammpstrj
		d['filename_lammpstrj'] = filename_input
		d['sampling_frame']= sampling_frame
		
		
		# Save the edited data
		prefix = filename_output
		suffix = 'sigma_2'
		utils.save(self.dir_edited_data, prefix, suffix, d)
		
		
if __name__ == '__main__':
	
	obj = EditData()
	obj.boundary_conditions2() #  conc_dependence(), valency_length(), valency_length_CG()
	obj.run()
	
	
	
	

