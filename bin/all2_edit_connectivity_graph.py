
import os, sys, glob, pickle, pprint
import numpy as np


import lib.utils as utils
import lib.utils_graph as utils_graph
from specification_datasets import SpecDatasets


class EditConnectivityGraph(SpecDatasets):
	def __init__( self ):
		
		pass
		
		
	def run( self ):
		
		for filename_input, filename_output in zip(self.filenames_lammpstrj, self.filenames_edited):
			self.edit_a_dataset( filename_input, filename_output )
		
		
	def edit_a_dataset( self, filename_input, filename_output ):
		
		# Check sampling point
		print("\n"+filename_input)
		sampling_frame = utils.get_num_frames(self.dir_lammpstrj, filename_input)
		# print("The last timeframe was sampled: ", sampling_frame )
		
		
		# Load data
		types, positions_,ids_molecule, mc_step = \
			utils.load_lammpstrj( self.dir_lammpstrj, filename_input, sampling_frame )
		bp = utils.load_lammpstrj_binding_partners( self.dir_lammpstrj, filename_input, sampling_frame )
		
		
		# Centering
		center_of_mass = utils.get_center_of_mass(types, positions_) # reference_molecule_for_centering = 'CaMKII'
		positions_real_coord = utils.centering(positions_, center_of_mass)
		
		
		# Generate graph
		multi_graph = utils_graph.get_multi_graph(ids_molecule, types, bp, positions_real_coord)
		d = {}
		d['multi_graph']   = multi_graph
		d['dir_lammpstrj'] = self.dir_lammpstrj
		d['filename']      = filename_input
		d['sampling_frame']= sampling_frame
		
		for species in ['STG','GluN2B', 'PSD95','CaMKII']:
			d[species] = {}
			for type_analysis in ['average', 'distribution']:
				d[species][type_analysis] = utils_graph.get_connection_statistics(multi_graph, species, type_analysis)
		
		species = 'GluN2B'
		for type_analysis in ['CaMKII', 'PSD95']:
			d[species][type_analysis] = utils_graph.get_connection_statistics(multi_graph, species, type_analysis)
		
		species = 'PSD95'
		type_analysis = 'ratio'
		d[species][type_analysis] = utils_graph.get_connection_statistics(multi_graph, species, type_analysis)
		
		species = 'PSD95'
		type_analysis = 'ratio'
		type_analysis_= 'ratio_condensate'
		d[species][type_analysis_] = utils_graph.get_connection_statistics_condensate(multi_graph, species, type_analysis)
		
		_, ids_bead_shared_PSD = utils.get_ids_PSD95_shared_by_STG_GluN2B(multi_graph, shared_or_unshared = 'shared')
		d['ids_bead_shared_PSD'] = ids_bead_shared_PSD
		
		
		d['condensate_CaMKII'] = utils_graph.get_surface_condensate_CaMKII(types, positions_real_coord, ids_molecule)
		
		# Save the edited data
		prefix = filename_output
		suffix = 'connectivity_graph'
		utils.save(self.dir_edited_data, prefix, suffix, d)
		
	
	
if __name__ == '__main__':
	
	obj = EditConnectivityGraph()
	obj.boundary_conditions2() #  conc_dependence(), valency_length(), valency_length_CG()
	obj.run()
	
	
	
	
