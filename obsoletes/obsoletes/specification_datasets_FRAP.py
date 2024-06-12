
import os, sys, glob, pickle, pprint
import numpy as np

from specification_datasets import SpecDatasets


class SpecDatasetsFRAP(SpecDatasets):
	
	def set_frames_photobleach_colony2( self, valency, length ):
		
		if length == 0 and valency in [10, 12]:
			num_frames_after_photobleach  = 900
			num_frames_before_photobleach = 90
		elif length == 0 and valency == 8:
			num_frames_after_photobleach  = 900
			num_frames_before_photobleach = 90
		elif length >= 4:
			num_frames_after_photobleach  = 50
			num_frames_before_photobleach = 20
		else:
			num_frames_after_photobleach  = 300
			num_frames_before_photobleach = 100

		return num_frames_before_photobleach, num_frames_after_photobleach
		
		
	def set_frames_photobleach_colony3( self, valencies,lengths ):
		
		num_frames_after_photobleach  = 100
		num_frames_before_photobleach = 30
		
		return num_frames_before_photobleach, num_frames_after_photobleach
		
		
	def __init__( self, target_dataset ):
		
		# Small colony 2
		if target_dataset  == 'small_colony2':
			
			self.valency_length_small_colony2()
			self.valencies = range(4,14,2)
			self.lengths   = range(7) # [1, 2, 3, 4, 5, 6, 9]
			self.set_frames_before_after_photobleach = self.set_frames_photobleach_colony2
		
		elif target_dataset == 'small_colony3':
			
			self.valency_length_small_colony3()
			self.valencies = range(4,14,2)
			self.lengths   = range(7) # [1, 2, 3, 4, 5, 6, 9]
			self.set_frames_before_after_photobleach = self.set_frames_photobleach_colony3
			
		elif target_dataset == 'CG_valency_length_only_local_move':
			
			self.CG_valency_length_only_local_move()
			self.valencies = range(4,14,2)
			self.lengths   = range(7) # [1, 2, 3, 4, 5, 6, 9]
			self.set_frames_before_after_photobleach = self.set_frames_photobleach_colony3
