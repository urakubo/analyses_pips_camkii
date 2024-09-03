

import os, sys, glob, pickle, pprint
import numpy as np


class SamplingFramePhotobleach():
	def colony2( valency, length ):
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
	
	
	def colony2_tmp( valency, length ):
		if length == 0 and valency in [10, 12]:
			num_frames_after_photobleach  = 900
			num_frames_before_photobleach = 90
		elif length == 0 and valency == 8:
			num_frames_after_photobleach  = 900
			num_frames_before_photobleach = 90
		elif length >= 6:
			num_frames_after_photobleach  = 50
			num_frames_before_photobleach = 20
		else:
			num_frames_after_photobleach  = 500
			num_frames_before_photobleach = 100
		
		return num_frames_before_photobleach, num_frames_after_photobleach
	
	
	def colony3( valencies,lengths ):
		
		num_frames_after_photobleach  = 100
		num_frames_before_photobleach = 30
		
		return num_frames_before_photobleach, num_frames_after_photobleach
	
	
	def colony4( valency, length ):
	
		if length == 0 and valency in [10, 12]:
			num_frames_after_photobleach  = 900
			num_frames_before_photobleach = 90
		elif length == 0 and valency == 8:
			num_frames_after_photobleach  = 900
			num_frames_before_photobleach = 90
		elif length >= 6:
			num_frames_after_photobleach  = 100
			num_frames_before_photobleach = 40
		else:
			num_frames_after_photobleach  = 4000
			num_frames_before_photobleach = 400
		
		return num_frames_before_photobleach, num_frames_after_photobleach
	
	
	def colony5( valency, length ):
		
		num_frames_after_photobleach  = 9500
		num_frames_before_photobleach = 500
		
		return num_frames_before_photobleach, num_frames_after_photobleach
	
	
	def colony6( valency, length ):
		if length == 6:
			num_frames_after_photobleach  = 1000
			num_frames_before_photobleach = 100
		elif (valency == 4) and (length in [3,4,5]):
			num_frames_after_photobleach  = 1000
			num_frames_before_photobleach = 100
		elif (valency == 6) and (length in [3,4,5]):
			num_frames_after_photobleach  = 2000
			num_frames_before_photobleach = 100
		else:
			num_frames_after_photobleach  = 9500
			num_frames_before_photobleach = 500
		
		return num_frames_before_photobleach, num_frames_after_photobleach
	
	
