
import os, sys, glob, pickle, pprint
import numpy as np
import lib.parameters as p


from sampling_frame_photobleach import SamplingFramePhotobleach

class SpecDatasets():
		
	def examples( self ):
		'''
		Binary and quaternary mixtures.
		Molecular number: standard.
		Monte Carlo move: full set.
		Fig. 1
		'''
		self.filenames_lammpstrj = ['{}.lammpstrj'.format(f) for f in ['PS', 'CG', 'CGPS'] ]
		self.filenames_edited    = ['PS', 'CG', 'CGPS']
		
		self.dir_lammpstrj    = os.path.join('example_data','lammpstrj')
		self.dir_edited_data  = os.path.join('example_data','edited')
		self.dir_imgs_root    = os.path.join('example_data','imgs')
		
		
	def __init__( self ):
		
		pass
		
		
	def _set_directories( self, i, dir_target ):
		self.dir_lammpstrj    = os.path.join( '..', 'lammpstrj{}'.format(i), dir_target )
		self.dir_edited_data  = os.path.join( '..', 'analyses_data', 'data{}'.format(i), dir_target )
		self.dir_imgs_root    = os.path.join( '..', 'analyses_data', 'imgs{}'.format(i), dir_target )
		os.makedirs(self.dir_edited_data, exist_ok=True)
		
		
	def _filename_lammpstrj_matrix_valency_length2(self, valency, length):
		subdir    = 'val_{}'.format(valency)
		filename  = 'R2_{}.lammpstrj'.format( str(length).zfill(3) )
		return os.path.join(subdir, filename)
		
		
	def _filename_lammpstrj_matrix_valency_length3(self, valency, length):
		subdir    = 'val_{}'.format(valency)
		filename  = 'R3_{}.lammpstrj'.format( str(length).zfill(3) )
		return os.path.join(subdir, filename)
		

		
		
	def boundary_conditions( self ):
		'''
		Binary, ternary, and quaternary mixtures.
		Molecular number: standard.
		Monte Carlo move: full set.
		Fig. 1, Supplementary fig 1
		'''
		
		dirs_fnms_lammpstrj = [
			('CPG','000'),('SP' ,'000'),('SPG','000'),
			('binary_CG','000'),('binary_CG','001'),('binary_CG','002'),('binary_CG','003')]
		self.filenames_lammpstrj = [os.path.join(d,'R2_{}.lammpstrj'.format(f)) for d, f in dirs_fnms_lammpstrj ]
		self.filenames_edited    = ['CPG', 'SP', 'SPG', 'CG_000','CG_001','CG_002','CG_003']
		
		dir_target = 'boundary_conditions'
		self._set_directories( 4, dir_target )
		
		
	def boundary_conditions2( self ):
		'''
		Binary mixtures.
		Molecular number: standard.
		Monte Carlo move: full set.
		Supplementary fig 1
		'''
		self.filenames_edited = [ 'PG_{}'.format( str(i).zfill(3) ) for i in range(6) ]
		self.filenames_lammpstrj = [ os.path.join( 'PG','R2_{}.lammpstrj'.format(str(i).zfill(3)) ) for i in range(6) ]
		dir_target = 'boundary_conditions'
		self._set_directories( 4, dir_target )
		
		
	def boundary_conditions_video( self ):
		'''
		Quaternary mixtures.
		Molecular number: standard.
		Monte Carlo move: full set.
		Supplementary video 1
		Supplementary fig 1
		'''
		
		dirs_fnms_lammpstrj = [
			('CG' ,'000'),('CPG','000'),
			('PG' ,'000'),('PG' ,'001'),('PG','002') ,('PG','003') ,('PG','004') ,('PG' ,'005'),
			('SP' ,'000'),
			('SPG','001'),('SPG','002'),('SPG','003'),('SPG','004'),('SPG','005')
			]
		self.filenames_lammpstrj = [os.path.join(d,'R2_{}.lammpstrj'.format(f)) for d, f in dirs_fnms_lammpstrj]
		self.filenames_edited    = [
			'CG','CPG',\
			'PG000','PG001','PG002','PG003','PG004','PG005', \
			'SP', \
			'SPG001', 'SPG002', 'SPG003', 'SPG004', 'SPG005']
		
		dir_target = 'boundary_conditions'
		self._set_directories( 5, dir_target )
		
		
		
	def conc_dependence( self ):
		'''
		Quaternary mixture of CaMKII, GluN2B, STG, PSD95.
		Molecular number: standard.
		Monte Carlo move: full set.
		Fig 2, Supplmentary figs 2, 3
		'''
		
		self.filenames_edited     = [str(i).zfill(3) for i in range(81) ]
		self.filenames_lammpstrj  = ['R2_{}.lammpstrj'.format(f) for f in self.filenames_edited ]
		dir_target  = 'conc_dependence'
		self._set_directories( 4, dir_target )
		
		
		
	def conc_dependence_033( self ):
		'''
		Quaternary mixture of CaMKII, GluN2B, STG, PSD95.
		Molecular number: standard.
		Monte Carlo move: full set.
		Fig 2, Supplmentary figs 2, 3
		Additional set of concentration
		'''
		
		self.filenames_edited = [str(i).zfill(3) for i in range(9) ]
		self.filenames_lammpstrj  = ['R2_{}.lammpstrj'.format(f) for f in self.filenames_edited ]
		dir_target          = 'conc_dependence_0.33'
		self._set_directories( 4, dir_target )
		
		
		
	def conc_dependence_merged( self, sub = False ):
		'''
		Quaternary mixture of CaMKII, GluN2B, STG, PSD95.
		Molecular number: standard.
		Monte Carlo move: full set.
		Fig 2, Supplmentary figs 2, 3
		Merge of 'conc_dependence' and 'conc_dependence_033'
		'''
		
		self.stgs    = range(10)
		self.glun2bs = range(9)
		
		self.filenames_edited = [str(stg).zfill(2)+'_'+str(glun2b).zfill(2) for stg in self.stgs for glun2b in self.glun2bs ]
		dir_target  = 'conc_dependence_merged'
		self._set_directories( 4, dir_target )
		
		self.filename_edited_matrix = lambda glun2b, stg: str(stg).zfill(2)+'_'+str(glun2b).zfill(2)
		
		
		STGs    = [108,216,432,576,864,1728,2592,3456,4320,5184]
		GluN2Bs = [270,540,1080,2160,4320,6480,8640,12960,17280]
		volume  = np.prod(p.space_np)
		self.concs_stg    = [ s / volume * 1000 for s in STGs    ]
		self.concs_glun2b = [ n / volume * 1000 for n in GluN2Bs ]
		
		if sub == True:
			self.stgs    = [1, 3, 5, 7]
			self.glun2bs = [1, 4, 6, 8]
		self.sub = sub
		# sub_stgs    [216, 576, 1728, 3456] 
		# sub_glun2bs [540, 4320, 8640,17280]
		
		
	def boundary_conditions( self ):
		'''
		Binary, ternary, and quaternary mixtures.
		Molecular number: standard.
		Monte Carlo move: full set.
		Fig. 1, Supplementary fig 1
		'''
		
		dirs_fnms_lammpstrj = [
			('CPG','000'),('SP' ,'000'),('SPG','000'),
			('binary_CG','000'),('binary_CG','001'),('binary_CG','002'),('binary_CG','003')]
		self.filenames_lammpstrj = [os.path.join(d,'R2_{}.lammpstrj'.format(f)) for d, f in dirs_fnms_lammpstrj ]
		self.filenames_edited    = ['CPG', 'SP', 'SPG', 'CG_000','CG_001','CG_002','CG_003']
		
		dir_target = 'boundary_conditions'
		self._set_directories( 4, dir_target )
		
		
	def boundary_conditions2( self ):
		'''
		Binary mixtures.
		Molecular number: standard.
		Monte Carlo move: full set.
		Supplementary fig 1
		'''
		self.filenames_edited = [ 'PG_{}'.format( str(i).zfill(3) ) for i in range(6) ]
		self.filenames_lammpstrj = [ os.path.join( 'PG','R2_{}.lammpstrj'.format(str(i).zfill(3)) ) for i in range(6) ]
		dir_target = 'boundary_conditions'
		self._set_directories( 4, dir_target )
		
		
	def boundary_conditions5( self ):
		'''
		Quaternary mixtures.
		Molecular number: standard.
		Monte Carlo move: full set.
		Supplementary video 1
		Supplementary fig 1
		'''
		
		dirs_fnms_lammpstrj = [
			('CG' ,'000'),('CPG','000'),
			('PG' ,'000'),('PG' ,'001'),('PG','002') ,('PG','003') ,('PG','004') ,('PG' ,'005'),
			('SP' ,'000'),
			('SPG','001'),('SPG','002'),('SPG','003'),('SPG','004'),('SPG','005')
			]
		self.filenames_lammpstrj = [os.path.join(d,'R2_{}.lammpstrj'.format(f)) for d, f in dirs_fnms_lammpstrj]
		self.filenames_edited    = [
			'CG','CPG',\
			'PG000','PG001','PG002','PG003','PG004','PG005', \
			'SP', \
			'SPG001', 'SPG002', 'SPG003', 'SPG004', 'SPG005']
		
		dir_target = 'boundary_conditions'
		self._set_directories( 5, dir_target )
		
		
	def valency_length( self, sub = False ):
		'''
		Quaternary mixture of CaMKII, GluN2B, STG, PSD95.
		Molecular number: standard.
		Monte Carlo move: full set.
		Fig 4, Supplmentary fig 5
		'''
		
		self.valencies    = range(2,14,2)
		self.lengths      = range(1,7)
		self.real_lengths = [2,3,4,5,6,9]
		self.sub = sub
		if sub == True:
			self.valencies = [2, 4, 6, 8, 12] 
			self.lengths      = [1, 2, 3, 4, 5]
			self.real_lengths = [2, 3, 4, 5, 6]
		
		
		self.filename_edited_matrix    = lambda valency, length: str(valency).zfill(2)+'_'+str(length).zfill(3)
		self.filename_lammpstrj_matrix = self._filename_lammpstrj_matrix_valency_length2
		
		self.filenames_lammpstrj  = [self.filename_lammpstrj_matrix(v,l) for v in self.valencies for l in self.lengths]
		self.filenames_edited     = [self.filename_edited_matrix(v,l) for v in self.valencies for l in self.lengths]
		
		dir_target       = 'valency_length'
		self._set_directories( 4, dir_target )
		
		
		
	def CG_valency_length( self, sub = False,  sub2 = False  ):
		'''
		Binary mixture of CaMKII and GluN2B.
		Molecular number: standard.
		Monte Carlo move: full set.
		Fig 5, Supplementary Fig 6
		'''
		
		self.valencies = range(2,14,2)
		self.lengths   = range(1,7) # 
		self.real_lengths = [2, 3, 4, 5, 6, 9] # [1,2,3,4,5,6,9]
		self.sub = sub
		if   sub == True:
			self.valencies    = [2, 4, 6, 8, 12] 
			self.lengths      = [1, 2, 3, 4, 5]
			self.real_lengths = [2, 3, 4, 5, 6]
		elif sub2 == True:
			self.valencies    = [4, 6, 8, 12] 
			self.lengths      = [2, 3, 4, 5]
			self.real_lengths = [3, 4, 5, 6]
		
		self.filename_edited_matrix    = lambda valency, length: str(valency).zfill(2)+'_'+str(length).zfill(3)
		self.filename_lammpstrj_matrix = self._filename_lammpstrj_matrix_valency_length2
		
		self.filenames_lammpstrj  = [self.filename_lammpstrj_matrix(v,l) for v in self.valencies for l in self.lengths]
		self.filenames_edited     = [self.filename_edited_matrix(v,l) for v in self.valencies for l in self.lengths]
		
		dir_target = 'CG_valency_length'
		self._set_directories( 4, dir_target )
		
		# Surface tension.
		self.real_linker_length_from_filename = \
			{self.filename_edited_matrix(v, l): real_length \
				for v in self.valencies \
				for l, real_length in zip(self.lengths, self.real_lengths) }
		
		self.surface_tension_valencies     = range(4,14,2)
		self.surface_tension_lengths       = range(1,7)
		self.surface_tension_real_lengths  = [2,3,4,5,6,9]
		

		
	def CG_valency_length_only_local_move( self, frap = False ):
		'''
		Binary mixture of CaMKII and GluN2B.
		Molecular number: standard.
		Monte Carlo move: local move only.
		Fig 6, Supplementary Figs 7, 8
		'''
		
		self.valencies = range(2,14,2)
		self.lengths   = range(7)
		self.real_lengths = [1, 2, 3, 4, 5, 6, 9]
		
		if frap == True:
			'''
			self.valencies = range(4,14,2)
			self.lengths   = range(1,7) # [1, 2, 3, 4, 5, 6, 9]
			self.real_lengths = [2, 3, 4, 5, 6, 9]
			self.set_frames_before_after_photobleach = SamplingFramesPhotobleach().colony2()
			'''
			self.valencies 	  = [12]
			self.lengths      = [5]
			self.real_lengths = [6]
			self.set_frames_before_after_photobleach =  SamplingFramesPhotobleach().colony2_tmp()
			
			
		self.filename_edited_matrix    = lambda valency, length: str(valency).zfill(2)+'_'+str(length).zfill(3)
		self.filename_lammpstrj_matrix = self._filename_lammpstrj_matrix_valency_length2
		
		dir_target       = 'CG_valency_length_only_local_move'
		self._set_directories( 5, dir_target )
		
		
		self.filenames_lammpstrj  = [self.filename_lammpstrj_matrix(v,l) for v in self.valencies for l in self.lengths]
		self.filenames_edited     = [self.filename_edited_matrix(v,l) for v in self.valencies for l in self.lengths]
		self.set_frames_before_after = [self.set_frames_before_after_photobleach(v, l) for v in self.valencies for l in self.lengths]
		
		
		
	def CG_valency_length_only_local_move_fine_sampling( self, frap = False ):
		'''
		Binary mixture of CaMKII and GluN2B.
		Molecular number: standard.
		Monte Carlo move: local move only.
		Fig 6
		Sample rate: fine
		'''
		self.valencies = range(2,14,2)
		self.lengths   = range(7)
		self.real_lengths = [1, 2, 3, 4, 5, 6, 9]
		
		if frap == True:
			self.valencies = range(4,14,2)
			self.lengths   = range(1,7) # [1, 2, 3, 4, 5, 6, 9]
			self.real_lengths = [2, 3, 4, 5, 6, 9]
			
		
		self.filename_edited_matrix    = lambda valency, length: str(valency).zfill(2)+'_'+str(length).zfill(3)
		self.filename_lammpstrj_matrix = self._filename_lammpstrj_matrix_valency_length3
		
		dir_target       = 'CG_valency_length_only_local_move_fine_sampling'
		self._set_directories( 5, dir_target )
		
		self.set_frames_before_after_photobleach = SamplingFramesPhotobleach().colony5()
		
		self.filenames_lammpstrj  = [self.filename_lammpstrj_matrix(v,l) for v in self.valencies for l in self.lengths]
		self.filenames_edited     = [self.filename_edited_matrix(v,l) for v in self.valencies for l in self.lengths]
		self.set_frames_before_after = [self.set_frames_before_after_photobleach(v, l) for v in self.valencies for l in self.lengths]
		
		
		
	def C_valency_length_FRAP_Control_fine_sampling( self, frap = False ):
		'''
		CaMKII only.
		Molecular number: ?
		Monte Carlo move: local move only.
		Supplementary Fig 9
		'''
		
		self.valencies = range(4,14,2)         # range(2,14,2)
		self.lengths   = range(1,7)            # range(7)
		self.real_lengths = [2, 3, 4, 5, 6, 9] # [1, 2, 3, 4, 5, 6, 9]
		
		
		self.filename_edited_matrix    = lambda valency, length: str(valency).zfill(2)+'_'+str(length).zfill(3)
		self.filename_lammpstrj_matrix = self._filename_lammpstrj_matrix_valency_length3
		
		dir_target       = 'C_valency_length_FRAP_Control_fine_sampling'
		self._set_directories( 5, dir_target )
		
		self.set_frames_before_after_photobleach = SamplingFramesPhotobleach().colony6()
		
		self.filenames_lammpstrj  = [self.filename_lammpstrj_matrix(v,l) for v in self.valencies for l in self.lengths]
		self.filenames_edited     = [self.filename_edited_matrix(v,l) for v in self.valencies for l in self.lengths]
		self.set_frames_before_after = [self.set_frames_before_after_photobleach(v, l) for v in self.valencies for l in self.lengths]
		
		
	def conc_dependence( self ):
		'''
		Quaternary mixture of CaMKII, GluN2B, STG, PSD95.
		Molecular number: standard.
		Monte Carlo move: full set.
		Fig 2, Supplmentary figs 2, 3
		'''
		
		self.filenames_edited     = [str(i).zfill(3) for i in range(81) ]
		self.filenames_lammpstrj  = ['R2_{}.lammpstrj'.format(f) for f in self.filenames_edited ]
		dir_target  = 'conc_dependence'
		self._set_directories( 4, dir_target )
		
		
		
	def conc_dependence_033( self ):
		'''
		Quaternary mixture of CaMKII, GluN2B, STG, PSD95.
		Molecular number: standard.
		Monte Carlo move: full set.
		Fig 2, Supplmentary figs 2, 3
		Additional set of concentration
		'''
		
		self.filenames_edited = [str(i).zfill(3) for i in range(9) ]
		self.filenames_lammpstrj  = ['R2_{}.lammpstrj'.format(f) for f in self.filenames_edited ]
		dir_target          = 'conc_dependence_0.33'
		self._set_directories( 4, dir_target )
		
		
		
	def conc_dependence_merged( self, sub = False ):
		'''
		Quaternary mixture of CaMKII, GluN2B, STG, PSD95.
		Molecular number: standard.
		Monte Carlo move: full set.
		Fig 2, Supplmentary figs 2, 3
		Merge of 'conc_dependence' and 'conc_dependence_033'
		'''
		
		self.stgs    = range(10)
		self.glun2bs = range(9)
		
		self.filenames_edited = [str(stg).zfill(2)+'_'+str(glun2b).zfill(2) for stg in self.stgs for glun2b in self.glun2bs ]
		dir_target  = 'conc_dependence_merged'
		self._set_directories( 4, dir_target )
		
		self.filename_edited_matrix = lambda glun2b, stg: str(stg).zfill(2)+'_'+str(glun2b).zfill(2)
		
		
		STGs    = [108,216,432,576,864,1728,2592,3456,4320,5184]
		GluN2Bs = [270,540,1080,2160,4320,6480,8640,12960,17280]
		volume  = np.prod(p.space_np)
		self.concs_stg    = [ s / volume * 1000 for s in STGs    ]
		self.concs_glun2b = [ n / volume * 1000 for n in GluN2Bs ]
		
		if sub == True:
			self.stgs    = [1, 3, 5, 7]
			self.glun2bs = [1, 4, 6, 8]
		self.sub = sub
		# sub_stgs    [216, 576, 1728, 3456] 
		# sub_glun2bs [540, 4320, 8640,17280]
		
		
	def conc_dependence_merged_sub( self ):
		self.filenames_edited = [str(stg).zfill(2)+'_'+str(glun2b).zfill(2) for stg in range(8,10) for glun2b in range(4,6) ]
		dir_target  = 'conc_dependence_merged'
		self._set_directories( 4, dir_target )
		
