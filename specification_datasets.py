
import os, sys, glob, pickle, pprint
import numpy as np


def set_frames_photobleach_colony2( valency, length ):
	
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
	
	
def set_frames_photobleach_colony3( valencies,lengths ):
	
	num_frames_after_photobleach  = 100
	num_frames_before_photobleach = 30
	
	return num_frames_before_photobleach, num_frames_after_photobleach
	
	
class SpecDatasets():
		
	def __init__( self ):
		pass
		
		
	def _shared4( self, dir_target ):
		
		self.dir_lammpstrj    = os.path.join( '..', 'lammpstrj4', dir_target )
		self.dir_edited_data  = os.path.join( 'data4', dir_target )
		self.dir_imgs_root    = os.path.join( 'imgs4', dir_target )
		os.makedirs(self.dir_edited_data, exist_ok=True)
		
		
	def _shared5( self, dir_target ):
		
		self.dir_lammpstrj    = os.path.join( '..', 'lammpstrj5', dir_target )
		self.dir_edited_data  = os.path.join( 'data5', dir_target )
		self.dir_imgs_root    = os.path.join( 'imgs5', dir_target )
		os.makedirs(self.dir_edited_data, exist_ok=True)
		
		
	def filename_lammpstrj_matrix_valency_length(self, valency, length):
		subdir    = 'val_{}'.format(valency)
		filename  = 'R2_{}.lammpstrj'.format( str(length).zfill(3) )
		return os.path.join(subdir, filename)
		
		
	def valency_length_small_colony2( self ):
		
		self.valencies = range(2,14,2)
		self.lengths   = range(7)
		subdirs    = ['val_{}'.format(i) for i in self.valencies]
		filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in self.lengths]
		self.filenames_lammpstrj = [ os.path.join(d, f) for d in subdirs for f in filenames]
		self.filenames_edited    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in self.valencies for id_f in self.lengths]
		
		dir_target       = 'small_colony2'
		self._shared4( dir_target )
		
		self.filename_edited_matrix    = lambda valency, length: str(valency).zfill(2)+'_'+str(length).zfill(3)
		self.filename_lammpstrj_matrix = self.filename_lammpstrj_matrix_valency_length
		
		self.valencies_frap = range(4,14,2)
		self.lengths_frap   = range(7) # [1, 2, 3, 4, 5, 6, 9]
		self.set_frames_before_after_photobleach = set_frames_photobleach_colony2
		
		
		
	def valency_length_small_colony3( self ):
		
		self.valencies = range(2,14,2)
		self.lengths   = range(7)
		subdirs    = ['val_{}'.format(i) for i in self.valencies]
		filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in self.lengths]
		self.filenames_lammpstrj = [ os.path.join(d, f) for d in subdirs for f in filenames]
		self.filenames_edited    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in self.valencies for id_f in self.lengths]
		
		dir_target       = 'small_colony3'
		self._shared4( dir_target )
		
		self.filename_edited_matrix    = lambda valency, length: str(valency).zfill(2)+'_'+str(length).zfill(3)
		self.filename_lammpstrj_matrix = self.filename_lammpstrj_matrix_valency_length
		
		self.valencies_frap = range(4,14,2)
		self.lengths_frap   = range(7) # [1, 2, 3, 4, 5, 6, 9]
		self.set_frames_before_after_photobleach = set_frames_photobleach_colony3
		
		
		
	def valency_length( self ):
		
		self.valencies = range(2,16,2)
		self.lengths   = range(7)
		subdirs    = ['val_{}'.format(i) for i in self.valencies]
		filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in self.lengths]
		self.filenames_lammpstrj  = [ os.path.join(d, f) for d in subdirs for f in filenames]
		self.filenames_edited = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in self.valencies for id_f in self.lengths ]
		dir_target       = 'valency_length'
		self._shared4( dir_target )
		
		self.filename_edited_matrix    = lambda valency, length: str(valency).zfill(2)+'_'+str(length).zfill(3)
		self.filename_lammpstrj_matrix = self.filename_lammpstrj_matrix_valency_length
		
		
	def CG_valency_length( self ):
		
		self.valencies = range(2,14,2)
		self.lengths   = range(7)
		
		subdirs    = ['val_{}'.format(i) for i in self.valencies]
		filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in self.lengths]
		self.filenames_lammpstrj = [ os.path.join(d, f) for d in subdirs for f in filenames]
		self.filenames_edited    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in self.valencies for id_f in self.lengths ]
		dir_target = 'CG_valency_length'
		self._shared4( dir_target )
		
		self.filename_edited_matrix    = lambda valency, length: str(valency).zfill(2)+'_'+str(length).zfill(3)
		self.filename_lammpstrj_matrix = self.filename_lammpstrj_matrix_valency_length
		
		# Plot surface tension.
		self.real_lengths  = [1,2,3,4,5,6,9]
		self.real_linker_length_from_filename = \
			{self.filename_edited_matrix(v, l): real_length \
				for v in self.valencies \
				for l, real_length in zip(self.lengths, self.real_lengths) }
		
		# matrix pyvista
		self.sub_valencies = [2, 4, 6, 8, 12] 
		self.sub_lengths   = [1, 3, 4, 5, 6]
		
		
		# matrix pyvista
		self.surface_tension_valencies = range(4,14,2)
		self.surface_tension_lengths   = range(7)
		self.surface_tension_real_lengths  = [1,2,3,4,5,6,9]
		
		
	def CG_valency_length_only_local_move( self ):
		
		self.valencies = range(2,16,2)
		self.lengths   = range(7)
		subdirs    = ['val_{}'.format(i) for i in self.valencies]
		filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in self.lengths]
		self.filenames_lammpstrj  = [ os.path.join(d, f) for d in subdirs for f in filenames]
		self.filenames_edited = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in self.valencies for id_f in self.lengths ]
		dir_target       = 'CG_valency_length_only_local_move'
		self._shared5( dir_target )
		
		self.filename_edited_matrix    = lambda valency, length: str(valency).zfill(2)+'_'+str(length).zfill(3)
		self.filename_lammpstrj_matrix = self.filename_lammpstrj_matrix_valency_length
		
		self.valencies_frap = range(4,14,2)
		self.lengths_frap   = range(7) # [1, 2, 3, 4, 5, 6, 9]
		self.set_frames_before_after_photobleach = set_frames_photobleach_colony2
		
		
	def C_valency_length_FRAP_Control( self ):
		
		self.valencies = range(2,16,2)
		self.lengths   = range(7)
		subdirs    = ['val_{}'.format(i) for i in self.valencies]
		filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in self.lengths]
		self.filenames_lammpstrj  = [ os.path.join(d, f) for d in subdirs for f in filenames]
		self.filenames_edited = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in self.valencies for id_f in self.lengths ]
		dir_target       = 'C_valency_length_FRAP_Control'
		self._shared5( dir_target )
		
		self.filename_edited_matrix    = lambda valency, length: str(valency).zfill(2)+'_'+str(length).zfill(3)
		self.filename_lammpstrj_matrix = self.filename_lammpstrj_matrix_valency_length
		
		self.valencies_frap = range(4,14,2)
		self.lengths_frap   = range(7) # [1, 2, 3, 4, 5, 6, 9]
		self.set_frames_before_after_photobleach = set_frames_photobleach_colony2
		
		
	def conc_dependence( self ):
		
		self.filenames_edited     = [str(i).zfill(3) for i in range(81) ]
		self.filenames_lammpstrj  = ['R2_{}.lammpstrj'.format(f) for f in self.filenames_edited ]
		dir_target  = 'conc_dependence'
		self._shared4( dir_target )
		
		
	def inhibitor( self ):
		
		self.filenames_edited    = [str(i).zfill(3) for i in range(3) ]
		self.filenames_lammpstrj = ['R2_{}.lammpstrj'.format(f) for f in self.filenames_edited ]
		dir_target  = 'cam2kn1'
		self._shared4( dir_target )
		
		
	def conc_dependence_033( self ):
		
		self.filenames_lammpstrj = [str(i).zfill(3) for i in range(9) ]
		self.filenames_edited    = ['R2_{}.lammpstrj'.format(f) for f in self.filenames_edited ]
		dir_target          = 'conc_dependence_0.33'
		self._shared4( dir_target )
		
		
	def conc_dependence_merged( self ):
		
		self.stgs    = range(10)
		self.glun2bs = range(9)
		
		self.filenames_edited = [str(stg).zfill(2)+'_'+str(glun2b).zfill(2) for stg in self.stgs for glun2b in self.glun2bs ]
		dir_target  = 'conc_dependence_merged'
		self._shared4( dir_target )
		
		self.filename_edited_matrix    = lambda stg, glun2b: str(stg).zfill(2)+'_'+str(glun2b).zfill(2)
		
		# All
		#STGs    = [108,216,432,576,864,1728,2592,3456,4320,5184]
		#GluN2Bs = [270,540,1080,2160,4320,6480,8640,12960,17280]
		#
		# Subset
		# self.stgs   [216, 576, 1728, 3456] 
		# sub_glun2bs [540, 4320, 8640,17280]
		
		self.sub_stgs    = [1, 3, 5, 7]
		self.sub_glun2bs = [1, 4, 6, 8]
		
		
	def conc_dependence_merged_sub( self ):
		self.filenames_edited = [str(stg).zfill(2)+'_'+str(glun2b).zfill(2) for stg in range(8,10) for glun2b in range(4,6) ]
		dir_target  = 'conc_dependence_merged'
		self._shared4( dir_target )
		
		
	def boundary_conditions4( self ):
		
		self.filenames_edited = ['CPG', 'SP', 'SPG','CG_000','CG_001','CG_002','CG_003']
		self.filenames_lammpstrj = [os.path.join('CPG','R2_000.lammpstrj'), \
							os.path.join('SP','R2_000.lammpstrj'), \
							os.path.join('SPG','R2_000.lammpstrj'), \
							os.path.join('binary_CG','R2_000.lammpstrj'), \
							os.path.join('binary_CG','R2_001.lammpstrj'), \
							os.path.join('binary_CG','R2_002.lammpstrj'), \
							os.path.join('binary_CG','R2_003.lammpstrj') \
							]
		dir_target = 'boundary_conditions'
		self._shared4( dir_target )
		
		
	def boundary_conditions5( self ):
		
		self.filenames_edited = ['CG','CPG',\
			'PG000','PG001','PG002','PG003','PG004','PG005', \
			'SP', \
			'SPG001', 'SPG002', 'SPG003', 'SPG004', 'SPG005']
		self.filenames_lammpstrj = [os.path.join('CG','R2_000.lammpstrj'), \
							os.path.join('CPG','R2_000.lammpstrj'), \
							os.path.join('PG','R2_000.lammpstrj'), \
							os.path.join('PG','R2_001.lammpstrj'), \
							os.path.join('PG','R2_002.lammpstrj'), \
							os.path.join('PG','R2_003.lammpstrj'), \
							os.path.join('PG','R2_004.lammpstrj'), \
							os.path.join('PG','R2_005.lammpstrj'), \
							os.path.join('SP','R2_000.lammpstrj'), \
							os.path.join('SPG','R2_001.lammpstrj'), \
							os.path.join('SPG','R2_002.lammpstrj'), \
							os.path.join('SPG','R2_003.lammpstrj'), \
							os.path.join('SPG','R2_004.lammpstrj'), \
							os.path.join('SPG','R2_005.lammpstrj') \
							]
		dir_target = 'boundary_conditions'
		self._shared5( dir_target )
		
		
	def boundary_conditions2( self ):
		
		self.filenames_edited = ['SPG0','SPG1','SPG2','SPG3']
		self.filenames_lammpstrj = [os.path.join('SPG', 'R{}_trj.lammpstrj'.format(i)) for i in [0,1,2,3]]
		dir_target = 'boundary_conditions'
		self._shared4( dir_target )
		
		
	def boundary_conditions3( self ):
		
		self.filenames_edited = [ 'PG_{}'.format( str(i).zfill(3) ) for i in range(6) ]
		self.filenames_lammpstrj = [ os.path.join( 'PG','R2_{}.lammpstrj'.format(str(i).zfill(3)) ) for i in range(6) ]
		dir_target = 'boundary_conditions'
		self._shared4( dir_target )
		
		
		
		
