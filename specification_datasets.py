
import os, sys, glob, pickle, pprint
import numpy as np

class SpecDatasets():
		
	def __init__( self ):
		pass
		
		
	def _shared( self, dir_target ):
		
		self.dir_lammpstrj    = os.path.join( '..', 'lammpstrj4', dir_target )
		self.dir_edited_data  = os.path.join( 'data4', dir_target )
		self.dir_imgs_root    = os.path.join( 'imgs4', dir_target )
		os.makedirs(self.dir_edited_data, exist_ok=True)
		
				
	def valency_length_small_colony2( self ):
		
		valencies = range(2,14,2)
		lengths   = range(7)
		subdirs    = ['val_{}'.format(i) for i in valencies]
		filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in lengths]
		self.filenames_lammpstrj = [ os.path.join(d, f) for d in subdirs for f in filenames]
		self.filenames_edited    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in valencies for id_f in lengths]
		
		dir_target       = 'small_colony2'
		self._shared( dir_target )
		
		
	def valency_length_small_colony3( self ):
		
		valencies = range(2,14,2)
		lengths   = range(7)
		subdirs    = ['val_{}'.format(i) for i in valencies]
		filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in lengths]
		self.filenames_lammpstrj = [ os.path.join(d, f) for d in subdirs for f in filenames]
		self.filenames_edited    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in valencies for id_f in lengths]
		
		dir_target       = 'small_colony3'
		self._shared( dir_target )
		
		
	def valency_length( self ):
		
		valencies = range(2,16,2)
		lengths   = range(7)
		subdirs    = ['val_{}'.format(i) for i in valencies]
		filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in lengths]
		self.filenames_lammpstrj  = [ os.path.join(d, f) for d in subdirs for f in filenames]
		self.filenames_edited = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in valencies for id_f in lengths ]
		dir_target       = 'valency_length'
		self._shared( dir_target )
		
		
	def conc_dependence( self ):
		
		self.filenames_edited     = [str(i).zfill(3) for i in range(81) ]
		self.filenames_lammpstrj  = ['R2_{}.lammpstrj'.format(f) for f in self.filenames_edited ]
		dir_target  = 'conc_dependence'
		self._shared( dir_target )
		
		
	def inhibitor( self ):
		
		self.filenames_edited    = [str(i).zfill(3) for i in range(3) ]
		self.filenames_lammpstrj = ['R2_{}.lammpstrj'.format(f) for f in self.filenames_edited ]
		dir_target  = 'cam2kn1'
		self._shared( dir_target )
		
		
	def conc_dependence_033( self ):
		
		self.filenames_lammpstrj = [str(i).zfill(3) for i in range(9) ]
		self.filenames_edited    = ['R2_{}.lammpstrj'.format(f) for f in self.filenames_edited ]
		dir_target          = 'conc_dependence_0.33'
		self._shared( dir_target )
		
		
	def conc_dependence_merged( self ):
		
		self.filenames_edited = [str(stg).zfill(2)+'_'+str(glun2b).zfill(2) for stg in range(10) for glun2b in range(9) ]
		dir_target  = 'conc_dependence_merged'
		self._shared( dir_target )
		
		
	def boundary_conditions1( self ):
		
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
		self._shared( dir_target )
		
		
	def boundary_conditions2( self ):
		
		self.filenames_edited = ['SPG0','SPG1','SPG2','SPG3']
		self.filenames_lammpstrj = [os.path.join('SPG', 'R{}_trj.lammpstrj'.format(i)) for i in [0,1,2,3]]
		dir_target = 'boundary_conditions'
		self._shared( dir_target )
		
		
	def valency_length_CG( self ):
		subdirs    = ['val{}'.format(i) for i in range(2,14,2)]
		filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
		self.filenames_lammpstrj = [ os.path.join(d, f) for d in subdirs for f in filenames]
		self.filenames_edited    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
		dir_target = 'CG_valency_length'
		self._shared( dir_target )
		
		
		#filenames_input  = [filenames_input[7*5+2] , filenames_input[7*5+6]]
		#filenames_output = [filenames_output[7*5+2], filenames_output[7*5+6]]
	

