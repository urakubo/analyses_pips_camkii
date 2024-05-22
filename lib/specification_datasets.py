
import os, sys, glob, pickle, pprint
import numpy as np

class SpecDatasets():
		
	def valency_length( self ):
		
		valencies = range(2,16,2)
		lengths   = range(7)
		subdirs    = ['val_{}'.format(i) for i in valencies]
		filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in lengths]
		self.filenames_lammpstrj  = [ os.path.join(d, f) for d in subdirs for f in filenames]
		self.filenames_edited = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in valencies for id_f in lengths ]
		self.dir_target       = 'valency_length'
		
		
	def conc_dependence( self ):
		
		filenames_output = [str(i).zfill(3) for i in range(81) ]
		filenames_input  = ['R2_{}.lammpstrj'.format(f) for f in filenames_output ]
		dir_lammpstrj    = 'conc_dependence'
		dir_edited_data  = 'conc_dependence'
		
		
	def conc_dependence_033( self ):
		
		self.filenames_lammpstrj = [str(i).zfill(3) for i in range(9) ]
		self.filenames_edited    = ['R2_{}.lammpstrj'.format(f) for f in filenames_output ]
		self.dir_target          = 'conc_dependence_0.33'
		
		
	def conc_dependence_merged( self ):
		
		self.filenames_edited = [str(stg).zfill(2)+'_'+str(glun2b).zfill(2) for stg in range(10) for glun2b in range(9) ]
		self.dir_target  = 'conc_dependence_merged'
		
		
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
		self.dir_target = 'boundary_conditions'
		
		
	def boundary_conditions2( self ):
		
		self.filenames_edited = ['SPG0','SPG1','SPG2','SPG3']
		self.filenames_lammpstrj = [os.path.join('SPG', 'R{}_trj.lammpstrj'.format(i)) for i in [0,1,2,3]]
		self.dir_target = 'boundary_conditions'
		
		
	def valency_length_CG( self ):
		subdirs    = ['val{}'.format(i) for i in range(2,14,2)]
		filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
		self.filenames_lammpstrj = [ os.path.join(d, f) for d in subdirs for f in filenames]
		self.filenames_edited    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
		self.dir_target = 'CG_valency_length'
		
		
		#filenames_input  = [filenames_input[7*5+2] , filenames_input[7*5+6]]
		#filenames_output = [filenames_output[7*5+2], filenames_output[7*5+6]]
	

