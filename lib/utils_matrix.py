
import os, sys, glob, pickle, pprint, copy, pickle
import numpy as np
import matplotlib.pyplot as plt


import lib.utils as utils
import lib.parameters as p


plt.rcParams.update(p.rc_param)


class MatrixValencyLength():
	def __init__( self ):
		plt.rcParams.update( {'font.size': 6} )
		
		
	def run( self, frap = False ):
		
		if frap == True:
			self.valencies = self.valencies_frap
			self.lengths   = self.lengths_frap
		
		self.num_rows			= len( self.valencies )
		self.num_columns			= len( self.lengths  )
		
		vals = {}
		# vals = np.zeros([self.num_rows, self.num_columns])
		self.fig  = plt.figure(figsize=(8, 8), tight_layout=True)
		#fig.subplots_adjust(wspace=0.4,  hspace=0.6)
		
		for i, v in enumerate( self.valencies ):
			for j, l in enumerate( self.lengths ):
				# Load data
				filename_edited_prefix = self.filename_edited_matrix(v, l)
				row    = self.num_rows-i-1
				column = j+1
				print('Target file: ', filename_edited_prefix, ', column: ', column, ', row: ', row)
				d      = utils.load(self.dir_edited_data, \
							self.filename_edited_matrix(v, l), \
							self.suffix)
				
				legend = (row == 0) and (column == 1)
				title = filename_edited_prefix # prefix, None
				vv, _ = self.plot_a_graph(row, column, d, title, legend = legend)
				vals[filename_edited_prefix] = vv
		return vals
		
	def save( self ):
		# Shared init
		dir_imgs = os.path.join(self.dir_imgs_root, 'matrix')
		os.makedirs(dir_imgs, exist_ok=True)
		self.fig.savefig( os.path.join(dir_imgs, self.basename + '.svg' ) )
		self.fig.savefig( os.path.join(dir_imgs, self.basename + '.png' ), dpi=150 )
		plt.show()
		plt.clf()
		plt.close(fig=self.fig)



class MatrixConcDependence():
	def __init__( self ):
		plt.rcParams.update( {'font.size': 6} )
		
		
	def run( self ):
		
		self.num_rows		= len( self.glun2bs )
		self.num_columns	= len( self.stgs  )
		
		vals = np.zeros([self.num_rows, self.num_columns])
		self.fig  = plt.figure(figsize=(10, 10), tight_layout=True)
		#fig.subplots_adjust(wspace=0.4,  hspace=0.6)
		for i, stg in enumerate(self.stgs):
			for j, glun2b in enumerate(self.glun2bs):
				# Load data
				filename_edited_prefix = self.filename_edited_matrix(stg, glun2b)
				
				title = filename_edited_prefix
				row    = self.num_rows-j-1
				column = i+1
				print('Target file: ', title, ', column: ', column, ', row: ', row)
				d = utils.load(self.dir_edited_data, filename_edited_prefix, self.suffix)
				title = filename_edited_prefix # prefix, None
				vv, _ = self.plot_a_graph(row, column, d, title)
				vals[row, column-1] = vv
		return vals
		
		
	def save( self ):
		dir_imgs = os.path.join(self.dir_imgs_root, 'matrix')
		os.makedirs(dir_imgs, exist_ok=True)
		self.fig.savefig( os.path.join( dir_imgs, self.basename + '.svg' ) )
		self.fig.savefig( os.path.join( dir_imgs, self.basename + '.png' ), dpi=150 )
		plt.show()
		plt.clf()
		plt.close(fig=self.fig)
