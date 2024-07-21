
import os, sys, glob, pickle, pprint, copy, pickle
import numpy as np
import matplotlib.pyplot as plt


import lib.utils as utils
import lib.parameters as p


plt.rcParams.update(p.rc_param)



class Matrix():
	def __init__( self ):
		plt.rcParams.update( {'font.size': 8} )
		#plt.rcParams.update( {'font.size': 6} )
		#pass
		
	def run( self ):
		
		
		if hasattr( self, 'valencies') and ( self, 'lengths'):
			rows	= self.valencies
			columns	= self.lengths
			figsize = (8, 8)
			figsize = (6, 6)
		elif hasattr( self, 'stgs') and ( self, 'glun2bs'):
			rows	= self.glun2bs
			columns	= self.stgs 
			figsize = (10, 10)
		else:
			ValueError("Neither conc_dependence nor valnency_length.")
		
		
		self.fig  = plt.figure(figsize=figsize, tight_layout=True)
		#fig.subplots_adjust(wspace=0.4,  hspace=0.6)
		
		self.num_rows		= len( rows )
		self.num_columns	= len( columns )
		print('num_column, num_row :', self.num_columns, self.num_rows)
		
		vals = {}
		self.vals = np.zeros([self.num_rows, self.num_columns])
		for i, column in enumerate( columns ):
			for j, row in enumerate( rows ):
				# Load data
				filename_edited_prefix = self.filename_edited_matrix(row, column)
				
				print('Target file: ', filename_edited_prefix, ', column: ', column, ', row: ', row)
				d      = utils.load(self.dir_edited_data, \
							filename_edited_prefix, \
							self.suffix)
				
				column_ = i+1
				row_    = self.num_rows-j-1
				print(' column_: ', column_, ', row_: ', row_)

				legend = (row_ == 0) and (column_ == 1)
				title = filename_edited_prefix
				vv, _ = self.plot_a_graph(row_, column_, d, title, legend = legend)
				vals[filename_edited_prefix] = vv
				self.vals[j, i] = vv
		return vals
		
		
	def save( self ):
		
		dir_imgs = os.path.join(self.dir_imgs_root, 'matrix')
		os.makedirs(dir_imgs, exist_ok=True)
		self.fig.savefig( os.path.join(dir_imgs, self.basename + '.svg' ) )
		self.fig.savefig( os.path.join(dir_imgs, self.basename + '.png' ), dpi=150 )
		plt.show()
		plt.clf()
		plt.close(fig=self.fig)

