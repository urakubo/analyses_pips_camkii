
import os, glob, pickle, pprint, copy
import numpy as np
import pyvista

import lib.utils as utils
import lib.parameters as p
import lib.utils_pyvista as utils_pyvista

from specification_datasets import SpecDatasets


class PlotMatrixValencyLengthPyvista(SpecDatasets):
	def __init__( self ):
		pass
		
	def run( self, sub = False ):
		
		if sub == True:
			valencies = self.sub_valencies
			lengths   = self.sub_lengths
			filename  = 'matrix_pyvista_sub.png'
		else:
			valencies = self.valencies
			lengths   = self.lengths
			filename  = 'matrix_pyvista_all.png'
		
		num_rows		= len( valencies )
		num_columns		= len( lengths )
		
		suffix = 'sigma_2'
		
		pl = pyvista.Plotter(window_size=[1500,1500], shape=(num_rows, num_columns), border=False)
		
		
		for i, v in enumerate(valencies):
			for j, l in enumerate(lengths):
				# Load data
				filename_edited_prefix = self.filename_edited_matrix(v, l)
				d      = utils.load(self.dir_edited_data, filename_edited_prefix, suffix)
				print('Target: {}'.format(filename_edited_prefix) )
				
				# Specify row and column
				row    = num_rows - i - 1
				column = j
				
				pl.subplot(row, column)
				#utils_pyvista.plot_a_pre_rotated_condensate_pyvista(d, pl)
				utils_pyvista.plot_a_condensate_pyvista(d, pl, rotation=False)

				pl.camera_position = [(150.0, 0.0, 0.0),\
					 (7.0, -0.2, 2.5),\
					 (0.0, -1.0, 0.0)]
				
				#pl.camera.Zoom(1)
				
		pl.show(interactive=False, auto_close=True) # off_screen = True
		#pl.show(interactive=True, auto_close=False)
		dir_imgs = os.path.join(self.dir_imgs_root, 'matrix_3d_pyvista')
		os.makedirs(dir_imgs, exist_ok=True)
		filename = os.path.join(dir_imgs, filename)
		pl.screenshot(filename)


class PlotMatrixConcPyvista(SpecDatasets):
	def __init__( self ):
		pass
		
	def run(self, sub = False): 
		
		if sub == True:
			glun2bs   = self.sub_glun2bs
			stgs      = self.sub_stgs
			filename  = 'matrix_pyvista_sub.png'
		else:
			glun2bs   = self.glun2bs
			stgs      = self.stgs
			filename  = 'matrix_pyvista_all.png'
		
		num_rows		= len( glun2bs )
		num_columns		= len( stgs    )
		
		suffix = 'sigma_2'
		
		
		pl = pyvista.Plotter(window_size=[1500,1500], \
			shape=(num_rows, num_columns),\
			border=False)
		
		for i, stg in enumerate(stgs):
			for j, glun2b in enumerate(glun2bs):
				# Load data
				filename_edited_prefix = self.filename_edited_matrix(stg, glun2b)
				d      = utils.load(self.dir_edited_data, filename_edited_prefix, suffix)
				print('Target: ', filename_edited_prefix)
				row    = i
				column = num_rows-j-1
				pl.subplot(column, row)
				#utils.plot_a_condensate_pyvista(d, pl)
				utils_pyvista.plot_a_pre_rotated_condensate_pyvista(d, pl)
				
				pl.camera_position = [(150.0, 0.0, 0.0),\
					 (7.0, -0.2, 2.5),\
					 (0.0, -1.0, 0.0)]
				
				#pl.camera.Zoom(1)
				
		pl.show(interactive=False, auto_close=True) # off_screen = True
		#pl.show(interactive=True, auto_close=False) # off_screen = True
		
		dir_imgs = os.path.join(self.dir_imgs_root, 'matrix_3d_pyvista')
		os.makedirs(dir_imgs, exist_ok=True)
		filename = os.path.join(dir_imgs, filename)
		pl.screenshot(filename)

	
	
if __name__ == '__main__':
	pass
	
