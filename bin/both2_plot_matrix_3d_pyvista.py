
import os, sys, glob, pickle, pprint, copy
import numpy as np
import pyvista

import lib.utils as utils
import lib.parameters as p
import lib.utils_pyvista as utils_pyvista

from specification_datasets import SpecDatasets


class PlotMatrixPyvista(SpecDatasets):
	def __init__( self ):
		self.suffix = 'sigma_2'
		
	def run(self, non_rotated = False): 
		
		if hasattr( self, 'sub') and self.sub == True:
			filename  = 'matrix_pyvista_sub.png'
		else:
			filename  = 'matrix_pyvista_all.png'
		
		if hasattr( self, 'valencies') and ( self, 'lengths'):
			rows		= self.valencies
			columns		= self.lengths
			print('valencies and lengths')
		elif hasattr( self, 'stgs') and ( self, 'glun2bs'):
			rows		= self.glun2bs
			columns		= self.stgs
			print('glun2b and stg-conc dependence')
		else:
			print("Neither conc_dependence nor valnency_length.")
			sys.exit(1)
		
		num_rows		= len( rows )
		num_columns		= len( columns )
		
		
		pl = pyvista.Plotter(window_size=[1500,1500], \
			shape=(num_rows, num_columns),\
			border=False)
		
		for i, column in enumerate( columns ):
			for j, row in enumerate( rows ):
				# Load data
				filename_edited_prefix = self.filename_edited_matrix(row, column)
				d      = utils.load(self.dir_edited_data, filename_edited_prefix, self.suffix)
				print('Target: ', filename_edited_prefix)
				pl.subplot(num_rows-j-1, i)
				if non_rotated == True:
					utils_pyvista.plot_a_condensate_pyvista(d, pl)
				else:
					utils_pyvista.plot_a_pre_rotated_condensate_pyvista(d, pl)
				
				pl.camera_position = [(150.0, 0.0, 0.0),\
					 (7.0, -0.2, 2.5),\
					 (0.0, -1.0, 0.0)]
				
				#pl.camera.Zoom(1)
				
		pl.show(interactive=False, auto_close=True) # off_screen = True
		
		dir_imgs = os.path.join(self.dir_imgs_root, 'matrix_3d_pyvista')
		os.makedirs(dir_imgs, exist_ok=True)
		filename = os.path.join(dir_imgs, filename)
		pl.screenshot(filename)

	
	
if __name__ == '__main__':
	pass
	
