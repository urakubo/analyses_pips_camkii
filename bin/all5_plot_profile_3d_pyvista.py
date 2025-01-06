
import os, glob, pickle, pprint, copy
import numpy as np

import pyvista

import lib.utils as utils
import lib.utils_pyvista as utils_pyvista
from specification_datasets import SpecDatasets



def save_a_plot(d, dir_img, prefix, suffix):
	pl = pyvista.Plotter(window_size=[400,1000], shape=(3, 1), border=False)
	#pl.add_text( '{}_{}'.format(prefix, suffix), position='lower_left', color='k', font='arial', font_size=10)
	
	pl.subplot(0, 0)
	utils_pyvista.plot_a_condensate_pyvista(d, pl, rotation=False)
	pl.add_mesh(utils_pyvista.square_yz(), color='black', style='wireframe')
	pl.view_yz()
	pl.camera.roll -= 90
	pl.camera.Zoom(1)
	
	pl.subplot(1, 0)
	utils_pyvista.plot_a_condensate_pyvista(d, pl, rotation=False)
	pl.add_mesh(utils_pyvista.square_zx(), color='black', style='wireframe')
	pl.view_zx()
	pl.camera.roll -= 90
	pl.camera.Zoom(1)
	
	pl.subplot(2, 0)
	utils_pyvista.plot_a_condensate_pyvista(d, pl, rotation=False)
	pl.add_mesh(utils_pyvista.square_xy(), color='black', style='wireframe')
	pl.view_xy()
	pl.camera.roll -= 90
	pl.camera.Zoom(1)
	
	#pl.show(interactive=True, auto_close=False) 
	pl.show(interactive=False, auto_close=True) # off_screen = True
	filename = os.path.join(dir_img, '{}_{}.png'.format(prefix, suffix))
	pl.screenshot(filename)
	
	
class Plot3dPyvista(SpecDatasets):
	def __init__( self ):
		
		self.suffix = 'sigma_2'
		pass
		
	def run( self ):
		
		# Shared init
		self.dir_imgs = os.path.join(self.dir_imgs_root, 'profiles_3d_pyvista' )
		os.makedirs(self.dir_imgs, exist_ok=True)
		
		for filename in self.filenames_edited:
			self.plot_a_dataset( filename )
		
		
	def plot_a_dataset( self, filename ):
		
		# Load data
		print('Target: ', filename)
		prefix = filename

		d      = utils.load(self.dir_edited_data, prefix, self.suffix)
		save_a_plot(d, self.dir_imgs, prefix, self.suffix)
	
	
if __name__ == '__main__':
	
	obj = Plot3dPyvista()
	obj.conc_dependence_merged() #  conc_dependence(), valency_length(), valency_length_CG(), boundary_conditions2()
	obj.run()
	
	
