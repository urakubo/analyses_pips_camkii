
import os, glob, pickle, pprint, copy
import numpy as np

import utils
import pyvista



	
def save_a_plot(d, dir_img, prefix, suffix):
	pl = pyvista.Plotter(window_size=[400,1000], shape=(3, 1), border=False)
	#pl.add_text( '{}_{}'.format(prefix, suffix), position='lower_left', color='k', font='arial', font_size=10)
	
	pl.subplot(0, 0)
	utils.plot_a_condensate_pyvista(d, pl, rotation=False)
	pl.add_mesh(utils.square_yz(), color='black', style='wireframe')
	pl.view_yz()
	pl.camera.roll -= 90
	pl.camera.Zoom(1)
	
	pl.subplot(1, 0)
	utils.plot_a_condensate_pyvista(d, pl, rotation=False)
	pl.add_mesh(utils.square_zx(), color='black', style='wireframe')
	pl.view_zx()
	pl.camera.roll -= 90
	pl.camera.Zoom(1)
	
	pl.subplot(2, 0)
	utils.plot_a_condensate_pyvista(d, pl, rotation=False)
	pl.add_mesh(utils.square_xy(), color='black', style='wireframe')
	pl.view_xy()
	pl.camera.roll -= 90
	pl.camera.Zoom(1)
	
	#pl.show(interactive=True, auto_close=False) 
	pl.show(interactive=False, auto_close=True) # off_screen = True
	filename = os.path.join(dir_imgs, '{}_{}.png'.format(prefix, suffix))
	pl.screenshot(filename)
	
	
def save_plots_matrix(dir_data, dir_imgs, sigma): 
	#STG    = [1, 500, 1000, 1500, 2000, 2500, 3000] 
	#GluN2B = [1, 500, 1000, 2000, 4000, 6000, 8000, 12000, 16000, 20000]
	
	STG    = [540 , 1620, 2160,  2700,  3240, 4520] 
	GluN2B = [1080, 4320, 8640, 10800, 12960]
	
	suffix = 'sigma_{}'.format(sigma)
	
	pl = pyvista.Plotter(window_size=[1500,1500], shape=(len(GluN2B), len(STG)), border=False)
	
	for i, stg in enumerate(STG):
		for j, glun in enumerate(GluN2B):
			# Load data
			id = i + j * len(STG)
			prefix = str(id).zfill(3)
			d      = utils.load(dir_data, prefix, suffix)
			print('Target: {}, sigma: {}'.format(prefix, sigma))
			
			row    = i
			column = len(GluN2B)-j-1
			pl.subplot(column, row)
			#utils.plot_a_condensate_pyvista(d, pl)
			utils.plot_a_pre_rotated_condensate_pyvista(d, pl)

			if i == 0 and j == 0:
				pl.set_position( [100.0, 0.0, 0.0], reset=True )
				pl.view_yz()
				pl.camera.roll -= 90
				pos = pl.camera_position
				print('pl.camera_position', pl.camera_position)
				pl.camera_position =  [(150.0, 0.0, 0.0),\
					 (7.0, -0.2, 2.5),\
					 (0.0, -1.0, 0.0)]
			else:
				pl.camera_position = [(150.0, 0.0, 0.0),\
					 (7.0, -0.2, 2.5),\
					 (0.0, -1.0, 0.0)]
			
			#pl.camera.Zoom(1)
			
	pl.show(interactive=False, auto_close=True) # off_screen = True
	#pl.show(interactive=True, auto_close=False) # off_screen = True
	
	filename = os.path.join(dir_imgs, 'summary_{}.png'.format(suffix))
	pl.screenshot(filename)
	
	
if __name__ == '__main__':

	# Profiles
	'''
	# Files
	dir_target = 'conc_dependence'
	dir_edited_data 		= os.path.join('data2', dir_target)
	filenames_edited_data 	= [str(i).zfill(3) for i in [50, 52, 67] ]
	dir_imgs = os.path.join('imgs2', dir_target,'3d_condensate')
	os.makedirs(dir_imgs, exist_ok=True)
	
	sigma = 2
	for filename in filenames_edited_data:
		# Load data
		prefix = filename
		suffix = 'sigma_{}'.format(sigma)
		d      = utils.load(dir_edited_data, prefix, suffix)	
		print('Target: {}, sigma: {}'.format(filename, sigma))
		save_a_plot(d, dir_imgs, prefix, suffix)
	'''
	
	
	# Matrix
	#'''
	# Input files
	dir_target = 'conc_dependence_0.2'
	dir_edited_data 		= os.path.join('data2', dir_target)
	# Output files
	dir_imgs = os.path.join('imgs2', dir_target,'3d_condensate')
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	sigma    = 2 # 2, 3, or 4
	save_plots_matrix(dir_edited_data, dir_imgs, sigma)
	#'''
	

