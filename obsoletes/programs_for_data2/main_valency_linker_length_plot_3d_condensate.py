
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
	
	
	valency = list(range(2,14,2)) 
	linker_length  = [1, 3, 5, 7, 9]
	
	fnames_valency       = { v: str(v).zfill(2) for v in valency }
	fnames_linker_length = {ll: str(i).zfill(3) for i,ll in enumerate(linker_length) }
	
	num_rows		= len( valency )
	num_columns		= len( linker_length )
	
	suffix = 'sigma_{}'.format(sigma)
	
	pl = pyvista.Plotter(window_size=[1500,1500], shape=(num_rows, num_columns), border=False)
	
	
	for i, v in enumerate(valency):
		for j, ll in enumerate(linker_length):
			# Load data
			prefix = fnames_valency[v]+'_'+fnames_linker_length[ll]
			suffix = 'sigma_{}'.format(sigma)
			d      = utils.load(dir_edited_data, prefix, suffix)
			print('Target: {}, sigma: {}'.format(prefix, sigma))
			
			# Specify row and column
			row    = num_rows - i - 1
			column = j
			
			pl.subplot(row, column)
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
	
	# File
	dir_target = 'valency_linker_length'
	
	
	dir_edited_data	= os.path.join('data2', dir_target)
	dir_imgs = os.path.join('imgs2', dir_target,'3d_condensate')
	os.makedirs(dir_imgs, exist_ok=True)
	sigma = 2
	
	
	'''
	# Plot the snapshots of each result from the three directions.
	filenames_edited_data = [str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(5) ]
	
	for filename in filenames_edited_data:
		# Load data
		prefix = filename
		suffix = 'sigma_{}'.format(sigma)
		d      = utils.load(dir_edited_data, prefix, suffix)	
		print('Target: {}, sigma: {}'.format(filename, sigma))
		save_a_plot(d, dir_imgs, prefix, suffix)
	'''
	
	
	# Plot the matrix
	#'''
	save_plots_matrix(dir_edited_data, dir_imgs, sigma)
	#'''
	

