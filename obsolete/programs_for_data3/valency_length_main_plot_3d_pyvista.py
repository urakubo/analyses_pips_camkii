
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
	
	
def save_plots_matrix(dir_data, dir_imgs, sigma, sub = True): 
	
	
	valency = list(range(2,14,2)) 
	linker_length  = [1, 2, 3, 4, 5, 6, 9]
	
	if sub == True:
		sub_valency = [2, 4, 6, 8, 12]
		sub_length  = [1, 3, 4, 5, 6]
		filename = 'summary_3d_sub.png'
	else:
		sub_valency = valency
		sub_length  = linker_length
		filename = 'summary_3d_all.png'
	
	fnames_valency       = { v: str(v).zfill(2) for v in valency if v in sub_valency}
	fnames_linker_length = {ll: str(i).zfill(3) for i,ll in enumerate(linker_length) if ll in sub_length}
	
	
	
	num_rows		= len( sub_valency )
	num_columns		= len( sub_length )
	
	suffix = 'sigma_{}'.format(sigma)
	
	pl = pyvista.Plotter(window_size=[1500,1500], shape=(num_rows, num_columns), border=False)
	
	
	for i, v in enumerate(sub_valency):
		for j, ll in enumerate(sub_length):
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

			pl.camera_position = [(150.0, 0.0, 0.0),\
				 (7.0, -0.2, 2.5),\
				 (0.0, -1.0, 0.0)]
			
			#pl.camera.Zoom(1)
			
	pl.show(interactive=False, auto_close=True) # off_screen = True
	#pl.show(interactive=True, auto_close=False) # off_screen = True
	
	filename = os.path.join(dir_imgs, filename)
	pl.screenshot(filename)
	
	
if __name__ == '__main__':
	
	# File
	dir_target = 'valency_length'
	
	
	dir_edited_data	= os.path.join('data3', dir_target)
	dir_imgs = os.path.join('imgs3', dir_target,'3d_condensate')
	os.makedirs(dir_imgs, exist_ok=True)
	sigma = 2
	
	
	# Plot the matrix
	#'''
	save_plots_matrix(dir_edited_data, dir_imgs, sigma)
	#'''
	

