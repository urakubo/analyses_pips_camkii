
import os, glob, pickle, pprint, copy
import numpy as np

import utils
import pyvista



def plot_a_image(d, pl, rotation=True): 
	
	flipz = False
	# Generate mesh
	r_CaMKII   = d['region_condensate_in_grid_mesh']['CaMKII'].astype(float)
	r_STG      = d['region_condensate_in_grid_mesh']['STG'].astype(float)
	mesh_CaMKII = utils.generate_mesh(r_CaMKII, flipz = flipz)
	mesh_STG    = utils.generate_mesh(r_STG   , flipz = flipz)
	
	if rotation == True:
		utils.rotate(mesh_CaMKII, mesh_STG)
	
	# Add cube
	'''
	cube  = pyvista.Cube(center=(0,0,0), \
		x_length=utils.space[0], y_length=utils.space[1], z_length=utils.space[2])
	pl.add_mesh(cube, color='black', style='wireframe')
	'''
	
	pl.add_mesh(mesh_CaMKII, color='green', show_edges=False,  opacity=0.4)
	pl.add_mesh(mesh_STG   , color='red', show_edges=False,  opacity=0.4)
	pl.set_background('white')
	
	
def save_a_plot(d, dir_img, prefix, suffix):
	pl = pyvista.Plotter(window_size=[300,800], shape=(3, 1), border=False)
	#pl.add_text( '{}_{}'.format(prefix, suffix), position='lower_left', color='k', font='arial', font_size=10)
	
	pl.subplot(0, 0)
	plot_a_image(d, pl, rotation=False)
	pl.add_mesh(utils.square_yz(), color='black', style='wireframe')
	pl.view_yz()
	pl.camera.roll -= 90
	pl.camera.Zoom(1)
	
	pl.subplot(1, 0)
	plot_a_image(d, pl, rotation=False)
	pl.add_mesh(utils.square_zx(), color='black', style='wireframe')
	pl.view_zx()
	pl.camera.roll -= 90
	pl.camera.Zoom(1)
	
	pl.subplot(2, 0)
	plot_a_image(d, pl, rotation=False)
	pl.add_mesh(utils.square_xy(), color='black', style='wireframe')
	pl.view_xy()
	pl.camera.roll -= 90
	pl.camera.Zoom(1)
	
	#pl.show(interactive=True, auto_close=False) 
	pl.show(interactive=False, auto_close=True) # off_screen = True
	filename = os.path.join(dir_imgs, '{}_{}.png'.format(prefix, suffix))
	pl.screenshot(filename)
	
	
def save_plots_matrix(dir_data, dir_imgs, sigma): 
	STG    = [500, 1000, 2000, 3000, 4000]
	GluN2B = [500, 2000, 4000, 6000, 8000, 12000]
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
			
			plot_a_image(d, pl)
			pl.view_yz()
			pl.camera.roll -= 90
			
			#plot_a_image(d, pl, rotation=False)
			#pl.view_xy()
			#pl.camera.roll -= 90

			pl.camera.Zoom(1)
			
	pl.show(interactive=False, auto_close=True) # off_screen = True
	filename = os.path.join(dir_imgs, 'summary_{}.png'.format(suffix))
	pl.screenshot(filename)
	
	
if __name__ == '__main__':

	# Target 1
	#'''
	# Input files
	dir_edited_data 		= 'data/conc_dependence'
	filenames_edited_data 	= [str(i).zfill(3) for i in range(3) ] # 30
	dir_imgs = os.path.join('imgs', 'conc_dependence','3d_condensate')
	sigma    = 2 # 2, 3, or 4
	
	os.makedirs(dir_imgs, exist_ok=True)
	
	sigmas = [2,4]
	for sigma in sigmas:
		for filename in filenames_edited_data:
			# Load data
			prefix = filename
			suffix = 'sigma_{}'.format(sigma)
			d      = utils.load(dir_edited_data, prefix, suffix)	
			print('Target: {}, sigma: {}'.format(filename, sigma))
			save_a_plot(d, dir_imgs, prefix, suffix)
	#'''
	
	
	# Target 2
	'''
	
	# Input files
	dir_edited_data 		= 'data/conc_dependence'
	filenames_edited_data 	= [str(i).zfill(3) for i in range(3) ] # 30
	# Output files
	dir_imgs = os.path.join('imgs', 'conc_dependence','3d_condensate')
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	sigma    = 2 # 2, 3, or 4
	save_plots_matrix(dir_edited_data, dir_imgs, sigma)
	'''
	

