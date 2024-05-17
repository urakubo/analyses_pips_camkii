
import os, glob, pickle, pprint, copy
import numpy as np

import lib.utils as utils
import lib.utils_pyvista as utils_pyvista


import pyvista



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
	filename = os.path.join(dir_imgs, '{}_{}.png'.format(prefix, suffix))
	pl.screenshot(filename)
	
	
	
if __name__ == '__main__':
	
	# Conc dependence
	'''
	dir_target = 'conc_dependence'
	filenames_edited_data 	= [str(i).zfill(3) for i in range(81) ]
	'''
	
	dir_target  = 'conc_dependence_merged'
	filenames_edited_data 	= [str(stg).zfill(2)+'_'+str(glun2b).zfill(2) for stg in range(10) for glun2b in range(9)  ]
	
	
	# Valency-length
	'''
	valencies = range(2,16,2)
	lengths   = range(7)
	filenames_edited_data = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in valencies for id_f in lengths ]
	dir_target  = 'valency_length'
	'''
	
	# 
	dir_edited_data 		= os.path.join('data4', dir_target)
	dir_imgs = os.path.join('imgs4', dir_target,'profiles_3d_pyvista')
	os.makedirs(dir_imgs, exist_ok=True)
	
	sigma = 2
	for filename in filenames_edited_data:
		# Load data
		prefix = filename
		suffix = 'sigma_{}'.format(sigma)
		d      = utils.load(dir_edited_data, prefix, suffix)	
		print('Target: {}, sigma: {}'.format(filename, sigma))
		save_a_plot(d, dir_imgs, prefix, suffix)
	
	
