
import os, glob, pickle, pprint, copy
import numpy as np

import parameters as p
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
	
	
def save_plots_matrix(dir_data, dir_imgs, sub = True): 
	
	STG    = [540, 1620, 2160, 2700, 3240, 4520] 
	GluN2B = [570, 1080, 4320, 6480, 8640, 10800, 12960, 17280]
	
	sub_STG    = [540, 1620, 2160, 2700, 3240] 
	sub_GluN2B = [1080, 4320, 6480, 10800, 17280]
	
	'''
	volume = np.prod(p.space_np)
	STG    = [ s / volume for s in STG    ]
	GluN2B = [ n / volume for n in GluN2B ]
	'''
	
	if sub == True:
		sampling_STG    = {num: i for i, num in enumerate(sub_STG)}
		sampling_GluN2B = {num: i for i, num in enumerate(sub_GluN2B)}
		filename = 'summary_3d_sub.png'
	else:
		sampling_STG    = {num: i for i, num in enumerate(STG)}
		sampling_GluN2B = {num: i for i, num in enumerate(GluN2B)}
		filename = 'summary_3d_all.png'
	
	
	pl = pyvista.Plotter(window_size=[1500,1500], \
		shape=(len(sampling_GluN2B), len(sampling_STG)),\
		border=False)
	
	for i, stg in enumerate(STG):
		for j, glun in enumerate(GluN2B):
			# Load data
			if (stg in sampling_STG.keys()) and (glun in sampling_GluN2B.keys()):
				# Load data
				id = i + j * len(STG)
				prefix = str(id).zfill(3)
				suffix = 'sigma_2'
				d      = utils.load(dir_data, prefix, suffix)
				print('Target: {}'.format(prefix))
				sub_i = sampling_STG[stg]
				sub_j = sampling_GluN2B[glun]
				row    = sub_i
				column = len(sampling_GluN2B)-sub_j-1
				pl.subplot(column, row)
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


	# Matrix
	'''
	# Input files
	dir_target = 'conc_dependence'
	dir_edited_data = os.path.join('data3', dir_target)
	# Output files
	dir_imgs = os.path.join('imgs3', dir_target,'3d_condensate')
	os.makedirs(dir_imgs, exist_ok=True)
	
	save_plots_matrix(dir_edited_data, dir_imgs)
	'''
	

