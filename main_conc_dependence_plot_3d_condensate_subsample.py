
import os, glob, pickle, pprint, copy
import numpy as np

import utils
import pyvista


	
	
def save_plots_matrix(dir_data, dir_imgs, sigma=2): 
	STG    = [1, 500, 1000, 1500, 2000, 2500, 3000] 
	GluN2B = [1, 500, 1000, 2000, 4000, 6000, 8000, 12000, 16000, 20000]
	
	sub_STG    = [500, 1000, 1500, 2000, 2500, 3000] 
	sub_GluN2B = [1000, 4000, 8000, 12000, 16000, 20000]
	
	sub_STG    = {num: i for i, num in enumerate(sub_STG)}
	sub_GluN2B = {num: i for i, num in enumerate(sub_GluN2B)}
	
	suffix = 'sigma_{}'.format(sigma)
	pl = pyvista.Plotter(window_size=[1500,1500], shape=(len(sub_GluN2B), len(sub_STG)), border=False)
	
	for i, stg in enumerate(STG):
		for j, glun in enumerate(GluN2B):
			# Load data
			if (stg in sub_STG.keys()) and (glun in sub_GluN2B.keys()):
				id = i + j * len(STG)
				prefix = str(id).zfill(3)
				d      = utils.load(dir_data, prefix, suffix)
				
				sub_i = sub_STG[stg]
				sub_j = sub_GluN2B[glun]
				
				row    = sub_i
				column = len(sub_GluN2B)-sub_j-1
				pl.subplot(column, row)
				utils.plot_a_pre_rotated_condensate_pyvista(d, pl)
				'''
				if i == 0 and j == 0:
					pl.view_yz()
					pl.camera.roll -= 90
				'''
				pl.camera_position =  [(150.0, 0.0, 0.0),\
					 (7.0, -0.2, 2.5),\
					 (0.0, -1.0, 0.0)]
	
	pl.show(interactive=False, auto_close=True)
	#pl.show(interactive=True, auto_close=False) # off_screen = True
	
	filename = os.path.join(dir_imgs, 'subsample_{}.png'.format(suffix))
	pl.screenshot(filename)
	
	
if __name__ == '__main__':

	# Profiles
	
	# Input files
	dir_edited_data 		= os.path.join('data', 'conc_dependence')
	filenames_edited_data 	= [str(i).zfill(3) for i in range(3) ] # 30
	# Output files
	dir_imgs = os.path.join('imgs', 'conc_dependence','3d_condensate')
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	save_plots_matrix(dir_edited_data, dir_imgs)
	

