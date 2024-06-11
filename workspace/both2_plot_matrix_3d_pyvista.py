
import os, glob, pickle, pprint, copy
import numpy as np

import lib.utils as utils
import lib.parameters as p
import lib.utils_pyvista as utils_pyvista


import pyvista


def plot_valency_length_save_it(dir_data, dir_imgs, sub = False): 
	
	
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
	
	fnames_valency       = {v: str(v).zfill(2) for v in valency if v in sub_valency}
	fnames_linker_length = {l: str(i).zfill(3) for i,l in enumerate(linker_length) if l in sub_length}
	
	
	
	num_rows		= len( sub_valency )
	num_columns		= len( sub_length )
	
	suffix = 'sigma_2'
	
	pl = pyvista.Plotter(window_size=[1500,1500], shape=(num_rows, num_columns), border=False)
	
	
	for i, v in enumerate(sub_valency):
		for j, ll in enumerate(sub_length):
			# Load data
			prefix = fnames_valency[v]+'_'+fnames_linker_length[ll]
			d      = utils.load(dir_edited_data, prefix, suffix)
			print('Target: {}'.format(prefix))
			
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
	#pl.show(interactive=True, auto_close=False) # off_screen = True
	
	filename = os.path.join(dir_imgs, filename)
	pl.screenshot(filename)


	
def plot_conc_dependence_save_it(dir_data, dir_imgs, sub = False): 
	
	#STGs    = [108,216,432,576,864,1728,2592,3456,4320,5184]
	#GluN2Bs = [270,540,1080,2160,4320,6480,8640,12960,17280]
	
	
	# Subset of concs
	sub_STGs    = [216, 576, 1728, 3456] 
	sub_GluN2Bs = [540, 4320, 8640,17280]
	volume      = np.prod(p.space_np)
	sub_STGs    = [ s / volume for s in sub_STGs    ]
	sub_GluN2Bs = [ n / volume for n in sub_GluN2Bs ]
	suffix = 'sigma_2'
	
	if sub == True:
		sampling_STG    = {i: STG    for i, STG    in enumerate(p.STGs)    if STG in sub_STGs}
		sampling_GluN2B = {i: GluN2B for i, GluN2B in enumerate(p.GluN2Bs) if GluN2B in sub_GluN2Bs}
		save_img_filename = 'summary_3d_sub.png'
	else:
		sampling_STG    = {i: STG    for i, STG    in enumerate(p.STGs)}
		sampling_GluN2B = {i: GluN2B for i, GluN2B in enumerate(p.GluN2Bs)}
		save_img_filename = 'summary_3d_all.png'
	
	
	pl = pyvista.Plotter(window_size=[1500,1500], \
		shape=(len(sampling_GluN2B), len(sampling_STG)),\
		border=False)
	
	for i_fig, (i_file, stg) in enumerate(sampling_STG.items()):
		for j_fig, (j_file, glun) in enumerate(sampling_GluN2B.items()):
			# Load data
			prefix = str(i_file).zfill(2)+'_'+str(j_file).zfill(2)
			d      = utils.load(dir_data, prefix, suffix)
			print('Target: {}'.format(prefix))
			row    = i_fig
			column = len(sampling_GluN2B)-j_fig-1
			pl.subplot(column, row)
			#utils.plot_a_condensate_pyvista(d, pl)
			utils_pyvista.plot_a_pre_rotated_condensate_pyvista(d, pl)
			
			pl.camera_position = [(150.0, 0.0, 0.0),\
				 (7.0, -0.2, 2.5),\
				 (0.0, -1.0, 0.0)]
			
			#pl.camera.Zoom(1)
			
	pl.show(interactive=False, auto_close=True) # off_screen = True
	#pl.show(interactive=True, auto_close=False) # off_screen = True
	
	save_img_filename = os.path.join(dir_imgs, save_img_filename)
	pl.screenshot(save_img_filename)
	
	
	
if __name__ == '__main__':
	
	# Valency length
	#'''
	#dir_target = 'valency_length'
	dir_target = 'CG_valency_length'
	func = plot_valency_length_save_it
	#'''
	
	'''
	dir_target = 'conc_dependence_merged'
	func       = plot_conc_dependence_save_it
	'''
	sub        = False # False
	
	dir_edited_data	= os.path.join('data4', dir_target)
	dir_imgs = os.path.join('imgs4', dir_target, 'matrix_3d_pyvista')
	os.makedirs(dir_imgs, exist_ok=True)
	func(dir_edited_data, dir_imgs, sub = sub)
	
