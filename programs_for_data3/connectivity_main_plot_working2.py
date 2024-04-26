
# https://stackoverflow.com/questions/59821151/plot-the-dendrogram-of-communities-found-by-networkx-girvan-newman-algorithm
# https://qiita.com/s-wakaba/items/a93f03f27137cff4a26c

import networkx as nx
from itertools import chain, combinations
import matplotlib.pyplot as plt
from matplotlib import patches

from scipy.cluster.hierarchy import dendrogram
import os, sys

import numpy as np

import utils
import parameters as p
import colormap as c

plt.rcParams.update(p.rc_param)

import os, glob

import numpy as np
import matplotlib.pyplot as plt

import utils
import colormap as c
import parameters as p

import pyvista



sphere = pyvista.Sphere(radius=2.0, phi_resolution=10, theta_resolution=10)
sphere_big = pyvista.Sphere(radius=20.0, center=(0, 20, 0), phi_resolution=40, theta_resolution=40)

def plot_pyvista(pl, bounds0, bounds1, xyzlabel, t, type, position, i, sphere_color, sphere_opacity): 

	pl.subplot(0, 0)
	text = "Frame:"+str(i).rjust(4)
	print(text)
	a0     = pl.add_text(text, position='lower_left', color='k', font='arial', font_size=10)
	actors = [a0]

	for _t in t.values():
		p_target     = position[ type == _t['id']   , :]
		p_foreground = p_target[ p_target[:,0] >= 0 , :]
		p_background = p_target[ p_target[:,0] <  0 , :]
		
		#
		pl.subplot(0, 0)
		a = pl.add_mesh(sphere_big, color=sphere_color, opacity=sphere_opacity)
		pl.subplot(0, 1)
		b = pl.add_mesh(sphere_big, color=sphere_color, opacity=sphere_opacity)
		actors.append([a, b])
		# print('p_background: ', p_background)
		
		if len(p_background) != 0:
			p_v_background = pyvista.PolyData(p_background).glyph(scale=False, geom=sphere, orient=False)
			pl.subplot(0, 0)
			a2 = pl.add_mesh(p_v_background, color=_t['c'])
			pl.subplot(0, 1)
			a3 = pl.add_mesh(p_v_background, color=_t['c'])
			actors.append([a2, a3])
		
		if len(p_foreground) != 0:
			p_v_foreground = pyvista.PolyData(p_foreground).glyph(scale=False, geom=sphere, orient=False)
			pl.subplot(0, 0)
			a1 = pl.add_mesh(p_v_foreground, color=_t['c'], opacity=0.04)
			actors.append([a1])

	pl.subplot(0, 0)
	a4 = pl.show_grid(bounds=bounds0, color='black', grid=True, xlabel=xyzlabel, ylabel=xyzlabel, zlabel=xyzlabel)
	actors.append(a4)

	pl.subplot(0, 1)
	a5 = pl.show_grid(bounds=bounds1, color='black', grid=True, xlabel=xyzlabel, ylabel=xyzlabel, zlabel=xyzlabel)
	actors.append([a4, a5])
	return actors




def plot_snapshots(data_all, num_frames, dir_imgs, fig_title, ids_untarget_molecules):

	rect     = utils.create_rectangle_yz_for_pyvista([m,m])
	xyzlabel = 'Distance from the center (grids)'
	
	

	#pl = pyvista.Plotter(window_size=[1000,500], shape=(1, 2), border=False, off_screen=True)
	pl = pyvista.Plotter(window_size=[1000,500], shape=(1, 2), border=False)
	pl.background_color = 'w'
	pl.subplot(0, 0)
	pl.add_mesh(rect, color='red', opacity=0.2)
	pl.camera.Zoom(0.5)	
	
	
	
	for i, (type, position) in enumerate(zip(types, positions)):
		
		# Plot data
		actors = utils.plot_pyvista(pl, bounds0, bounds1, xyzlabel, table_molecules, \
				type, position, time_frame[i], sphere_color, sphere_opacity)
		if i == 0:
			pl.subplot(0, 1)
			pl.view_yz()
			camera = pl.camera.position
			f = False
		else:
			pl.camera.position = camera
		filename = os.path.join(dir_imgs, str(i).zfill(4)+ '.png')
		pl.show(interactive=False, auto_close=False) # Comment it out if off_screen = True.
		pl.screenshot(filename)
		
		# Remove plots
		for a in actors:
			pl.remove_actor(a)



if __name__ == '__main__':

	## Init
	dir_target  =  'small_colony'
	prefix = '01_008'
	#prefix = '01_004'

	# lengths_clusters  [1503, 881, 699, 447, 274, 1, 1, 1, 1, 1]
	prefix = '01_004'
	nth_largest = 1

	# lengths_clusters  [845, 838, 793, 443, 372, 368, 1, 1, 1, 1]
	prefix = '00_004'
	nth_largest = 0 #0


	fig_title = '{}_{}'.format(nth_largest, prefix)
	print(fig_title)
	
	
	
	subdirs    = ['CG_con', 'CG_len9', 'CG_lin']
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(10)]
	filenames_lammpstrj = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_edited    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(3) for id_f in range(10)]
	
	filename_lammpstrj  = filenames_lammpstrj[ filenames_edited.index(prefix) ]
	
	
	
	
	dir_lammpstrj    = os.path.join('..','lammpstrj3', dir_target)
	dir_edited_data  = os.path.join('data3', dir_target)
	dir_imgs         = os.path.join('imgs3', dir_target, 'connectivity_movie')
	os.makedirs(dir_imgs, exist_ok=True)

	# Load graph.
	d2 = utils.load(dir_edited_data, prefix, 'cluster_dendrogram')
	
	Z          = d2['Z']
	labels     = d2['labels']
	node_order = d2['node_order']
	blocks     = d2['blocks'] 
	partition  = d2['partition'] 
	leaves_color_list  = d2['leaves_color_list'] 
	multi_graph_CaMKII = d2['multi_graph_CaMKII']
	
	# Load lammpstrj data
	num_frames = utils.get_num_frames(dir_lammpstrj, filename_lammpstrj)
	
	print('num_frames ', num_frames)
	
	
	
	data_all   = import_file(os.path.join(dir_lammpstrj, filename_lammpstrj), input_format= "lammps/dump" )
	
	plot_snapshots(data_all, num_frames, dir_imgs, fig_title, ids_untarget_molecules)
	
	