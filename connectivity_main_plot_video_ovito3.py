
# https://stackoverflow.com/questions/59821151/plot-the-dendrogram-of-communities-found-by-networkx-girvan-newman-algorithm
# https://qiita.com/s-wakaba/items/a93f03f27137cff4a26c



import os, sys, glob, pprint,itertools
import numpy as np
os.environ['OVITO_GUI_MODE'] = '1' # Request a session with OpenGL support

import networkx as nx
import matplotlib.pyplot as plt


import subprocess as s


from ovito.io import import_file
from ovito.modifiers import SliceModifier, SelectTypeModifier, AssignColorModifier
from ovito.pipeline import ModifierInterface
from ovito.vis import *
from ovito.pipeline import *
from ovito.qt_compat import QtCore


import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
#import lib.utils_graph as utils_graph


def delete_particle(frame, data, ids):
	transparency = np.ones((data.particles.count), dtype='int')
	transparency[ids] = 0

class DeleteParticle(ModifierInterface):
	def modify(self, data, frame, **kwargs):
		transparency = np.ones((data.particles.count), dtype='int')
		transparency[self.ids] = 0
		data.particles_.delete_elements(transparency)
		#print('transparency ', transparency)

class MakeTransparent(ModifierInterface):
	def modify(self, data, frame, **kwargs):
		transparency = np.ones((data.particles.count), dtype='float')*0.95
		transparency[self.ids] = 0.0
		data.particles_.create_property('Transparency', data=transparency)

class SetColorsById(ModifierInterface):
	def modify(self, data, frame, **kwargs):
		color_values = np.ones((data.particles.count, 3), dtype='float')
		for id, col in self.ids_colors.items():
			color_values[id,:] = col
		data.particles_.create_property('Color', data=color_values)
	
	
def plot_snapshots(data_all, dir_imgs, fig_title):
	
	'''
	for k, v in p.molecules_without_all.items():
		data_all.modifiers.append(SelectTypeModifier(types=set(v['id'])))
		data_all.modifiers.append(AssignColorModifier(color=c.cmap_universal_ratio[k] ))
	'''
	
	#data_all.modifiers.append(compute_particle_transparency)
	
	

	vp = Viewport()
	vp.type = Viewport.Type.Perspective
	vp.fov = np.radians(10.0) # math.radians
	vp.camera_pos = (850, 60, 60)
	vp.camera_dir = (-1, 0, 0)
	
	
	# Calc colors
	from cycler import cycler
	cols_matplotlib = plt.rcParams['axes.prop_cycle'].by_key()['color']
	cols = [get_ratio_code(hex_code) for hex_code in cols_matplotlib]
	
	# Intial
	rcm, multi_graph, ref_mc_step = load_graph_data(dir_edited_data, filename_edited, 0)
	ids_col  = {id: col for ids_node, col in zip(rcm, cols) for id_node in ids_node for id in multi_graph.nodes[id_node]['id_bead_all']}
	rcm_prev = rcm
	
	max_backward_frames_for_sampling = 40
	num_skip_frames_for_sampling     = 1
	
	num_frames = utils.get_num_frames(dir_lammpstrj, filename_lammpstrj)
	
	data   = data_all.compute()
	data_all.add_to_scene()
	
	for i in range(0, max_backward_frames_for_sampling):
		
		sampling_frame   = num_frames + (i - max_backward_frames_for_sampling)*num_skip_frames_for_sampling
		rcm, multi_graph, mc_step = load_graph_data(dir_edited_data, filename_edited, i )
		t =  ( mc_step - ref_mc_step ) / 1e9
		
		ids_col = {}
		cols_next = []
		for i_current in rcm:
			shared_hubs = [len(set(i_prev) & set(i_current)) for i_prev in rcm_prev]
			max_index   = shared_hubs.index(max(shared_hubs))
			ids_col.update( {id: cols[max_index] for id_node in i_current for id in multi_graph.nodes[id_node]['id_bead_all']} )
			cols_next.append(cols[max_index])
		
		cols = cols_next
		rcm_prev = rcm
		
		modifier = SetColorsById()
		modifier.ids_colors = ids_col
		data_all.modifiers.append(modifier)
		
		modifier_delete = DeleteParticle()
		modifier_delete.ids = list(ids_col.keys())
		data_all.modifiers.append(modifier_delete)
		
		
		timelabel = TextLabelOverlay( text = 'Time: {:.5f} G MC step'.format(t),\
			font_size = 0.03, \
			text_color = (0,0,0) )
		timelabel.alignment = QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignBottom
		vp.overlays.append(timelabel)
		
		filename = os.path.join(dir_imgs, fig_title+'__{}.png'.format( str(i).zfill(4)) )
		print( filename )
		i += 1
		vp.render_image(size=(800,800), filename=filename, background=(1,1,1),frame=sampling_frame, renderer=OpenGLRenderer())
		vp.overlays.remove(timelabel)
		del data_all.modifiers[-1]
		# print('dir(data_all) ', dir(data_all))
		# data_all.delete(modifier_delete)
	
	
	return


def load_graph_data(dir_edited_data, filename_edited, i ):
	prefix_connect_graph = '{}_{}'.format(i, filename_edited)
	suffix_connect_graph ='louvain'
	d_connect_graph = utils.load(dir_edited_data, prefix_connect_graph, suffix_connect_graph)
	rcm             = d_connect_graph['rcm']
	multi_graph     = d_connect_graph['multi_graph_CaMKII']
	rcm             = [list(r) for r in rcm]
	mc_step         = d_connect_graph['mc_step']
	return rcm, multi_graph, mc_step


def get_ratio_code(hex_code):
	hex_code   = hex_code.lstrip("#")
	ratio_code = [int(hex_code[i:i+2], 16)/255 for i in range(0, 6, 2)]
	return ratio_code

if __name__ == '__main__':
	
	
	# Input file
	
	# Small colony 2
	subdirs    = ['val_{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
	filenames_lammpstrj = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_edited    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7)]
	dir_target  = 'small_colony2'
	
	
	i = 7*5+2 # val_12\R2_002
	#i = 7*4+2 # val_10\R2_002
	#i = 7*3+2 # val_08\R2_002
	#i = 7*2+2 # val_06\R2_002
	filename_lammpstrj = filenames_lammpstrj[i]
	filename_edited    = filenames_edited[i]
	
	
	
	
	## Shared init
	#fig_title = '{}_{}'.format(nth_largest, prefix)
	#print(fig_title)
	
	dir_lammpstrj    = os.path.join('..','lammpstrj4', dir_target)
	dir_edited_data  = os.path.join('data4', dir_target)
	dir_imgs         = os.path.join('imgs4', dir_target, 'connectivity_imgs_for_movie')
	dir_videos       = os.path.join('imgs4', dir_target, 'connectivity_movie')
	os.makedirs(dir_imgs, exist_ok=True)
	os.makedirs(dir_videos, exist_ok=True)
	
	
	data_all = import_file(os.path.join(dir_lammpstrj, filename_lammpstrj), input_format= "lammps/dump" )
	plot_snapshots(data_all, dir_imgs, filename_edited)
	
	
	ffname = os.path.join(dir_imgs, '{}__%04d.png'.format(filename_edited))
	targ_name = os.path.join(dir_videos,'{}.mp4'.format(filename_edited) )
	com = ['ffmpeg','-r', '5', \
		'-i', ffname,\
		'-vf','scale=trunc(iw/2)*2:trunc(ih/2)*2',\
		'-pix_fmt', 'yuv420p', targ_name, '-y']
	print(' '.join(com))
	s.call(com)

	
	