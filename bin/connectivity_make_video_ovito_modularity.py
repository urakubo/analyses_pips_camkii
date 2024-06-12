
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
import lib.utils_ovito as utils_ovito
from specification_datasets import SpecDatasets




def plot_snapshots(data_all, sampling_frames, dir_imgs, fig_title, ids_col, ref_mc_step):
	
	
	modifier = utils_ovito.SetColorsById()
	modifier.ids_colors = ids_col
	data_all.modifiers.append(modifier)
	
	
	modifier = utils_ovito.DeleteParticle()
	modifier.ids = list(ids_col.keys())
	data_all.modifiers.append(modifier)
	
	
	data   = data_all.compute()
	data_all.add_to_scene()
	
	
	vp = Viewport()
	vp.type = Viewport.Type.Perspective
	#vp.fov = math.radians(40.0)
	#vp.camera_pos = (250+60, 60, 60)
	vp.fov = np.radians(10.0) # math.radians
	vp.camera_pos = (850, 60, 60)
	vp.camera_dir = (-1, 0, 0)
	
	
	i = 0
	for time_frame in sampling_frames:
		
		# Time label
		t =  ( data_all.compute(time_frame).attributes['Timestep'] - ref_mc_step ) / 1e9
		timelabel = TextLabelOverlay( text = 'Time: {:.5f} 10^9 MC step'.format(t),\
			font_size = 0.03, \
			text_color = (0,0,0) )
		timelabel.alignment = QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignBottom
		vp.overlays.append(timelabel)
		
		# Set filename
		filename = os.path.join(dir_imgs, fig_title+'_{}.png'.format( str(i).zfill(4)) )
		print( filename )
		i += 1
		
		# Plot
		vp.render_image(size=(800,800), filename=filename, background=(1,1,1),frame=time_frame, renderer=OpenGLRenderer())
		vp.overlays.remove(timelabel)
	
	return
	
	
def assign_color_for_each_bead(partitions, rcm, multi_graph_CaMKII):
	from cycler import cycler
	cols   = cycler( color=c.cols_list_ratio )
	colors = [c['color'] for p, c in zip(partitions, cols())]
	
	# ids_col = {id: c for ids_node, c in zip(rcm, colors) for id_node in ids_node for id in multi_graph_CaMKII.nodes[id_node]['id_bead_all']}
	ids_col = {}
	for ids_node, col in zip(rcm, colors):
		for id_node in ids_node:
			for id in multi_graph_CaMKII.nodes[id_node]['id_bead_all']:
				ids_col[id] = col
	return ids_col
	
	
	
	
class MakeOvitoVideoModularity(SpecDatasets):
		
	def __init__( self ):
		
		self.max_backward_frames_for_sampling = 80
		self.num_skip_frames_for_sampling     = 1
		
		
	def inspect( self ):
		for i, (f_lammpstrj, f_edited) in enumerate( zip(self.filenames_lammpstrj, self.filenames_edited) ):
			print('ID {}: {}, {}'.format(i, f_lammpstrj, f_edited))
		print()
		
		
		
	def run( self, i ):
		
		filename_lammpstrj = self.filenames_lammpstrj[i]
		filename_edited    = self.filenames_edited[i]
		print('ID {}: {}, {}'.format(i, filename_lammpstrj, filename_edited))
		
		
		dir_imgs         = os.path.join( self.dir_imgs_root,'connectivity_imgs_for_movie',filename_edited )
		os.makedirs(dir_imgs, exist_ok=True)
		
		
		# Load connect_graph.
		time_frame = 0
		prefix_connect_graph = '{}_{}'.format(time_frame, filename_edited)
		suffix_connect_graph ='louvain'
		
		d_connect_graph = utils.load(self.dir_edited_data, prefix_connect_graph, suffix_connect_graph)
		rcm                = d_connect_graph['rcm']
		multi_graph_CaMKII = d_connect_graph['multi_graph_CaMKII']
		partitions         = d_connect_graph['partitions']
		rcm                = [list(r) for r in rcm]
		ref_mc_step = d_connect_graph['mc_step']
		
		
		# Sampling frames
		num_frames = utils.get_num_frames(self.dir_lammpstrj, filename_lammpstrj)
		sampling_frames = []
		for i in range(0, self.max_backward_frames_for_sampling):
			step = num_frames + (i - self.max_backward_frames_for_sampling)*self.num_skip_frames_for_sampling
			sampling_frames.append( step )
		
		
		# Generate images
		ids_col  = assign_color_for_each_bead(partitions, rcm, multi_graph_CaMKII)
		data_all = import_file(os.path.join(self.dir_lammpstrj, filename_lammpstrj), input_format= "lammps/dump" )
		plot_snapshots(data_all, sampling_frames, dir_imgs, filename_edited, ids_col, ref_mc_step)
		
		
		
	def make_a_video( self, i ):
		
		filename_edited    = self.filenames_edited[i]
		dir_imgs         = os.path.join( self.dir_imgs_root,'connectivity_imgs_for_movie',filename_edited )
		dir_videos       = os.path.join( self.dir_imgs_root,'connectivity_movie' )
		os.makedirs(dir_videos, exist_ok=True)
		
		
		ffname = os.path.join(dir_imgs, '{}_%04d.png'.format(filename_edited))
		targ_name = os.path.join(dir_videos,'{}.mp4'.format(filename_edited) )
		com = ['ffmpeg','-r', '5', \
			'-i', ffname,\
			'-vf','scale=trunc(iw/2)*2:trunc(ih/2)*2',\
			'-pix_fmt', 'yuv420p', targ_name, '-y']
		print(' '.join(com))
		s.call(com)
		
		
if __name__ == '__main__':
	
	
	obj = MakeOvitoVideoModularity()
	obj.valency_length_small_colony2()
	obj.inspect()
	#i = 7*5+2 # val_12\R2_002
	i = 7*5+6 # val_12\R2_006
	obj.run( i = i )
	obj.make_a_video( i )
	
	
	