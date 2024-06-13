

import os, sys, glob, pprint,itertools, shutil
import numpy as np
os.environ['OVITO_GUI_MODE'] = '1' # Request a session with OpenGL support

#import networkx as nx
#import matplotlib.pyplot as plt

import subprocess as s


from ovito.io import import_file
from ovito.modifiers import SliceModifier, SelectTypeModifier, AssignColorModifier
from ovito.vis import Viewport, TextLabelOverlay, OpenGLRenderer
from ovito.pipeline import Pipeline, StaticSource
from ovito.qt_compat import QtCore
from ovito.data import DataCollection, ParticleType


import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph
import lib.utils_ovito as utils_ovito
from specification_datasets import SpecDatasets



def plot_snapshots(data_all, dir_imgs, \
		frame_num_after_time_zero  = 100, \
		frame_num_before_time_zero = 10, \
		num_skip_frames_for_sampling = 1 \
		):
	
	
	# Time frames
	num_frames = data_all.source.num_frames
	
	target_frames     = range(num_frames-frame_num_after_time_zero-frame_num_before_time_zero,\
						num_frames,\
						num_skip_frames_for_sampling)
	target_frames     = list(target_frames)
	frame_time_zero = num_frames - frame_num_after_time_zero
	time_time_zero  = data_all.compute(frame_time_zero).attributes['Timestep'] / 1e9
	target_times    = np.array([data_all.compute(t).attributes['Timestep'] for t in target_frames])/ 1e9 - time_time_zero
	
	
	print( 'target_frames' )
	print( target_frames )
	
	# label colors
	for k, v in p.molecules_without_all.items():
		data_all.modifiers.append(SelectTypeModifier(types=set(v['id'])))
		data_all.modifiers.append(AssignColorModifier(color=c.cmap_universal_ratio[k] ))
	
	
	# Centering
	types, positions, ids_molecule = utils.decode_data(data_all.compute(frame_time_zero))
	center = utils.get_center_of_mass(types, positions)
	for dim in [0,1,2]:
		center[dim] += - p.space[dim] * (center[dim]  >=  p.space[dim]) + p.space[dim] * (center[dim]  <  0)
	
	modifier = utils_ovito.CenteringModifier()
	modifier.center = center
	data_all.modifiers.append(modifier)
	
	
	# Slice data
	modifier = SliceModifier()
	modifier.normal   = (1.0, 0.0, 0.0)
	modifier.distance = 60
	data_all.modifiers.append(modifier)
	
	
	vp = Viewport()
	vp.type = Viewport.Type.Perspective
	vp.fov = np.radians(10.0) # math.radians
	vp.camera_pos = (850, 60, 60)
	vp.camera_dir = (-1, 0, 0)
	
	#data   = data_all.compute()
	data_all.add_to_scene()
	
	i = 0
	for target_frame, target_time in zip(target_frames, target_times):
		print('Time frame: ', target_frame, ', Time: ', target_time)
		
		# Set time
		timelabel = TextLabelOverlay( text = 'Time: {:.0f} G MC steps'.format(target_time),\
			font_size = 0.03, \
			text_color = (0,0,0) )
		timelabel.alignment = QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignBottom
		vp.overlays.append(timelabel)
		
		filename = os.path.join(dir_imgs, '{}.png'.format( str(i).zfill(4)) )
		print( filename )
		i += 1
		vp.render_image(size=(800,800), filename=filename, background=(1,1,1),frame=target_frame, renderer=OpenGLRenderer())
		vp.overlays.remove(timelabel)
	
	data_all.remove_from_scene()
	return
	
	
class MakeOvitoVideo(SpecDatasets):
		
	def __init__( self ):
		
		self.frame_num_after_time_zero  = 100
		self.frame_num_before_time_zero = 10
		self.num_skip_frames_for_sampling = 1
		
		
	def inspect( self ):
		for i, (f_lammpstrj, f_edited) in enumerate( zip(self.filenames_lammpstrj, self.filenames_edited) ):
			print('ID {}: {}, {}'.format(i, f_lammpstrj, f_edited))
		
		
	def run( self, i ):
		
		print('\n ID {}: {}, {}'.format(i, self.filenames_lammpstrj[i], self.filenames_edited[i]))
		# Shared init
		
		dir_imgs         = os.path.join( self.dir_imgs_root, 'for_movie_{}'.format( self.filenames_edited[i] ) )
		os.makedirs(dir_imgs, exist_ok=True)
		
		# Load lammpstrj file.
		data_all = import_file(os.path.join(self.dir_lammpstrj, self.filenames_lammpstrj[i]), input_format= "lammps/dump" )
		
		# Make images
		plot_snapshots(data_all, dir_imgs, \
			self.frame_num_after_time_zero,  
			self.frame_num_before_time_zero, 
			self.num_skip_frames_for_sampling
			)
		
		
	def make_a_video( self, i ):
		
		#
		dir_videos    = os.path.join( self.dir_imgs_root, 'movies' )
		dir_imgs = os.path.join( self.dir_imgs_root, 'for_movie_{}'.format( self.filenames_edited[i] ) )
		os.makedirs(dir_videos, exist_ok=True)
		
		ffname = os.path.join( dir_imgs, '%04d.png' )
		targ_name = os.path.join(dir_videos,'{}.mp4'.format( self.filenames_edited[i] ) )
		com = ['ffmpeg','-r', '5', \
			'-i', ffname,\
			'-vf','scale=trunc(iw/2)*2:trunc(ih/2)*2',\
			'-pix_fmt', 'yuv420p', targ_name, '-y']
		print(' '.join(com))
		s.call(com)
		
		
if __name__ == '__main__':
	
	obj = MakeOvitoVideo()
	obj.boundary_conditions2()
	obj.inspect()
	obj.run(2)
	obj.make_a_video(2)
	
	# Target: SPG
	#filename_lammpstrj = 'PIPS_activated_trj.lammpstrj'
	#filename_video     = 'PIPS'
	#dir_lammpstrj_sub  = 'SPG'
	#dir_target         = 'boundary_conditions'
	
	
