

import os, sys, glob, pprint,itertools, shutil
import numpy as np
os.environ['OVITO_GUI_MODE'] = '1' # Request a session with OpenGL support

import networkx as nx
import matplotlib.pyplot as plt


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



		

def plot_snapshots(data_all, dir_imgs, filename_edited, \
		frame_num_after_photobleach  = 300, \
		frame_num_before_photobleach = 10, \
		num_skip_frames_for_sampling = 1 \
		):
	
	
	# Time frames
	num_frames = data_all.source.num_frames
	
	target_frames     = range(num_frames-frame_num_after_photobleach-frame_num_before_photobleach,\
						num_frames,\
						num_skip_frames_for_sampling)
	target_frames     = list(target_frames)
	frame_photobleach = num_frames - frame_num_after_photobleach
	time_photobleach  = data_all.compute(frame_photobleach).attributes['Timestep'] / 1e9
	time_steps        = np.array([data_all.compute(t).attributes['Timestep'] for t in target_frames])/ 1e9 - time_photobleach
	
	# label colors
	for k, v in p.molecules_without_all.items():
		data_all.modifiers.append(SelectTypeModifier(types=set(v['id'])))
		data_all.modifiers.append(AssignColorModifier(color=c.cmap_universal_ratio[k] ))
		
	
	
	# Centering
	types, positions, ids_molecule = utils.decode_data(data_all.compute(num_frames))
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
	for target_frame, target_time in zip(target_frames, time_steps):
		print('Time frame: ', target_frame, ', Time: ', target_time)
		
		# Set time
		timelabel = TextLabelOverlay( text = 'Time: {:.0f} G MC steps'.format(target_time),\
			font_size = 0.03, \
			text_color = (0,0,0) )
		timelabel.alignment = QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignBottom
		vp.overlays.append(timelabel)
		
		filename = os.path.join(dir_imgs, filename_edited+'_{}.png'.format( str(i).zfill(4)) )
		print( filename )
		i += 1
		vp.render_image(size=(800,800), filename=filename, background=(1,1,1),frame=target_frame, renderer=OpenGLRenderer())
		vp.overlays.remove(timelabel)
	return
	
	
	
if __name__ == '__main__':
	
	
	frame_num_after_photobleach  = 300
	frame_num_before_photobleach = 100
	num_skip_frames_for_sampling = 1
	
	
	# Target: SPG
	filename_lammpstrj = 'PIPS_activated_trj.lammpstrj'
	filename_video     = 'PIPS'
	dir_lammpstrj_sub  = 'SPG'
	dir_target         = 'special_conditions'
	
	# Shared init
	dir_lammpstrj    = os.path.join('..','lammpstrj4', dir_lammpstrj_sub)
	dir_imgs         = os.path.join('imgs4', dir_target, 'for_movie')
	dir_videos       = os.path.join('imgs4', dir_target, 'movie')
	
	os.makedirs(dir_imgs, exist_ok=True)
	os.makedirs(dir_videos, exist_ok=True)
	
	# Load lammpstrj file.
	data_all = import_file(os.path.join(dir_lammpstrj, filename_lammpstrj), input_format= "lammps/dump" )
	
	# Make images
	plot_snapshots(data_all, dir_imgs, filename_video, \
		frame_num_after_photobleach,  
		frame_num_before_photobleach, 
		num_skip_frames_for_sampling
		)
	
	# Save a video
	ffname = os.path.join(dir_imgs, '{}_%04d.png'.format(filename_video))
	targ_name = os.path.join(dir_videos,'{}.mp4'.format(filename_video) )
	com = ['ffmpeg','-r', '5', \
		'-i', ffname,\
		'-vf','scale=trunc(iw/2)*2:trunc(ih/2)*2',\
		'-pix_fmt', 'yuv420p', targ_name, '-y']
	print(' '.join(com))
	s.call(com)
	
	#if os.path.isdir(dir_imgs):
	#	shutil.rmtree(dir_imgs)
