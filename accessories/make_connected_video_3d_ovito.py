

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
import lib.utils_ovito as utils_ovito
from specification_datasets import SpecDatasets



def plot_snapshots(data_all, dir_imgs, \
				center, \
				i_frame, \
				accum_mc_steps, \
				num_skip_frames_for_sampling):
	
	# Get sampling frames and times
	num_frames      = data_all.source.num_frames
	target_frames   = list( range(0, num_frames, num_skip_frames_for_sampling) )
	target_times    = np.array([data_all.compute(t).attributes['Timestep'] / 1e9 for t in target_frames]) + accum_mc_steps / 1e9
	
	print( 'target_times' )
	print( target_times )
	
	# label colors
	for k, v in p.molecules_without_all.items():
		data_all.modifiers.append(SelectTypeModifier(types=set(v['id'])))
		data_all.modifiers.append(AssignColorModifier(color=c.cmap_universal_ratio[k] ))
	
	
	# Centering
	'''
	modifier = utils_ovito.CenteringModifier()
	modifier.center = center
	data_all.modifiers.append(modifier)
	'''
	modifier = utils_ovito.CenteringEachFrameModifier()
	data_all.modifiers.append(modifier)
	
	
	# Slice data
	modifier = SliceModifier()
	modifier.normal   = (1.0, 0.0, 0.0)	# (0.0, 1.0, 0.0)
	modifier.distance = 60
	data_all.modifiers.append(modifier)
	
	vp = Viewport()
	vp.type = Viewport.Type.Perspective
	vp.fov = np.radians(10.0)
	vp.camera_pos = (850, 60, 60)		# (60, 850, 60)
	vp.camera_dir = (-1, 0, 0)			# (0, -1, 0)

	
	data_all.add_to_scene()
	
	for target_frame, target_time in zip(target_frames, target_times):
		print('Time frame: {}, Time: {}, Save ID: {}'.format(target_frame, target_time, i_frame) )
		
		# Set time
		timelabel = TextLabelOverlay(
			text = 'Time after the activation of CaMKII: {:.1f} G MC steps'.format(target_time),\
			font_size = 0.025 ,\
			text_color = (0,0,0) ,\
			alignment = QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignBottom ,\
			offset_x = 0.08 ,\
			offset_y = 0.005
			)
		vp.overlays.append(timelabel)
		
		filename = os.path.join(dir_imgs, '{}.png'.format( str(i_frame).zfill(4)) )
		print( filename )
		vp.render_image(size=(800,800), filename=filename, background=(1,1,1),frame=target_frame, renderer=OpenGLRenderer())
		vp.overlays.remove(timelabel)
		i_frame += 1
	
	data_all.remove_from_scene()
	return i_frame
	
	
	
class MakeOvitoConnectedVideo(SpecDatasets):
		
	def __init__( self ):
		
		self.frame_time_zero = 100
		self.num_skip_frames_for_sampling = 5
		
		
	def inspect( self ):
		
		for i, (f_lammpstrj, f_edited) in enumerate( zip(self.filenames_lammpstrj, self.filenames_edited) ):
			print('ID {}: {}, {}'.format(i, f_lammpstrj, f_edited))
		print('')
		
		
	def make_connected_snapshots( self, ids_target ):
		
		if not isinstance(ids_target, (list, tuple, range, set, frozenset)):
			print('Use list, tuple, range, set, or frozenset.')
			return
		
		self.folder_name = 'connected'
		dir_imgs         = os.path.join( self.dir_imgs_root, 'for_movie_{}'.format( self.folder_name ) )
		os.makedirs(dir_imgs, exist_ok=True)
		
		accum_frames    = 0
		accum_mc_steps  = 0
		accums_frames   = []
		accums_mc_steps = []
		
		# Calculate accumulated time for connection.
		print('Frame id for time zero and for centering: ', self.frame_time_zero)
		done_centering = False
		
		for i in ids_target:
			num_frames, mc_steps = utils.get_num_frames_mc_steps(self.dir_lammpstrj, self.filenames_lammpstrj[i])
			accums_frames.append( accum_frames )
			accums_mc_steps.append( accum_mc_steps )
			accum_frames   += num_frames
			accum_mc_steps += mc_steps
			
			print('ID {}, {}: frames {}, mc_steps {}, accum frames {}, accum mc_steps {}'.format(\
				i, self.filenames_lammpstrj[i], \
				num_frames, mc_steps, accum_frames, accum_mc_steps ) )
			
			# Centering 
			if (done_centering == False) & ( accum_frames > self.frame_time_zero ):
				done_centering  = True
				frame_centering = self.frame_time_zero - accums_frames[-1]
				types, positions,ids_molecule, mc_step = \
					utils.load_lammpstrj( self.dir_lammpstrj, self.filenames_lammpstrj[i], frame_centering )
				accum_mc_steps_centering = mc_step + accums_mc_steps[-1]
				
				center = utils.get_center_of_mass(types, positions)
				for dim in [0,1,2]:
					center[dim] += - p.space[dim] * (center[dim]  >=  p.space[dim]) + p.space[dim] * (center[dim]  <  0)
				print('Centering:')
				print('ID {}, {}: frames {}, mc_steps {}, accum frames {}, accum mc_steps {}'.format(\
					i, self.filenames_lammpstrj[i], \
					frame_centering,  mc_step, self.frame_time_zero, accum_mc_steps_centering) )
		
		
		accums_mc_steps = [a - accum_mc_steps_centering for a in accums_mc_steps]
		
		# Repeat run.
		i_frame = 0
		for i, accum_mc_steps  in zip(ids_target, accums_mc_steps):
			data_all = import_file(os.path.join(self.dir_lammpstrj, self.filenames_lammpstrj[i]), input_format= "lammps/dump" )
			i_frame = plot_snapshots(data_all, dir_imgs, \
					center, \
					i_frame, \
					accum_mc_steps, \
					self.num_skip_frames_for_sampling)
		
		
	def make_a_video( self ):
		#
		dir_videos    = os.path.join( self.dir_imgs_root, 'movies' )
		dir_imgs = os.path.join( self.dir_imgs_root, 'for_movie_{}'.format( self.folder_name ) )
		os.makedirs(dir_videos, exist_ok=True)
		
		ffname = os.path.join( dir_imgs, '%04d.png' )
		targ_name = os.path.join(dir_videos,'{}.mp4'.format( self.folder_name ) )
		com = ['ffmpeg','-r', '5', \
			'-i', ffname,\
			'-vf','scale=trunc(iw/2)*2:trunc(ih/2)*2',\
			'-pix_fmt', 'yuv420p', targ_name, '-y']
		print(' '.join(com))
		s.call(com)
		
		
if __name__ == '__main__':
	
	'''
	obj =  MakeOvitoConnectedVideo()
	obj.CaMKII_blocked4() 
	obj.inspect()
	ids_target = [0,1,2]
	obj.make_connected_snapshots(ids_target)
	obj.make_a_video()
	'''
	
	obj =  MakeOvitoConnectedVideo()
	obj.boundary_conditions5()
	obj.inspect()
	ids_target = [11,12,13]
	obj.make_connected_snapshots(ids_target)
	obj.make_a_video()	
	
	
	## R2_000  Spatial constraint, 10000 steps  
	## R2_001  3.2T to 1.0 T , 10000 steps
	# R2_003 SPG steady state and CaMKII inactivated, 1.0 T 
	# R2_004 CaMKII activated and the temperature increased o 1.2 T
	# R2_005 1.2 T, Long simulation to get the PIPS, this is still running, 
	
	
	