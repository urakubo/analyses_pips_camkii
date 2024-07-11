

import os, sys, glob, pprint,itertools, shutil
import numpy as np
os.environ['OVITO_GUI_MODE'] = '1' # Request a session with OpenGL support

import subprocess as s


from ovito.io import import_file
from ovito.modifiers import SliceModifier, SelectTypeModifier, AssignColorModifier
from ovito.vis import Viewport, TextLabelOverlay, OpenGLRenderer
from ovito.pipeline import Pipeline, StaticSource
from ovito.qt_compat import QtCore
from ovito.data import DataCollection


import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph
import lib.utils_ovito as utils_ovito

from specification_datasets import SpecDatasets


rm = 15
def get_molecular_ids_in_unbleached_area(position):
	x, y, z               = position[:,0], position[:,1], position[:,2]
	xc, yc, zc, rm_       = 0, rm, 0, rm
	ids_unbleached_moleules = ( (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc) > rm_ * rm_ )
	return ids_unbleached_moleules
	
	
	
def plot_snapshots(data_all, dir_edited_data, dir_imgs, \
		frame_num_after_photobleach  = 300, \
		frame_num_before_photobleach = 10, \
		num_skip_frames_for_sampling = 1, \
		target_molecule = 'Both'
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
	
	
	'''
	# Centering at the specified frame
	frame_time_zero = frame_photobleach
	types, positions, ids_molecule = utils.decode_data(data_all.compute( frame_time_zero ))
	center = utils.get_center_of_mass(types, positions)
	for dim in [0,1,2]:
		center[dim] += - p.space[dim] * (center[dim]  >=  p.space[dim]) + p.space[dim] * (center[dim]  <  0)
	modifier = utils_ovito.CenteringModifier()
	modifier.center = center
	data_all.modifiers.append(modifier)
	'''
	
	# View
	vp = Viewport()
	vp.type = Viewport.Type.Perspective
	vp.fov = np.radians(10.0) # math.radians
	vp.camera_pos = (850, 60, 60)
	vp.camera_dir = (-1, 0, 0)
	
	data_all.add_to_scene()	
	
	i        = 0
	f_bleach = 0
	for target_frame, target_time in zip(target_frames, time_steps):
		print('Time frame: ', target_frame, ', Time: ', target_time)
		
		if target_molecule == 'Both':
			ids_bead_cluster = utils_ovito.get_beads_in_largest_cluster(data_all, target_frame)
		elif target_molecule == 'CaMKII':
			ids_bead_cluster = utils_ovito.get_CaMKII_beads_in_largest_cluster(data_all, target_frame)
		elif target_molecule == 'GluN2B':
			ids_bead_cluster = utils_ovito.get_GluN2B_beads_in_largest_cluster(data_all, target_frame)
		elif target_molecule == 'Both_condensate_diluent':
			data   = data_all.compute(target_frame)
			ids_bead_cluster = np.arange(data.particles.count)
		else:
			sys.exit('target_molecule was inapproriate.')
		
		
		# Slice
		data      = data_all.compute(target_frame)
		positions = np.array( data.particles['Position'] )
		sliced    = set( np.flatnonzero( (positions[:,0] < 60+15) * (positions[:,0] > 60-15) ))
		ids_bead_cluster = list( set(ids_bead_cluster) & sliced )
		
		# Photobleach
		if (target_frame >= frame_photobleach) and (f_bleach == 0):
			#data   = data_all.compute(target_frame)
			#types, positions, ids_molecule = utils.decode_data(data)
			
			'''
			center = (60,60-15,60)
			'''
			
			types, positions, ids_molecule = utils.decode_data(data_all.compute( target_frame ))
			center    = utils.get_center_of_mass(types, positions)
			
			position_centered       = utils.centering(positions, center)
			unbleached_moleules     = get_molecular_ids_in_unbleached_area(position_centered)
			ids_unbleached_moleules = set( np.flatnonzero(unbleached_moleules) )
			
			pipeline = utils_ovito.add_sphere_yellow(center, rm)
			pipeline.add_to_scene()
			print('Photo bleach')
			f_bleach += 1
			#print('ids_unbleached_moleules ', ids_unbleached_moleules)
		
		# Remove the bleached molecules
		if target_frame >= frame_photobleach:
			ids_bead_cluster = list( set(ids_bead_cluster) & ids_unbleached_moleules & sliced )
			#print('ids_bead_cluster  ',ids_bead_cluster )
		
		
		
		# Exec remove
		modifier_del = utils_ovito.DeleteParticle()
		modifier_del.ids = ids_bead_cluster
		data_all.modifiers.append(modifier_del)
		
		
		# Set time
		'''
		timelabel = TextLabelOverlay( text = 'Time: {:.4f} G MC step'.format(target_time),\
			offset_x = 0.08 ,\
			offset_y = 0.005, \
			font_size = 0.03, \
			text_color = (0,0,0) )
		timelabel.alignment = QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignBottom
		vp.overlays.append(timelabel)
		filename = os.path.join(dir_imgs, '{}.png'.format( str(i).zfill(4)) )
		'''
		
		filename = os.path.join(dir_imgs, '{:+.4f}_G_MC_step.png'.format(target_time) )
		print( filename )
		i += 1
		vp.render_image(size=(800,800), filename=filename, background=(1,1,1),frame=target_frame, renderer=OpenGLRenderer())
		#vp.overlays.remove(timelabel)
		del data_all.modifiers[-1]
		
	
	return
	
	
class MakeOvitoVideoFRAP(SpecDatasets):
	
	def __init__( self ):
		
		self.num_frames_after_photobleach = 300
		self.num_frames_before_photobleach = 10
		self.num_skip_frames_for_sampling = 1
		
		
	def inspect( self ):
		for i, (f_lammpstrj, f_edited) in enumerate( zip(self.filenames_lammpstrj, self.filenames_edited) ):
			print('ID {}: {}, {}'.format(i, f_lammpstrj, f_edited))
		
		
	def run( self, i, target_molecule = 'CaMKII' ):
		
		print('\nID {}: {}, {}'.format(i, self.filenames_lammpstrj[i], self.filenames_edited[i]))
		# Shared init
		self.dir_imgs         = os.path.join( self.dir_imgs_root, 'for_movie_{}'.format( self.filenames_edited[i] ) )
		os.makedirs(self.dir_imgs, exist_ok=True)
		
		# Load lammpstrj file.
		data_all = import_file(os.path.join(self.dir_lammpstrj, self.filenames_lammpstrj[i]), input_format= "lammps/dump" )
		
		self.num_frames_before_photobleach, self.num_frames_after_photobleach = self.set_frames_before_after[i]
		
		# Make images
		plot_snapshots(data_all, self.dir_edited_data, self.dir_imgs, \
			self.num_frames_after_photobleach,  
			self.num_frames_before_photobleach, 
			self.num_skip_frames_for_sampling, 
			target_molecule = target_molecule
			)
		
		
	def make_a_video( self, i ):
		
		dir_videos    = os.path.join( self.dir_imgs_root, 'movies_FRAP' )
		self.dir_imgs = os.path.join( self.dir_imgs_root, 'for_movie_{}'.format( self.filenames_edited[i] ) )
		os.makedirs(dir_videos, exist_ok=True)
		
		ffname = os.path.join( self.dir_imgs, '%04d.png' )
		targ_name = os.path.join(dir_videos,'{}.mp4'.format( self.filenames_edited[i] ) )
		com = ['ffmpeg','-r', '5', \
			'-i', ffname,\
			'-vf','scale=trunc(iw/2)*2:trunc(ih/2)*2',\
			'-pix_fmt', 'yuv420p', targ_name, '-y']
		print(' '.join(com))
		s.call(com)
	
	
	
if __name__ == '__main__':
	
	
	# Small colony 2
	obj = MakeOvitoVideoFRAP()
	obj.valency_length_small_colony2()
	obj.inspect()
	i = 7*5+2 # val_12\R2_002
	target_molecule = 'CaMKII' # 'CaMKII', 'GluN2B', 'Both'
	
	obj.run(i, target_molecule)
	obj.make_a_video(i)
	
	
	'''
	# Small colony 3
	dir_target  = 'small_colony3'
	i = 7*5+6 # val_12\R2_006
	i = 7*5+2 # val_12\R2_006
	
	'''
	
	
	
	
#if os.path.isdir(dir_imgs):
#	shutil.rmtree(dir_imgs)

