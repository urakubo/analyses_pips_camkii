

import os, sys, glob, pprint,itertools, shutil
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
from ovito.data import *
from ovito.pipeline import *

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph



rm = 15
def get_molecular_ids_in_unbleached_area(position):
	x, y, z               = position[:,0], position[:,1], position[:,2]
	xc, yc, zc, rm_       = 0, rm, 0, rm
	ids_unbleached_moleules = ( (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc) > rm_ * rm_ )
	return ids_unbleached_moleules


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
	
def add_sphere_yellow(center):
	data_photobleach = DataCollection()
	particles = data_photobleach.create_particles()
	pos = [(center[0]+120, center[1]+120+rm, center[2]+120)]
	pos_prop  = particles.create_property('Position', data=pos)
	type_prop = particles.create_property('Particle Type')
	type_prop.types.append(ParticleType(id = 1, name = 'Photobleach', color = (1.0,1.0,0.0), radius=rm))
	type_prop[0] = 1
	particles.create_property('Transparency', data=0.6)
	pipeline = Pipeline(source = StaticSource(data = data_photobleach))
	return pipeline
	
	
def plot_snapshots(data_all, dir_edited_data, dir_imgs, filename_edited):
	
	
	# Create the data collection containing a Particles object:

	# Time frames
	num_frames = data_all.source.num_frames
	frame_num_after_photobleach  = 300
	frame_num_before_photobleach = 10
	num_skip_frames_for_sampling = 1
	
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
	
	
	vp = Viewport()
	vp.type = Viewport.Type.Perspective
	#vp.fov = math.radians(40.0)
	#vp.camera_pos = (250+60, 60, 60)
	vp.fov = np.radians(10.0) # math.radians
	vp.camera_pos = (850, 60, 60)
	vp.camera_dir = (-1, 0, 0)
	
	i        = 0
	f_bleach = 0
	for target_frame, target_time in zip(target_frames, time_steps):
		print('Time frame: ', target_frame, ', Time: ', target_time)
		# decode data
		data   = data_all.compute(target_frame)
		types, positions, ids_molecule = utils.decode_data(data)
		bp 	= np.array( data.particles['bp'] ).astype('int')
		
		
		# Get ids in the largest cluster
		multi_graph = utils_graph.get_multi_graph(ids_molecule, types, bp, positions)
		clusters = sorted(nx.connected_components(multi_graph), key=len, reverse=True)
		g_largest_cluster = nx.MultiGraph( multi_graph.subgraph(clusters[0]) )
		
		ids_GluN2B = [v['ids_bead'][0].tolist() for n, v in g_largest_cluster.nodes.items() if v['species'] == 'GluN2B' ]
		ids_CaMKII = [v['ids_bead'][0].tolist() for n, v in g_largest_cluster.nodes.items() if v['species'] == 'CaMKII' ]
		ids_GluN2B = utils.flatten(ids_GluN2B)
		ids_CaMKII = utils.flatten(ids_CaMKII)
		#ids_bead_cluster = ids_GluN2B + ids_CaMKII
		#ids_bead_cluster = ids_CaMKII
		ids_bead_cluster = ids_GluN2B
		
		# Photobleach
		
		# Photobleach
		if (target_frame >= frame_photobleach) and (f_bleach == 0):
			center    = utils.get_center_of_mass(types, positions)
			position_centered = utils.centering(positions, center)
			unbleached_moleules     = get_molecular_ids_in_unbleached_area(position_centered)
			ids_unbleached_moleules = set( np.flatnonzero(unbleached_moleules) )
			
			pipeline = add_sphere_yellow(center)
			pipeline.add_to_scene()
			
			print('Photo bleach')
			f_bleach += 1
			#print('ids_unbleached_moleules ', ids_unbleached_moleules)
		
		# Remove the bleached molecules
		if target_frame >= frame_photobleach:
			ids_bead_cluster = list( set(ids_bead_cluster) & ids_unbleached_moleules )
			#print('ids_bead_cluster  ',ids_bead_cluster )
		
		#print('ids_bead ', ids_bead)
		
		modifier_del = DeleteParticle()
		modifier_del.ids = ids_bead_cluster
		data_all.modifiers.append(modifier_del)
		
		data   = data_all.compute()
		data_all.add_to_scene()	
		
		
		# Set time
		timelabel = TextLabelOverlay( text = 'Time: {:.5f} G MC step'.format(target_time),\
			font_size = 0.03, \
			text_color = (0,0,0) )
		timelabel.alignment = QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignBottom
		vp.overlays.append(timelabel)
		
		filename = os.path.join(dir_imgs, filename_edited+'_{}.png'.format( str(i).zfill(4)) )
		print( filename )
		i += 1
		vp.render_image(size=(800,800), filename=filename, background=(1,1,1),frame=target_frame, renderer=OpenGLRenderer())
		vp.overlays.remove(timelabel)
		del data_all.modifiers[-1]
	
	return
	
	
	
if __name__ == '__main__':
	
	
	# Input file
	
	# Small colony 2
	subdirs    = ['val_{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
	filenames_lammpstrj = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_edited    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7)]
	dir_target  = 'small_colony2'
	
	max_backward_frames_for_sampling = 80
	num_skip_frames_for_sampling     = 1
	
	i = 7*5+2 # val_12\R2_002
	#i = 7*5+6 # val_12\R2_006
	
	filename_lammpstrj = filenames_lammpstrj[i]
	filename_edited    = filenames_edited[i]
	
	
	
	## Shared init
	#fig_title = '{}_{}'.format(nth_largest, prefix)
	#print(fig_title)
	
	dir_lammpstrj    = os.path.join('..','lammpstrj4', dir_target)
	dir_edited_data  = os.path.join('data4', dir_target)
	dir_imgs         = os.path.join('imgs4', dir_target, 'FRAP_imgs_for_movie',filename_edited)
	dir_videos       = os.path.join('imgs4', dir_target, 'FRAP_movie')
	
	#if os.path.isdir(dir_imgs):
	#	shutil.rmtree(dir_imgs)
	os.makedirs(dir_imgs, exist_ok=True)
	os.makedirs(dir_videos, exist_ok=True)
	
	# Load connect_graph.
	
	data_all = import_file(os.path.join(dir_lammpstrj, filename_lammpstrj), input_format= "lammps/dump" )
	plot_snapshots(data_all, dir_edited_data, dir_imgs, filename_edited)
	
	
	ffname = os.path.join(dir_imgs, '{}_%04d.png'.format(filename_edited))
	targ_name = os.path.join(dir_videos,'{}.mp4'.format(filename_edited) )
	com = ['ffmpeg','-r', '5', \
		'-i', ffname,\
		'-vf','scale=trunc(iw/2)*2:trunc(ih/2)*2',\
		'-pix_fmt', 'yuv420p', targ_name, '-y']
	print(' '.join(com))
	s.call(com)
	
	
