#
# I referred the following web page.
# https://www.ovito.org/docs/current/python/introduction/pipelines.html
# 

import os, glob
os.environ['OVITO_GUI_MODE'] = '1' # Request a session with OpenGL support

import numpy as np

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
from lib.specification_datasets import SpecDatasets


def plot_snapshots(data_all, time_frame, dir_imgs, filename_output):
	
	# Fix bounging box
	data_all.modifiers.append( utils_ovito.FixBoundingbox() )
	# Non boundary
	cell_vis = data_all.source.data.cell.vis
	#cell_vis.render_cell = False
	cell_vis.line_width = 0.5
	##cell_vis.rendering_color = (1.0, 1.0, 1.0)
	
	
	# Show particles only in the largest cluster.
	'''
	modifier_del     = utils_ovito.DeleteParticle()
	#modifier_del.ids = utils_ovito.get_beads_in_largest_cluster(data_all, time_frame)
	modifier_del.ids = utils_ovito.get_CaMKII_beads_in_largest_cluster(data_all, time_frame)
	data_all.modifiers.append(modifier_del)
	'''
	
	
	# Centering
	'''
	num_frames = data_all.source.num_frames
	types, positions, ids_molecule = utils.decode_data(data_all.compute(num_frames))
	center = utils.get_center_of_mass(types, positions)
	for dim in [0,1,2]:
		center[dim] += - p.space[dim] * (center[dim]  >=  p.space[dim]) + p.space[dim] * (center[dim]  <  0)
	
	modifier = utils_ovito.CenteringModifier()
	modifier.center = center
	data_all.modifiers.append(modifier)	
	'''
	
	
	# Slice data
	#'''
	modifier = SliceModifier()
	modifier.normal   = (1.0, 0.0, 0.0)
	modifier.distance = 60
	data_all.modifiers.append(modifier)
	#'''
	
	
	for k, v in p.molecules_without_all.items():
		data_all.modifiers.append(SelectTypeModifier(types=set(v['id'])))
		data_all.modifiers.append(AssignColorModifier(color=c.cmap_universal_ratio[k] ))
	
	
	vp = Viewport()
	vp.type = Viewport.Type.Perspective
	vp.fov = np.radians(10.0)
	vp.camera_pos = (850, 60, 60)
	vp.camera_dir = (-1, 0, 0)
	'''
	dist = 200/np.sqrt(3)+10
	vp.camera_pos = (60+dist, 60+dist, 60+dist)
	vp.camera_dir = (-1, -1, -1)
	'''
	data_all.add_to_scene()
	filename = os.path.join(dir_imgs, filename_output+'_{}.png'.format( str(time_frame).zfill(4)) )
	vp.render_image(size=(800,800), filename=filename, background=(1,1,1),frame=time_frame, renderer=OpenGLRenderer())
	
	return


	# Scale bar
	# https://ovito.org/manual/python/introduction/examples/overlays/scale_bar.html
	
		
class Plot3dOvito(SpecDatasets):
	def __init__( self ):
		
		pass
		
		
	def inspect( self ):
		for i, (f_lammpstrj, f_edited) in enumerate( zip(self.filenames_lammpstrj, self.filenames_edited) ):
			print('ID {}: {}, {}'.format(i, f_lammpstrj, f_edited))
		
		
	def run( self, i ):
		
		##
		dir_lammpstrj   = os.path.join('..'   , 'lammpstrj4', self.dir_target)
		dir_imgs        = os.path.join('imgs4', self.dir_target, 'profiles_3d_ovito' )
		os.makedirs(dir_imgs, exist_ok=True)
		##
		
		
		sampling_frame = utils.get_num_frames(dir_lammpstrj, self.filenames_lammpstrj[i])
		data_all   = import_file(os.path.join(dir_lammpstrj, self.filenames_lammpstrj[i]), input_format= "lammps/dump" )
		plot_snapshots(data_all, sampling_frame, dir_imgs, self.filenames_edited[i])
		print('filename_input: ', self.filenames_lammpstrj[i])
		print('sampling_frame: ', sampling_frame )
		
		
if __name__ == '__main__':
	
	# I manually ran each one of them,
	# because I do not know how to fully reset the ovito visualization system.
	
	obj = Plot3dOvito()
	obj.boundary_conditions2() #  conc_dependence(), valency_length(), valency_length_CG(), boundary_conditions2()
	# obj.inspect()
	obj.run(1)
	
	
