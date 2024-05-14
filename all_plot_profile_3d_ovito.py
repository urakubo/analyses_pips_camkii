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

  	
	
def plot_snapshots(data_all, time_frame, dir_imgs, filename_output):
	
	data   = data_all.compute()

	# Bounding box
	# https://ovito.org/manual/python/introduction/examples/modifiers/shrink_wrap_box.html
	#   (x_max-x_min  0            0            x_min)
	#   (0            y_max-y_min  0            y_min)
	#   (0            0            z_max-z_min  z_min)
	coords_min = np.array([0,0,0])
	matrix = np.empty((3,4))
	matrix[:,:3] = np.diag( p.space_np - coords_min )
	matrix[:, 3] = coords_min
	# Assign the cell matrix - or create whole new SimulationCell object in
	# the DataCollection if there isn't one already.
	data.create_cell(matrix, (False, False, False))

	# Slice data
	modifier = SliceModifier()
	modifier.normal   = (1.0, 0.0, 0.0)
	modifier.distance = 60
	data_all.modifiers.append(modifier)
	data_all.add_to_scene()

	for k, v in p.molecules_without_all.items():
		data_all.modifiers.append(SelectTypeModifier(types=set(v['id'])))
		data_all.modifiers.append(AssignColorModifier(color=c.cmap_universal_ratio[k] ))

	# No boundary
	cell_vis = data_all.source.data.cell.vis
	#cell_vis.render_cell = False
	cell_vis.line_width = 0.5
	##cell_vis.rendering_color = (1.0, 1.0, 1.0)
	
	
	vp = Viewport()
	vp.type = Viewport.Type.Perspective
	vp.fov = np.radians(10.0)
	vp.camera_pos = (850, 60, 60)
	vp.camera_dir = (-1, 0, 0)
	
	filename = os.path.join(dir_imgs, filename_output+'_{}.png'.format( str(time_frame).zfill(4)) )
	vp.render_image(size=(800,800), filename=filename, background=(1,1,1),frame=time_frame, renderer=OpenGLRenderer())
	
	'''
	dist = 200/np.sqrt(3)+10
	vp.camera_pos = (60+dist, 60+dist, 60+dist)
	vp.camera_dir = (-1, -1, -1)
	'''
	
	return data_all


	# Scale bar
	# https://ovito.org/manual/python/introduction/examples/overlays/scale_bar.html

if __name__ == '__main__':
	
	'''
	# I manually ran each one of them,
	# because I do not know how to fully reset the ovito visualization system.
	i = 6
	# Special conditions
	filenames_output 	= ['CPG', 'SP', 'SPG', 'CG_000','CG_001','CG_002','CG_003']
	filenames_lammpstrj = [os.path.join('CPG','R2_000.lammpstrj'), \
						os.path.join('SP','R2_000.lammpstrj'), \
						os.path.join('SPG','R2_000.lammpstrj'), \
						os.path.join('binary_CG','R2_000.lammpstrj'), \
						os.path.join('binary_CG','R2_001.lammpstrj'), \
						os.path.join('binary_CG','R2_002.lammpstrj'), \
						os.path.join('binary_CG','R2_003.lammpstrj') \
						]
	dir_target     = 'special_conditions'
	dir_lammpstrj  = '.'
	'''
	
	
	## Conc dependence
	filenames_output    = [str(i).zfill(3) for i in range(81) ]
	filenames_lammpstrj = ['R2_{}.lammpstrj'.format(f) for f in filenames_output ]
	dir_lammpstrj    = 'conc_dependence'
	dir_target       = 'conc_dependence'
	i = 9*6+4
	
	
	##
	dir_lammpstrj   = os.path.join('..'   , 'lammpstrj4',dir_lammpstrj)
	dir_imgs        = os.path.join('imgs4', dir_target, 'profiles_3d_ovito' )
	os.makedirs(dir_imgs, exist_ok=True)
	##
	
	
	
	print('filename_input: ', filenames_lammpstrj[i])
	sampling_frame = utils.get_num_frames(dir_lammpstrj, filenames_lammpstrj[i])
	print("sampling_frame ", sampling_frame )
	
	data_all   = import_file(os.path.join(dir_lammpstrj, filenames_lammpstrj[i]), input_format= "lammps/dump" )
	plot_snapshots(data_all, sampling_frame, dir_imgs, filenames_output[i])

