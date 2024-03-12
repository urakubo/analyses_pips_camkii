#
# I referred the following web page.
# https://www.ovito.org/docs/current/python/introduction/pipelines.html
# 

import os, glob
os.environ['OVITO_GUI_MODE'] = '1' # Request a session with OpenGL support

import numpy as np
import matplotlib.pyplot as plt
from ovito.io import import_file
from ovito.modifiers import SliceModifier, SelectTypeModifier, ColorByTypeModifier, AssignColorModifier, CoordinationPolyhedraModifier, ComputePropertyModifier, GenerateTrajectoryLinesModifier
import math
from ovito.vis import *
import ovito.data
from ovito.pipeline import *

import utils
import colormap as c
import parameters as p


def compute_particle_transparency(frame, data):
    transparency = np.ones((data.particles.count))*0.3
    data.particles_.create_property('Transparency', data=transparency)
  	
	
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
		#data_all.modifiers.append(compute_particle_transparency)

	# No boundary
	cell_vis = data_all.source.data.cell.vis
	#cell_vis.render_cell = False
	cell_vis.line_width = 0.5
	##cell_vis.rendering_color = (1.0, 1.0, 1.0)

	'''
	line_vis = data_all.source.data.line.vis
	line_vis.color = (0.0, 0.0, 0.0)
	line_vis.width = 0.5
	positions = np.array([[0,0,60],[0,120,60],[120,120,60],[120,0,60],[0,0,60] ])
	lines_data = ovito.data.Lines
	lines_data.create_line = positions
	lines_data.color = (0.0, 0.0, 0.0)
	lines_data.width = 0.5

	line = GenerateTrajectoryLinesModifier(only_selected=False)
	# line.generate()
	line.vis.color = (0.0, 0.0, 0.0)
	line.vis.width = 0.5
	data_all.modifiers.append(line)
	'''
	
	
	vp = Viewport()
	vp.type = Viewport.Type.Perspective
	#vp.fov = math.radians(40.0)
	#vp.camera_pos = (250+60, 60, 60)
	vp.fov = math.radians(10.0)
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
	
	# I manually ran each of them,
	# because I do not know how to reset the ovito visualization completely .
	i = 2 # 0-3 for dataset 1, 0-2 for dataset 2, 0-2 for dataset 3
	
	
	

	## Dataset 1
	filenames_lammpstrj = ['CG_trj.lammpstrj', 'CPG_trj.lammpstrj', 'SP_trj.lammpstrj' ]
	dir_lammpstrj    = 'Feb_Sub_Figure1'
	filenames_output = ['CG','CPG','SP']
	dir_output = 'types_mixture'
	
	##
	dir_lammpstrj   = os.path.join('..'   , 'lammpstrj2',dir_lammpstrj)
	dir_imgs        = os.path.join('imgs2', dir_output, 'snapshot' )
	os.makedirs(dir_imgs, exist_ok=True)
	##
	
	
	
	print('filename_input: ', filenames_lammpstrj[i])
	sampling_frame = utils.get_num_frames(dir_lammpstrj, filenames_lammpstrj[i])
	print("sampling_frame ", sampling_frame )
	
	data_all   = import_file(os.path.join(dir_lammpstrj, filenames_lammpstrj[i]), input_format= "lammps/dump" )
	plot_snapshots(data_all, sampling_frame, dir_imgs, filenames_output[i])

