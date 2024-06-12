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


def compute_particle_transparency(frame, data):
    transparency = np.ones((data.particles.count))*0.3
    data.particles_.create_property('Transparency', data=transparency)
  	
	
def plot_snapshots(data_all, time_frames, dir_imgs, filename_output):
	
	data   = data_all.compute()

	# Bounding box
	# https://ovito.org/manual/python/introduction/examples/modifiers/shrink_wrap_box.html
	#   (x_max-x_min  0            0            x_min)
	#   (0            y_max-y_min  0            y_min)
	#   (0            0            z_max-z_min  z_min)
	coords_min = np.array([0,0,0])
	matrix = np.empty((3,4))
	matrix[:,:3] = np.diag( utils.space_np - coords_min )
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

	for k, v in utils.molecules_without_all.items():
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
	
	for time_frame in time_frames:
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

	## Dataset 1
	dir_input       = './../lammpstrj/mix_two_three_components'
	filenames_input = [		'CG.lammpstrj',\
							'SP.lammpstrj',\
							'SPG1.lammpstrj',\
							'SPG2.lammpstrj']
	filenames_output = [	'CG',\
							'SP',\
							'SPG1',\
							'SPG2']
	sampling_interval = 20
	sampling_time_frames = list(range(0,200+sampling_interval,sampling_interval))
	

	
	i = 2 # 0-3 , I do not know how to completely reset the ovito visualization.

	## Dataset 2
	dir_input       = './../lammpstrj/self_affinity'
	filenames_input = [		'OnlyCaMKIICaMKIIAffinity.lammpstrj',\
							'OnlyGluN2BGluN2B.lammpstrj',\
							'onlySTGSTG.lammpstrj']
	filenames_output = [	'CaMKIIalone',\
							'GluN2Balone',\
							'STGalone']	
	
	sampling_time_framess = [ [0, 100, 200], [0, 100, 187] , [0, 100, 195] ]
	sampling_time_frames = sampling_time_framess[i]
	###
	
	
	
	print('Recording time frame: ', sampling_time_frames)
	
	# for filename_input, filename_output in zip(filenames_input, filenames_output):
	
	print('filename_input: ', filenames_input[i])
	data_all   = import_file(os.path.join(dir_input, filenames_input[i]), input_format= "lammps/dump" )
	dir_imgs   = os.path.join( 'imgs', 'time_series', 'snapshot', filenames_output[i] )
	os.makedirs(dir_imgs, exist_ok=True)

	plot_snapshots(data_all, sampling_time_frames, dir_imgs, filenames_output[i])

