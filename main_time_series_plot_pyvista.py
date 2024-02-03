#
# I referred the following web page.
# https://www.ovito.org/docs/current/python/introduction/pipelines.html
# 

import os, glob
os.environ['OVITO_GUI_MODE'] = '1' # Request a session with OpenGL support

import numpy as np
import matplotlib.pyplot as plt
from ovito.io import import_file
from ovito.modifiers import SliceModifier, SelectTypeModifier, ColorByTypeModifier, AssignColorModifier, CoordinationPolyhedraModifier, ComputePropertyModifier
import math
from ovito.vis import *
import utils
import colormap as c


def compute_particle_transparency(frame, data):
    transparency = np.ones((data.particles.count))*0.5
    data.particles_.create_property('Transparency', data=transparency)
  	
	
def plot_snapshots(data_all, time_frames, dir_imgs, filename_output):
	
	data   = data_all.compute()
	
	modifier = SliceModifier()
	modifier.normal   = (1.0, 0.0, 0.0)
	modifier.distance = 60
	data_all.modifiers.append(modifier)
	data_all.add_to_scene()

	for k, v in utils.molecules_without_all.items():
		data_all.modifiers.append(SelectTypeModifier(types=set(v['id'])))
		data_all.modifiers.append(AssignColorModifier(color=c.cmap_universal_ratio[k] ))
		#data_all.modifiers.append(compute_particle_transparency)


	vp = Viewport()
	vp.type = Viewport.Type.Perspective
	vp.fov = math.radians(60.0)

	vp.camera_pos = (170+60, 60, 60)
	vp.camera_dir = (-1, 0, 0)
	
	for time_frame in time_frames:
		filename = os.path.join(dir_imgs, filename_output+'_{}.png'.format( str(time_frame).zfill(4)) )
		vp.render_image(size=(1200,1200), filename=filename, background=(1,1,1), frame=time_frame, renderer=OpenGLRenderer())

	'''
	dist = 200/np.sqrt(3)+10
	vp.camera_pos = (60+dist, 60+dist, 60+dist)
	vp.camera_dir = (-1, -1, -1)
	'''

	return data_all


if __name__ == '__main__':


	dir_input       = './../lammpstrj/mix_two_three_components'
	filenames_input = [		'CG.lammpstrj',\
							'SP.lammpstrj',\
							'SPG1.lammpstrj',\
							'SPG2.lammpstrj']
	filenames_output = [	'CG',\
							'SP',\
							'SPG1',\
							'SPG2']
	
	time_frames = list(range(0,200,20))
	print('Recording time frame: ', time_frames)
	
	for filename_input, filename_output in zip(filenames_input, filenames_output):
		print('filename_input: ', filename_input)
		data_all   = import_file(os.path.join(dir_input, filename_input), input_format= "lammps/dump" )
		dir_imgs   = os.path.join( 'imgs', 'time_series', 'snapshot', filename_output )
		os.makedirs(dir_imgs, exist_ok=True)
	
		plot_snapshots(data_all, time_frames, dir_imgs, filename_output)

