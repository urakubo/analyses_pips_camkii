##
##
##

import os, glob, sys
os.environ['OVITO_GUI_MODE'] = '1' # Request a session with OpenGL support

import numpy as np
import matplotlib.pyplot as plt
from ovito.io import import_file
from ovito.modifiers import SliceModifier, SelectTypeModifier, ColorByTypeModifier, AssignColorModifier, CoordinationPolyhedraModifier, ComputePropertyModifier, GenerateTrajectoryLinesModifier, DeleteSelectedModifier
import math
from ovito.vis import *
import ovito.data
from ovito.pipeline import *

import networkx as nx

import utils
import colormap as c
import parameters as p

id_shared_PSD95 = 20

def plot_snapshots(data_all, time_frame, dir_imgs, filename_output, ids_beads_shared_PSD95):
	
	data   = data_all.compute()
	
	# Label shared PSD95
	for id in ids_bead_both:
		data.particles_.particle_types_[id] = id_shared_PSD95
	
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
	
	
	def insert_shared_PSD95(frame, data):
	    for i in ids_beads_shared_PSD95:
	    	data.particles_.particle_types_[i] = id_shared_PSD95
	    
	data_all.modifiers.append(insert_shared_PSD95) 
	
	molecule_ids = {k: v['id'] for k, v in p.molecules_without_all.items()}
	molecule_ids['Shared PSD95'] = [id_shared_PSD95]
	for k, v in molecule_ids.items():
		data_all.modifiers.append(SelectTypeModifier(types=set(v)))
		if k not in ['PSD95', 'Shared PSD95']:
			data_all.modifiers.append(DeleteSelectedModifier())
		elif k == 'PSD95':
			data_all.modifiers.append(AssignColorModifier(color=(191/255,228/255,255/255) ))
		elif k == 'Shared PSD95':
			data_all.modifiers.append(AssignColorModifier(color=c.cmap_universal_ratio[k] ))
	
	
	#data_all.modifiers.append(SelectTypeModifier(types=[id_shared_PSD95]))
	# data_all.modifiers.append(AssignColorModifier(color=c.cmap_universal_ratio['Shared PSD95']) )
	#data_all.modifiers.append(AssignColorModifier(color=c.cmap_universal_ratio['STG']) )

	#data_all.modifiers.append(SelectTypeModifier(types=[1]))
	#data_all.modifiers.append(DeleteSelectedModifier())


	# Slice data
	modifier = SliceModifier()
	modifier.normal   = (1.0, 0.0, 0.0)
	modifier.distance = 60
	data_all.modifiers.append(modifier)
	data_all.add_to_scene()

	# No boundary
	cell_vis = data_all.source.data.cell.vis
	#cell_vis.render_cell = False
	cell_vis.line_width = 0.5
	##cell_vis.rendering_color = (1.0, 1.0, 1.0)
	
	
	vp = Viewport()
	vp.type = Viewport.Type.Perspective
	#vp.fov = math.radians(40.0)
	#vp.camera_pos = (250+60, 60, 60)
	vp.fov = math.radians(10.0)
	vp.camera_pos = (850, 60, 60)
	vp.camera_dir = (-1, 0, 0)
	
	'''
	dist = 850/np.sqrt(3)
	vp.camera_pos = (60+dist, 60+dist, 60+dist)
	vp.camera_dir = (-1, -1, -1)
	'''
	filename = os.path.join(dir_imgs, filename_output+'_{}.png'.format( str(time_frame).zfill(4)) )
	vp.render_image(size=(800,800), filename=filename, background=(1,1,1),frame=time_frame, renderer=OpenGLRenderer())
	
	
	return data_all




	
def get_ids_bead_shared_psd(multi_graph):
	
	species = 'PSD95'
	
	ids = [ i for i, attr in multi_graph.nodes('species') if attr == species ]
	
	edges_from_molecules = [ multi_graph.edges(id, keys=True, data=True) for id in ids ]
	connections_from_one_molecules = [ [e[3]['type_connection'] for e in es] for es in edges_from_molecules ]
	
	ids_both = [id for id, c in zip(ids, connections_from_one_molecules) if ('STG_PSD95' in c) and ('GluN2B_PSD95' in c)]
	ids_bead_both = [multi_graph.nodes[id]['ids_bead'] for id in ids_both]
	ids_bead_both = np.ravel(ids_bead_both).tolist()
	
	ids_all = [id for id, c in zip(ids, connections_from_one_molecules)]
	ids_bead_all = [multi_graph.nodes[id]['ids_bead'] for id in ids_all]
	ids_bead_all = np.ravel(ids_bead_all).tolist()
	
	return ids_bead_both, ids_bead_all
	
	
if __name__ == '__main__':
	
	# Target file
	'''
	target = 'CGSP' #  ['CGSP', 'CG','CPG','PG','SP'] # ,'SPG'
	filename_lammpstrj 	= '{}.lammpstrj'.format(target)
	dir_input       = 'mixtures'
	filename_edited = target
	filename_output = target
	'''

	target = '023' #  ['CGSP', 'CG','CPG','PG','SP'] # ,'SPG'
	filename_lammpstrj 	= 'R2_{}.lammpstrj'.format(target)
	dir_input       = 'conc_dependence'
	filename_edited = target
	filename_output = target

	
	# Shared part of initialization
	dir_lammpstrj    = os.path.join('..', 'lammpstrj3', dir_input)
	dir_edited_data  = os.path.join('data3',dir_input)
	dir_imgs         = os.path.join('imgs3',dir_input,'snapshot2')
	os.makedirs(dir_imgs, exist_ok=True)
	suffix = 'connectivity_graph'
	
	
	# Load graph
	d = utils.load(dir_edited_data, filename_edited, suffix)
	multi_graph = d['multi_graph']
	# Find ids_bead for the PSD95 shared by STG and GluN2B
	ids_bead_both, ids_bead_all = get_ids_bead_shared_psd(multi_graph)
	
	
	# Load lammpstrj data
	print("\n"+filename_lammpstrj)
	sampling_frame = utils.get_num_frames(dir_lammpstrj, filename_lammpstrj)
	data_all   = import_file(os.path.join(dir_lammpstrj, filename_lammpstrj), input_format= "lammps/dump" )
	
	#sys.exit(0)
	plot_snapshots(data_all, sampling_frame, dir_imgs, filename_output, ids_bead_both)
	
