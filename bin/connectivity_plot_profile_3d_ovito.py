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
from specification_datasets import SpecDatasets


def plot_snapshots(data_all, time_frame, dir_imgs, filename_output, ids_col):
	
	# Fix bounging box
	data_all.modifiers.append( utils_ovito.FixBoundingbox() )
	# Non boundary
	cell_vis = data_all.source.data.cell.vis
	#cell_vis.render_cell = False
	cell_vis.line_width = 0.5
	##cell_vis.rendering_color = (1.0, 1.0, 1.0)
	
	
	modifier = utils_ovito.SetColorsById()
	modifier.ids_colors = ids_col
	data_all.modifiers.append(modifier)
	
	
	modifier = utils_ovito.DeleteParticle()
	modifier.ids = list(ids_col.keys())
	data_all.modifiers.append(modifier)
	
	# Centering
	modifier = utils_ovito.CenteringEachFrameModifier()
	data_all.modifiers.append(modifier)
	
	
	# Slice data
	#'''
	modifier = SliceModifier()
	modifier.normal   = (1.0, 0.0, 0.0)
	modifier.distance = 60
	data_all.modifiers.append(modifier)
	#'''
	
	
	vp = Viewport()
	vp.type = Viewport.Type.Perspective
	vp.fov = np.radians(10.0)
	vp.camera_pos = (850, 60, 60)
	vp.camera_dir = (-1, 0, 0)
	
	data_all.add_to_scene()
	filename = os.path.join(dir_imgs, filename_output+'_{}.png'.format( str(time_frame).zfill(4)) )
	vp.render_image(size=(800,800), filename=filename, background=(1,1,1),frame=time_frame, renderer=OpenGLRenderer())
	data_all.remove_from_scene()
	
	return


	# Scale bar
	# https://ovito.org/manual/python/introduction/examples/overlays/scale_bar.html


def assign_color_for_each_bead(rcms, multi_graph_CaMKII):
	from cycler import cycler
	cols   = cycler( color=c.cols_list_ratio )
	colors = [c['color'] for p, c in zip(rcms, cols())]
	
	# ids_col = {id: c for ids_node, c in zip(rcm, colors) for id_node in ids_node for id in multi_graph_CaMKII.nodes[id_node]['id_bead_all']}
	ids_col = {}
	for ids_node, col in zip(rcms, colors):
		for id_node in ids_node:
			for id in multi_graph_CaMKII.nodes[id_node]['id_bead_all']:
				ids_col[id] = col
	return ids_col


class Plot3dOvitoConnectivity(SpecDatasets):
	def __init__( self ):
		pass
		
		
	def inspect( self ):
		for i, (f_lammpstrj, f_edited) in enumerate( zip(self.filenames_lammpstrj, self.filenames_edited) ):
			print('ID {}: {}, {}'.format(i, f_lammpstrj, f_edited))
		print()
		
	def plot_an_image( self, i ):
		
		filename_lammpstrj = self.filenames_lammpstrj[i]
		filename_edited    = self.filenames_edited[i]
		
		print('\n ID {}: {}, {}'.format(i, filename_lammpstrj, filename_edited))
		
		dir_imgs = os.path.join(self.dir_imgs_root, 'connectivity_3d_ovito' )
		os.makedirs(dir_imgs, exist_ok=True)
		
		
		# Load particles
		sampling_frame = utils.get_num_frames(self.dir_lammpstrj, filename_lammpstrj)
		data_all   = import_file(os.path.join(self.dir_lammpstrj, filename_lammpstrj), input_format= "lammps/dump" )
		
		
		# Load graph
		loadname_prefix = '{}_{}'.format(0, self.filename_edited )
		data = utils.load(self.dir_edited_data, loadname_prefix, 'cluster_dendrogram')
		
		#  Assign colors
		#color_list        = data['color_list']
		multi_graph_CaMKII = data['multi_graph_CaMKII']
		partition          = data['partition']
		ids_group = sorted(set(partition.values())) # List of group-ids
		rcms = [[k for k, v in partition.items() if v == i] for i in ids_group ] # List of the CaMKII-ids list
		ids_col = assign_color_for_each_bead(rcms, multi_graph_CaMKII)
		
		# print(ids_col)
		
		plot_snapshots(data_all, sampling_frame, dir_imgs, filename_edited, ids_col)
		print('filename_input: ', filename_lammpstrj)
		print('sampling_frame: ', sampling_frame )
		
		
if __name__ == '__main__':
	
	# I manually ran each one of them,
	# because I do not know how to fully reset the ovito visualization system.
	
	obj = Plot3dOvito()
	obj.boundary_conditions2() #  conc_dependence(), valency_length(), valency_length_CG(), boundary_conditions2()
	# obj.inspect()
	obj.run(1)
	
	
