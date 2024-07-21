#
# I referred the following web page.
# https://www.ovito.org/docs/current/python/introduction/pipelines.html
# 

import os, glob
os.environ['OVITO_GUI_MODE'] = '1' # Request a session with OpenGL support

import numpy as np
import networkx as nx


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


def assign_color_for_CaMKII_hub_beads(multi_graph):
	clusters = sorted(nx.connected_components(multi_graph), key=len, reverse=True)
	g_largest_cluster = nx.MultiGraph( multi_graph.subgraph(clusters[0]) )
	
	ids_CaMKII_hub = []
	ids_CaMKII_binding = []
	for id, v in g_largest_cluster.nodes.items():
		if v['species'] == 'CaMKII':
			# Make CaMKII hub beads in the new graphs
			id_CaMKII_hub  = np.nonzero(v['types_bead'] == p.subunits['CaMKII hub']['id'])
			id_hub_bead    = v['ids_bead'][0][id_CaMKII_hub[0][0]]
			ids_CaMKII_hub.append( id_hub_bead )
			# Make CaMKII interaction beads
			id_CaMKII_binding = np.nonzero(v['types_bead'] == p.subunits['CaMKII binding site']['id'])
			id_binding_bead  = v['ids_bead'][0][id_CaMKII_binding]
			ids_CaMKII_binding.append( id_binding_bead.tolist() )
	
	ids_CaMKII_hub     = utils.flatten(ids_CaMKII_hub)
	ids_CaMKII_binding = utils.flatten(ids_CaMKII_binding)
	
	ids_col_CaMKII_hub = {id :  (0,0,0) for id in ids_CaMKII_hub}
	ids_col_CaMKII_binding = {id : (3/255,175/255,122/255) for id in ids_CaMKII_binding}
	
	
	#print('ids_col_CaMKII_hub ', ids_col_CaMKII_hub)
	#print('ids_col_CaMKII_binding ', ids_col_CaMKII_binding)
	
	ids_col = {**ids_col_CaMKII_hub, **ids_col_CaMKII_binding}
	
	
	return ids_col
	


def assign_color_for_each_bead(partitioned, color_list, multi_graph_CaMKII):


	#import matplotlib.pyplot as plt
	from cycler import cycler
	#plt.rcParams["axes.prop_cycle"] = cycler( color=c.cols_list_ratio )
	#cols = plt.rcParams['axes.prop_cycle']
	#'''
	#from itertools import cycle
	#colors = cycle(prop_cycle.by_key()['color'])
	
	cols   = cycler( color=c.cols_list_ratio )
	#cols   = cycler( color=color_list )
	#color_codes = map('C{}'.format, cols(10))
	#colors = [c['color'] for p, c in zip(rcms+[0], cols())]
	#colors = colors[1:]
	#'''
	
	print('len(partitioned) ', len(partitioned))
	print('len(color_list) ', len(color_list))
	
	
	# ids_col = {id: c for ids_node, c in zip(rcm, colors) for id_node in ids_node for id in multi_graph_CaMKII.nodes[id_node]['id_bead_all']}
	ids_col = {}
	for ids_node, col in zip(partitioned, color_list):
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
		
	def plot_an_image( self, i, mode = 'clustering' ):
		
		filename_lammpstrj = self.filenames_lammpstrj[i]
		filename_edited    = self.filenames_edited[i]
		
		print('\n ID {}: {}, {}'.format(i, filename_lammpstrj, filename_edited))
		
		dir_imgs = os.path.join(self.dir_imgs_root, 'connectivity_3d_ovito' )
		os.makedirs(dir_imgs, exist_ok=True)
		
		
		# Load particles
		sampling_frame = utils.get_num_frames(self.dir_lammpstrj, filename_lammpstrj)
		data_all   = import_file(os.path.join(self.dir_lammpstrj, filename_lammpstrj), input_format= "lammps/dump" )
		
		
		
		#  Assign colors
		if mode == 'clustering':
			# Load graph
			loadname_prefix = '{}_{}'.format(0, filename_edited )
			data = utils.load(self.dir_edited_data, loadname_prefix, 'cluster_dendrogram')
			multi_graph_CaMKII = data['multi_graph_CaMKII']
			blocks             = data['blocks']
			color_list         = data['color_list']
			node_order         = data['node_order']
			
			partitioned = [node_order[s:e] for s,e  in blocks]
			
			
			# partition          = data['partition']
			#print('partition ', partition)
			# print('blocks     ', blocks)
			# print('node_order ', node_order)
			# ids_group = sorted(set(partition.values())) # List of group-ids
			# rcms = {i: [k for k, v in partition.items() if v == i] for i in ids_group } # List of the CaMKII-ids list
			
			ids_col = assign_color_for_each_bead(partitioned, color_list, multi_graph_CaMKII)
		elif mode == 'CaMKII_hub_beads':
			# Load graph
			data = utils.load(self.dir_edited_data, filename_edited,  'connectivity_graph')
			multi_graph = data['multi_graph']
			ids_col = assign_color_for_CaMKII_hub_beads(multi_graph)
		
		
		
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
	
	
