
# https://stackoverflow.com/questions/59821151/plot-the-dendrogram-of-communities-found-by-networkx-girvan-newman-algorithm
# https://qiita.com/s-wakaba/items/a93f03f27137cff4a26c

import networkx as nx
import matplotlib.pyplot as plt

import os, sys

import numpy as np

import utils
import parameters as p
import colormap as c

plt.rcParams.update(p.rc_param)

import os, glob
os.environ['OVITO_GUI_MODE'] = '1' # Request a session with OpenGL support

import numpy as np
import matplotlib.pyplot as plt
from ovito.io import import_file
from ovito.modifiers import SliceModifier, SelectTypeModifier, AssignColorModifier
from ovito.pipeline import ModifierInterface


import math
from ovito.vis import *
#import ovito.data
from ovito.pipeline import *

import utils
import colormap as c
import parameters as p

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


	
	
def plot_snapshots(data_all, num_frames, dir_imgs, filename_output, ids_col):
	
	'''
	for k, v in p.molecules_without_all.items():
		data_all.modifiers.append(SelectTypeModifier(types=set(v['id'])))
		data_all.modifiers.append(AssignColorModifier(color=c.cmap_universal_ratio[k] ))
	'''
	
	#data_all.modifiers.append(compute_particle_transparency)
	
	modifier = SetColorsById()
	modifier.ids_colors = ids_col
	data_all.modifiers.append(modifier)
	

	modifier = DeleteParticle()
	modifier.ids = list(ids_col.keys())
	data_all.modifiers.append(modifier)
	
	'''
	modifier = MakeTransparent()
	modifier.ids = list(ids_col.keys())
	data_all.modifiers.append(modifier)
	'''
	
	'''
	# Slice data
	modifier = SliceModifier()
	modifier.normal   = (-1.0, 0.0, 0.0)
	modifier.distance = 60
	data_all.modifiers.append(modifier)
	'''	
	
	
	
	
	data   = data_all.compute()
	data_all.add_to_scene()

	vp = Viewport()
	vp.type = Viewport.Type.Perspective
	#vp.fov = math.radians(40.0)
	#vp.camera_pos = (250+60, 60, 60)
	vp.fov = math.radians(10.0)
	vp.camera_pos = (850, 60, 60)
	vp.camera_dir = (-1, 0, 0)
	
	for time_frame in range( num_frames ):
		filename = os.path.join(dir_imgs, filename_output+'_{}.png'.format( str(time_frame).zfill(4)) )
		vp.render_image(size=(800,800), filename=filename, background=(1,1,1),frame=time_frame, renderer=OpenGLRenderer())
	
	return

def get_ratio_code(hex_code):
	hex_code   = hex_code.lstrip("#")
	ratio_code = [int(hex_code[i:i+2], 16)/255 for i in range(0, 6, 2)]
	return ratio_code

if __name__ == '__main__':

	## Init
	dir_target  =  'small_colony'
	prefix = '01_008'
	#prefix = '01_004'

	# lengths_clusters  [1503, 881, 699, 447, 274, 1, 1, 1, 1, 1]
	prefix = '01_004'
	nth_largest = 1

	# lengths_clusters  [845, 838, 793, 443, 372, 368, 1, 1, 1, 1]
	prefix = '00_004'
	nth_largest = 0 #0


	fig_title = '{}_{}'.format(nth_largest, prefix)
	print(fig_title)
	
	
	
	subdirs    = ['CG_con', 'CG_len9', 'CG_lin']
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(10)]
	filenames_lammpstrj = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_edited    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(3) for id_f in range(10)]
	
	
	filename_lammpstrj  = filenames_lammpstrj[ filenames_edited.index(prefix) ]
	
	
	dir_lammpstrj    = os.path.join('..','lammpstrj3', dir_target)
	dir_edited_data  = os.path.join('data3', dir_target)
	dir_imgs         = os.path.join('imgs3', dir_target, 'connectivity_movie')
	os.makedirs(dir_imgs, exist_ok=True)
	
	# Load graph.
	d2 = utils.load(dir_edited_data, prefix, 'cluster_dendrogram')
	
	Z          = d2['Z']
	labels     = d2['labels']
	node_order = d2['node_order']
	blocks     = d2['blocks'] 
	partition  = d2['partition'] 
	leaves_color_list  = d2['leaves_color_list'] 
	multi_graph_CaMKII = d2['multi_graph_CaMKII']
	
	# Load lammpstrj data
	num_frames = utils.get_num_frames(dir_lammpstrj, filename_lammpstrj)
	types, positions_grid_coord,ids_molecule, mc_step = utils.load_lammpstrj( dir_lammpstrj, filename_lammpstrj, num_frames )
	print('num_frames ', num_frames)
	
	
	#ids_molecule_target_CaMKII_hub = [multi_graph_CaMKII.nodes[i]['id_bead'] for i in node_order]
	ids_molecule_target_CaMKII_hub = [v['id_bead'] for i, v in multi_graph_CaMKII.nodes.items()]
	
	from cycler import cycler
	col = cycler( color=c.cols_list_ratio )
	cols_matlab = plt.rcParams['axes.prop_cycle'].by_key()['color']
	cols = {'C'+str(i+1): get_ratio_code(hex_code) for i, hex_code in enumerate(cols_matlab) }
	
	# ids_col = {multi_graph_CaMKII.nodes[i]['id_bead']: cols[c] for i, c in zip(node_order, leaves_color_list)}
		
	ids_col = {ii: cols[c] for i, c in zip(node_order, leaves_color_list) for ii in multi_graph_CaMKII.nodes[i]['id_bead_all']}
	
	# ids_untarget_molecules = list( set(ids_molecule).difference(set(ids_molecule_target_CaMKII_hub)) )
	
	
	#print('len(ids_molecule_target_CaMKII_hub) ', len(ids_molecule_target_CaMKII_hub))
	#print('len(ids_untarget_molecules) ', len(ids_untarget_molecules))
	
	data_all = import_file(os.path.join(dir_lammpstrj, filename_lammpstrj), input_format= "lammps/dump" )
	
	
	
	plot_snapshots(data_all, num_frames, dir_imgs, fig_title, ids_col)
	