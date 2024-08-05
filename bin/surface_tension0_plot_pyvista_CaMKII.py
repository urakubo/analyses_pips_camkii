
import os, sys, glob, pickle, pprint, copy
import numpy as np
import random
from scipy.spatial.transform import Rotation as R

import networkx as nx
import pyvista as pv
import trimesh
import copy

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_pyvista as utils_pyvista
from specification_datasets import SpecDatasets


rng = np.random.default_rng(seed=13)
sphere = pv.Sphere(radius=0.6, phi_resolution=10, theta_resolution=10)
sphere_small = pv.Sphere(radius=0.3, phi_resolution=10, theta_resolution=10)
sphere_large = pv.Sphere(radius=1.0, phi_resolution=10, theta_resolution=10)
line   = pv.Line()


def plot_3D_pvista_cond_CaMKII(radius, locs_hub_bead_CaMKII, locs_binding_beads_CaMKII, vertices, lines, rot, locs_hub_bead_all):
	
	
	# Prepare the figure
	pl = pv.Plotter(window_size=[1200,1200], shape=(1, 1), border=False)
	pl.subplot(0, 0)
	
	# Add lines from the hub beads to binding beads
	
	#'''
	lines = pv.PolyData(vertices, lines=lines)
	pl.add_mesh(lines, lighting=False, color="black")
	pl.add_mesh(lines, lighting=False, color=c.green_universal_uint)
	#'''
	
	
	'''
	gray = 0.6
	hub_bead_all  = pv.PolyData(locs_hub_bead_all) # shape (n_points, 3)
	glyph_hub_bead_all  = hub_bead_all.glyph(geom=sphere_small)
	pl.add_mesh(glyph_hub_bead_all , lighting=False, color=[gray, gray, gray]) # , ambient=0.6
	'''
	
	
	# Center
	gray = 1.0
	hub_bead_all  = pv.PolyData([0,0,0]) # shape (n_points, 3)
	glyph_hub_bead_all  = hub_bead_all.glyph(geom=sphere_large)
	pl.add_mesh(glyph_hub_bead_all , lighting=False, color=[gray, gray, gray] , ambient=0.6) #
	
	
	# Add the spheres of hubs and binding beads.
	hub_beads     = pv.PolyData(locs_hub_bead_CaMKII) # shape (n_points, 3)
	binding_beads = pv.PolyData(locs_binding_beads_CaMKII) # shape (n_points, 3)
	
	glyph_hub_beads     = hub_beads.glyph(geom=sphere)
	glyph_binding_beads = binding_beads.glyph(geom=sphere)
	
	lighting = True # True, False
	pl.add_mesh(glyph_binding_beads, lighting=lighting, color=c.green_universal_uint) # c.light_green_universal_uint
	pl.add_mesh(glyph_hub_beads    , lighting=lighting, color="black")
	
	
	# plot the contour of condensate	
	'''
	mesh = utils.generate_mesh(cond_CaMKII, flipz = False)
	mesh.vertices = rot.apply( mesh.vertices ) 
	pl.add_mesh(mesh, color='green', show_edges=False,  opacity=0.1)
	'''
	
	condensate = pv.Sphere(radius=radius, phi_resolution=30, theta_resolution=30)
	pl.add_mesh(condensate, lighting=lighting, color=c.light_green_universal_uint, show_edges=False,  opacity=0.15)
	
	# Plot the bounding box
	pl.add_mesh( utils_pyvista.square_yz_mag(magnification=2.0), color='black', style='wireframe')
	pl.add_title('Radius: {:.3f} l.s.'.format(radius / np.sqrt(3)))
	
	# Show the image
	pl.set_background('white')
	pl.view_yz()
	pl.camera.roll -= 90
	pl.camera.Zoom(1.3)
	#pl.camera.focal_point = (100, 0.0, 0.0)
	print('pl.camera.position ', pl.camera.position)
	print('pl.camera.focal_point ', pl.camera.focal_point)
	print('pl.camera.view_angle ', pl.camera.view_angle)
	#pl.camera.position = (70, 0.0, 0.0)
	#pl.camera.view_angle = 60
	
	return pl
	
	
def pickup_CaMKII_samples(g_largest_cluster, num_samples, random_seed, rot):
	
	locs_hub_bead_all = []
	
	
	
	
	locs_hub_bead =[]
	locs_binding_beads =[]
	vertices = []
	lines    = []
	
	ids_CaMKII = [n for n, v in g_largest_cluster.nodes.items() if v['species'] == 'CaMKII' ]
	random.seed(random_seed)
	random.shuffle(ids_CaMKII)
	
	
	for id in ids_CaMKII:
		pos = g_largest_cluster.nodes[id]['positions_grid_coord']
		pos = rot.apply(pos)
		x = pos[0,0]
		y = pos[0,1]
		z = pos[0,2]
		r = np.sqrt(x*x+y*y+z*z)
		if (abs(x) < 4): 
			locs_hub_bead_all.append(pos[0,:].tolist())
		
		
		
	N_regularity = 4
	accum_partitoned_num = np.array( np.ones(num_samples * N_regularity), dtype='int' )
	ids_target = [i for i in range(N_regularity*num_samples) if np.mod(i,N_regularity) == N_regularity-2] # N_regularity-1
	accum_partitoned_num[ids_target] = 0
	#ids_no_show = [i for i in range(num_samples * N_regularity) if i > num_samples * N_regularity / 2-2]
	ids_no_show = [i for i in range(num_samples * N_regularity) if (i < num_samples * N_regularity / 2)]
	accum_partitoned_num[ids_no_show] = 1
	print('accum_partitoned_num.shape ', accum_partitoned_num.shape)
	i  = 0
	while (1):
		i += 1
		pos = g_largest_cluster.nodes[ids_CaMKII[i]]['positions_grid_coord']
		pos = rot.apply(pos)
		x = pos[0,0]
		y = pos[0,1]
		z = pos[0,2]
		r = np.sqrt(y*y+z*z)
		
		if (abs(x) > 2.5) or (r < 14) : #  or (r > 12 *1.414)
			continue
		
		angle = (np.arctan2(y, z)+np.pi) / 2 /np.pi * (num_samples * N_regularity)
		id_angle = np.floor(angle).astype('int')
		#print( id_angle, end=' ' )
		if np.any(accum_partitoned_num[id_angle] >= 1):
			continue
		accum_partitoned_num[id_angle] += 1
		
		# current_num += 1
		locs_hub_bead.append(pos[0,:].tolist())
		locs_binding_beads.extend(pos[1:,:].tolist())
		id_end_last = len( vertices )
		vertices.extend(pos.tolist())
		id_end_current = len( vertices )
		lines.extend([[2, id_end_last, j] for j in range(id_end_last+1, id_end_current)])
		
		
		# print('accum_partitoned_num', accum_partitoned_num)
		if np.all(accum_partitoned_num >= 1):
			break
	
	locs_hub_bead_all  = np.array( locs_hub_bead_all )
	locs_hub_bead      = np.array( locs_hub_bead )
	locs_binding_beads = np.array( locs_binding_beads )
	vertices = np.array( vertices )
	lines    = np.ravel(np.array( lines ) )
	
	return locs_hub_bead, locs_binding_beads, vertices, lines, locs_hub_bead_all
	
	
	
class PlotPyvistaCaMKII(SpecDatasets):
	def __init__( self ):
		
		pass
		
	def plot_save( self, prefix = '12_006', random_seed = 0, num_samples = 30 ):
		
		dir_imgs = os.path.join(self.dir_imgs_root,'surface_tension_prof')
		os.makedirs(dir_imgs, exist_ok=True)
		print(prefix)
		
		
		## Load radius
		radiuses = utils.load(self.dir_edited_data, 'radiuses', 'CaMKII')
		radius   = radiuses[prefix]
		
		
		## Load graph
		suffix = 'connectivity_graph'
		d = utils.load(self.dir_edited_data, prefix, suffix)
		multi_graph = d['multi_graph']
		cond_CaMKII = d['condensate_CaMKII']['condensate_CaMKII_in_grid_mesh']
		
		
		## Pickup the largest cluster
		nodes_clusters = sorted(nx.connected_components(multi_graph), key=len, reverse=True)
		g_largest_cluster = nx.MultiGraph( multi_graph.subgraph(nodes_clusters[0]) )
		
		
		## Rotation
		r1 = R.from_euler('x', 80, degrees=True)
		r2 = R.from_euler('z', 10, degrees=True)
		rot = r2 * r1		
		
		
		## Pickup CaMKII samples
		locs_hub_bead, locs_binding_beads, vertices, lines, locs_hub_bead_all = \
			pickup_CaMKII_samples(g_largest_cluster, num_samples, random_seed, rot)
		
		
		'''
		print('locs_hub_bead.shape     ', locs_hub_bead.shape )
		print('locs_binding_beads.shape ', locs_binding_beads.shape )
		print('vertices.shape     ', vertices.shape )
		print('lines.shape ', lines.shape )
		'''
		print('locs_hub_bead_all.shape ', locs_hub_bead_all.shape )
		
		
		## Plot and save
		
		pl = plot_3D_pvista_cond_CaMKII(radius, locs_hub_bead, locs_binding_beads, vertices, lines, rot, locs_hub_bead_all )
		
		filename = os.path.join(dir_imgs, '{}_.png'.format(prefix))
		#pl.show(interactive=True, auto_close=False)
		pl.show(interactive=False, auto_close=True)
		pl.screenshot(filename)
		
		
		
if __name__ == '__main__':
	
	## Target file definition
	
	num_samples = 30
	prefix, random_seed = '12_002', 1
	prefix, random_seed = '12_006', 0
	#prefix, random_seed = '12_005', 2
	
	
	obj = PlotPyvistaCaMKII()
	obj.CG_valency_length()
	obj.plot_save( prefix = prefix, random_seed = random_seed, num_samples = 10 )
	
	
