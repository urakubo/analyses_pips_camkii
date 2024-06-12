
import os, sys, glob, pickle, pprint, copy
import numpy as np
import random
from scipy.spatial.transform import Rotation as R

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot, patches

import networkx as nx

import utils
import parameters as p
import colormap as c

import pyvista as pv
import trimesh

rng = np.random.default_rng(seed=13)
sphere = pv.Sphere(radius=0.3, phi_resolution=10, theta_resolution=10)
line   = pv.Line()


def square_yz(magnification):
	y = p.space[1]/2/magnification
	z = p.space[2]/2/magnification
	pointa = [0.0, -y,  z]
	pointb = [0.0, -y, -z]
	pointc = [0.0, y , -z]
	pointd = [0.0, y ,  z]
	return pv.Rectangle([pointa, pointb, pointc])


def plot_3D_pvista_CaMKII_interface_region(d):
	
	pl = pyvista.Plotter(window_size=[400,1000], shape=(3, 1), border=False)
	pl.subplot(0, 0)
	cond = d['condensate_CaMKII']['condensate_CaMKII_in_grid_mesh']
	plot_condensates_pyvista(pl, cond, color='green')
	#cond = d['condensate_CaMKII']['region_interface_CaMKII']
	cond = d['condensate_CaMKII']['region_interface_all']
	plot_condensates_pyvista(pl, cond, color='red')
	pl.add_mesh(utils.square_yz(), color='black', style='wireframe')
	pl.view_yz()
	pl.camera.roll -= 90
	pl.camera.Zoom(1)
	pl.show(interactive=True, auto_close=True) # off_screen = True
	#filename = os.path.join(dir_imgs, '{}_{}.png'.format(prefix, suffix))
	#pl.screenshot(filename)


def plot_3D_pvista_cond_CaMKII(r, locs_binding_beads_CaMKII, locs_hub_bead_CaMKII, vertices, lines, rot):
	
	
	# Prepare the figure
	pl = pv.Plotter(window_size=[1200,1200], shape=(1, 1), border=False)
	pl.subplot(0, 0)
	
	# Add lines from the hub beads to binding beads
	
	#'''
	lines = pv.PolyData(vertices, lines=lines)
	# pl.add_mesh(lines, lighting=False, color="black")
	pl.add_mesh(lines, lighting=False, color=c.green_universal_uint)
	#'''
	
	# Add the spheres of hubs and binding beads.
	binding_beads = pv.PolyData(locs_binding_beads_CaMKII) # shape (n_points, 3)
	hub_beads     = pv.PolyData(locs_hub_bead_CaMKII) # shape (n_points, 3)
	
	glyph_binding_beads = binding_beads.glyph(geom=sphere)
	glyph_hub_beads     = hub_beads.glyph(geom=sphere)
	
	lighting = True # True, False
	pl.add_mesh(glyph_binding_beads, lighting=lighting, color=c.green_universal_uint) # c.light_green_universal_uint
	pl.add_mesh(glyph_hub_beads    , lighting=lighting, color="black")
	
	
	# plot the contour of condensate	
	'''
	mesh = utils.generate_mesh(cond_CaMKII, flipz = False)
	mesh.vertices = rot.apply( mesh.vertices ) 
	pl.add_mesh(mesh, color='green', show_edges=False,  opacity=0.1)
	'''
	
	condensate = pv.Sphere(radius=r, phi_resolution=30, theta_resolution=30)
	pl.add_mesh(condensate, lighting=lighting, color=c.light_green_universal_uint, show_edges=False,  opacity=0.05)
	
	# Plot the bounding box
	pl.add_mesh( square_yz(magnification=2.0), color='black', style='wireframe')
	pl.add_title('Radius: {:.3f}'.format(r))
	
	# Show the image
	pl.set_background('white')
	pl.view_yz()
	pl.camera.roll -= 90
	pl.camera.Zoom(1.3)
	
	return pl
	
	
if __name__ == '__main__':
	
	# Target file
	
	
	'''
	dir_edited_data  =  'small_colony'
	prefix = '01_008'
	random_seed = 2
	prefix = '00_008'
	random_seed = 5
	num_samples = 20
	
	num_samples = 20
	prefix = '00_009'
	#prefix = '01_009'
	radius = {'00_008': 14.484, '00_009': 15.806, '01_008': 15.045, '01_009':15.606 }
	'''
	
	num_samples = 30
	dir_target      = 'CG_valency_length'
	dir_edited_data = os.path.join('data3', dir_target)

	prefix, random_seed = '12_002', 1
	prefix, random_seed = '12_006', 0
	#prefix, random_seed = '12_005', 2
	
	radius = utils.load(dir_edited_data, 'Radiuses', 'CaMKII_bead')
	print('radius: ' )
	pprint.pprint(radius)
	
	# Output
	dir_imgs  = os.path.join('imgs3', dir_target,'surface_tension_prof')
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	# Load graph
	print(prefix)
	suffix = 'connectivity_graph'
	d = utils.load(dir_edited_data, prefix, suffix)
	multi_graph = d['multi_graph']
	cond_CaMKII = d['condensate_CaMKII']['condensate_CaMKII_in_grid_mesh']
	
	
	
	# Pickup the largest cluster
	clusters = sorted(nx.connected_components(multi_graph), key=len, reverse=True)
	lengths_clusters = [len(c) for c in clusters]
	nodes_cluster = clusters[lengths_clusters.index(lengths_clusters[0])] # 0: max
	g_largest_cluster = nx.MultiGraph( multi_graph.subgraph(nodes_cluster) )
	
	
	targ = 'CaMKII'
	ids_CaMKII = [n for n, v in g_largest_cluster.nodes.items() if v['species'] == targ ]
	random.seed(random_seed)
	random.shuffle(ids_CaMKII)
	
	# Rotation
	r1 = R.from_euler('x', 80, degrees=True)
	r2 = R.from_euler('z', 10, degrees=True)
	rot = r2 * r1
	
	
	locs_hub_bead =[]
	locs_binding_beads =[]
	vertices = []
	lines    = []
	i = 0
	current_num = 0
	while (current_num < num_samples):
		pos = g_largest_cluster.nodes[ids_CaMKII[i]]['positions_grid_coord']
		pos = rot.apply(pos)
		x = np.abs( pos[0,0] )
		y = np.abs( pos[0,1] )
		z = np.abs( pos[0,2] )
		xx = pos[1:,0]
		if (abs(x) > 1): # abs(xx[0]) < 2
			i += 1
			continue
		else:
			locs_hub_bead.append(pos[0,:].tolist())
			locs_binding_beads.extend(pos[1:,:].tolist())
			id_end_last = len( vertices )
			vertices.extend(pos.tolist())
			id_end_current = len( vertices )
			lines.extend([[2, id_end_last, j] for j in range(id_end_last+1, id_end_current)])
			i += 1
			current_num += 1
			
	locs_hub_bead      = np.array( locs_hub_bead )
	locs_binding_beads = np.array( locs_binding_beads )
	vertices = np.array( vertices )
	lines    = np.ravel(np.array( lines ) )
	
	
	'''
	print('locs_hub_bead.shape     ', locs_hub_bead.shape )
	print('locs_binding_beads.shape ', locs_binding_beads.shape )
	print('vertices.shape     ', vertices.shape )
	print('lines.shape ', lines.shape )
	'''
	
	# pl = plot_3D_pvista_cond_CaMKII(cond_CaMKII, locs_binding_beads, locs_hub_bead, vertices, lines, rot )
	pl = plot_3D_pvista_cond_CaMKII(radius[prefix], locs_binding_beads, locs_hub_bead, vertices, lines, rot )
	
	filename = os.path.join(dir_imgs, '{}_.png'.format(prefix))
	pl.show(interactive=False, auto_close=True) # 
	pl.screenshot(filename)
	pl.show(interactive=True, auto_close=True) # off_screen = True
	

