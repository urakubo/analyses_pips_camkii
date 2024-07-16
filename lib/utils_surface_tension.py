
import os, sys, glob, pickle, pprint
import numpy as np
import math

import matplotlib
import matplotlib.pyplot as plt

import networkx as nx

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c


def arrange_graph_bar(ax, panel_size_x, panel_size_y):
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax_h, ax_w = ax.bbox.height, ax.bbox.width
	loc = ax.get_position()
	ax.set_position([loc.x0, loc.y0, panel_size_x, panel_size_y])
	ax.set_xlabel('Linker length')
	
def plot_graphs(ave_cos, pull_force, surface_tensions, pull_force_per_area):
	
	# Figure specification
	num_columns = 10
	num_rows    = 3
	fig = plt.figure(figsize=(20, 8)) # , tight_layout=False
	
	left, right, bottom, top = 0.0, 0.95, 0.10, 0.99
	wspace, hspace = 0.2, 0.1
	plt.subplots_adjust(left, bottom, right, top, wspace, hspace)
	panel_size = 0.15
	
	# Plot concs in condensates
	column = 3
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*1 )
	ax.set_title('Average cosine')
	ax.bar(*zip(*ave_cos.items()), width=0.6 , color='gray' ) # , color=colormap_conc_bargraph
	ymin = np.min([min(ave_cos.values()), 0])
	ax.set_ylim([ymin,1.2*max(ave_cos.values())])
	ax.set_ylabel('E[cos(theta)]')
	arrange_graph_bar(ax, panel_size/5, panel_size)
	
	column = 5
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*1 )
	ax.set_title('Contraction force')
	ax.bar(*zip(*pull_force.items()), width=0.6, color='gray' ) # , color=colormap_conc_bargraph
	arrange_graph_bar(ax, panel_size/5, panel_size)
	
	column = 6
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*1 )
	ax.set_title('Surface tension')
	ax.bar(*zip(*surface_tensions.items()), width=0.6, color='gray' ) # , color=colormap_conc_bargraph
	arrange_graph_bar(ax, panel_size/5, panel_size)
	
	column = 7
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*1 )
	ax.set_title('Contraction force per area')
	ax.bar(*zip(*pull_force_per_area.items()), width=0.6, color='gray' ) # , color=colormap_conc_bargraph
	arrange_graph_bar(ax, panel_size/5, panel_size)
	
	return fig
	
	
def get_angle(x, y):
	dot_xy = np.dot(x, y)
	norm_x = np.linalg.norm(x)
	norm_y = np.linalg.norm(y)
	cos = dot_xy / (norm_x*norm_y)
	#rad = np.arccos(cos)
	rad = np.emath.arccos(cos)
	#theta = rad * 180 / np.pi
	#return theta
	return rad
	
	
def flatten(sequence):
    result = []
    for item in sequence:
        if isinstance(item, (list, tuple, range, dict, set, frozenset)):
            result.extend(flatten(item))
        else:
            result.append(item)
    return result
	
	
def plot_polar_scatter(angles, distances_to_hub, max_legnth, dir_imgs, prefix):
	
	fig = plt.figure(figsize=(2, 2))
	ax = fig.add_subplot(111, projection='polar')
	ax.scatter(angles, distances_to_hub, color='k', s=0.5)
	ax.set_title('Polar coordinates',fontsize=18)
	
	ax.set_title('{}'.format(prefix) )
	ax.set_theta_offset(np.pi / 2.0 * 3)
	ax.set_rlabel_position(180-20)
	
	fig.savefig( os.path.join(dir_imgs, '{}_polar_scatter.svg'.format( prefix ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}_polar_scatter.png'.format( prefix ) ) , dpi=150)
	plt.show()
	
	return
	
	
def plot_polar_histogram(angles, distances_to_hub, max_linker_length, dir_imgs, prefix):
	
	
	rbins = np.arange(0, max_linker_length+1, 1)
	abins = np.linspace(0,2*np.pi, 30)
	#abins = np.linspace(0,2*np.pi, 20)
	abins_s = np.hstack( [abins[1:], abins[0]] )
	# https://en.wikipedia.org/wiki/Spherical_segment
	normal_abins = 2 * np.pi * np.abs( np.cos(abins[1:]) - np.cos(abins[:-1]) )
	normal_abins = normal_abins.reshape(normal_abins.shape[0],-1)
	normal_rbins = 1/3 * ( rbins[1:]**3 - rbins[:-1]**3 )
	
	hist, _, _ = np.histogram2d(angles, distances_to_hub, bins=(abins, rbins))
	
	hist = hist / normal_abins / normal_rbins / angles.shape[0]
	hist = hist + np.flipud(hist)
	A, R = np.meshgrid(abins, rbins)
	fig, ax = plt.subplots(figsize=(2, 2), subplot_kw=dict(projection="polar"))
	pc = ax.pcolormesh(A, R, hist.T, cmap="Greys", vmin=0)  # , vmax = 0.0005, vmax = 0.0005
	ax.set_rlim(0, max_linker_length)
	ax.set_title('{}'.format(prefix) )
	#fig.colorbar(pc)
	ax.set_theta_offset(np.pi / 2.0 * 3)
	ax.set_rlabel_position(180-20)
	
	fig.savefig( os.path.join(dir_imgs, '{}_polar_hist.svg'.format( prefix ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}_polar_hist.png'.format( prefix ) ) , dpi=150)
	plt.show()
	
	return fig, ax
	
	
def get_properties_beads_CaMKII(g_largest_cluster):
	
	CaMKII_binding_site = {}
	for id, v in g_largest_cluster.nodes.items():
		if v['species'] == 'CaMKII':
			# Make CaMKII hub beads in the new graphs
			id_CaMKII_hub = np.nonzero(v['types_bead'] == p.subunits['CaMKII hub']['id'])
			position_CaMKII_hub = np.ravel( v['positions_grid_coord'][id_CaMKII_hub, :] )
			#id_bead       = v['ids_bead'][0][id_CaMKII_hub[0][0]]
			
			# Make CaMKII interaction beads
			id_CaMKII_binding = np.nonzero(v['types_bead'] == p.subunits['CaMKII binding site']['id'])
			position_CaMKII_binding = v['positions_grid_coord'][id_CaMKII_binding, :][0]
			ids_bead       = v['ids_bead'][0][id_CaMKII_binding]
			
			#print('position_CaMKII_binding: ', position_CaMKII_binding, ', ids_bead: ',ids_bead)
			for id_bead_, loc in zip(ids_bead, position_CaMKII_binding):
				vector_to_hub = position_CaMKII_hub - loc
				CaMKII_binding_site[id_bead_] = {}
				CaMKII_binding_site[id_bead_]['loc'] = loc
				CaMKII_binding_site[id_bead_]['connect'] = 0
				CaMKII_binding_site[id_bead_]['id_molecule'] = id
				CaMKII_binding_site[id_bead_]['vector_to_hub'] = vector_to_hub
				CaMKII_binding_site[id_bead_]['distance_to_hub'] = np.linalg.norm( vector_to_hub ) / np.sqrt(3)
				
				CaMKII_binding_site[id_bead_]['angle_from_condensate_center'] = get_angle(-loc, vector_to_hub)
				
				CaMKII_binding_site[id_bead_]['angle_from_hub'] = get_angle(position_CaMKII_hub, vector_to_hub)
				
				#CaMKII_binding_site[id_bead_]['angle'] = angle(vector_to_hub, -loc)
	
	return CaMKII_binding_site
	
	
	
def get_properties_beads_CaMKII_tmp(g_largest_cluster):
	
	CaMKII_binding_site = {}
	for id, v in g_largest_cluster.nodes.items():
	
		if v['species'] == 'CaMKII':
			# Make CaMKII hub beads in the new graphs
			id_CaMKII_hub = np.nonzero(v['types_bead'] == p.subunits['CaMKII hub']['id'])
			position_CaMKII_hub = np.ravel( v['positions_grid_coord'][id_CaMKII_hub, :] )
			#id_bead       = v['ids_bead'][0][id_CaMKII_hub[0][0]]
			if np.linalg.norm( position_CaMKII_hub ) / np.sqrt(3) < 8:
				# Make CaMKII interaction beads
				id_CaMKII_binding = np.nonzero(v['types_bead'] == p.subunits['CaMKII binding site']['id'])
				position_CaMKII_binding = v['positions_grid_coord'][id_CaMKII_binding, :][0]
				ids_bead       = v['ids_bead'][0][id_CaMKII_binding]
				
				#print('position_CaMKII_binding: ', position_CaMKII_binding, ', ids_bead: ',ids_bead)
				for id_bead_, loc in zip(ids_bead, position_CaMKII_binding):
					vector_to_hub = position_CaMKII_hub - loc
					CaMKII_binding_site[id_bead_] = {}
					CaMKII_binding_site[id_bead_]['loc'] = loc
					CaMKII_binding_site[id_bead_]['connect'] = 0
					CaMKII_binding_site[id_bead_]['id_molecule'] = id
					CaMKII_binding_site[id_bead_]['vector_to_hub'] = vector_to_hub
					CaMKII_binding_site[id_bead_]['distance_to_hub'] = np.linalg.norm( vector_to_hub ) / np.sqrt(3)
					
					CaMKII_binding_site[id_bead_]['angle_from_condensate_center'] = get_angle(-loc, vector_to_hub)
					
					CaMKII_binding_site[id_bead_]['angle_from_hub'] = get_angle(position_CaMKII_hub, vector_to_hub)
					
					#CaMKII_binding_site[id_bead_]['angle'] = angle(vector_to_hub, -loc)
	
	return CaMKII_binding_site
	
	
def calc_angle_and_distance_to_hub(d):
	
	# Obtain graph
	multi_graph = d['multi_graph']
	# cond_CaMKII = d['condensate_CaMKII']['condensate_CaMKII_in_grid_mesh']
	
	# Select the largest cluster
	nodes_clusters = sorted(nx.connected_components(multi_graph), key=len, reverse=True)
	g_largest_cluster = nx.MultiGraph( multi_graph.subgraph(nodes_clusters[0]) )
	
	# Get the ids of CaMKII binding beads of the CaMKII condensate.
	#CaMKII_binding_sites = get_properties_beads_CaMKII(g_largest_cluster)
	CaMKII_binding_sites = get_properties_beads_CaMKII_tmp(g_largest_cluster)
	
	angles_from_condensate_center = [ v['angle_from_condensate_center'] for v in CaMKII_binding_sites.values() ]
	angles_from_hub = [ v['angle_from_hub'] for v in CaMKII_binding_sites.values() ]
	distance_to_hub = [v['distance_to_hub'] for v in CaMKII_binding_sites.values()]
	
	angles_from_condensate_center = np.array( angles_from_condensate_center )
	angles_from_hub = np.array( angles_from_hub )
	distance_to_hub = np.array( distance_to_hub )
	
	return  angles_from_condensate_center, angles_from_hub, distance_to_hub
	
	
def calc_contraction_force(angles, distance_to_hub, max_linker_length, radius_condensate):
	
	streched_linkers = distance_to_hub >= max_linker_length - 1 # - 1 -1/np.sqrt(3)##############
	stretched_linkers_mult_cos = np.sum( streched_linkers * np.cos(angles))
	
	ave_cos           = np.sum( np.cos(angles) ) / distance_to_hub.shape[0]
	contraction_force = stretched_linkers_mult_cos
	surface_tension   = stretched_linkers_mult_cos / (2*np.pi*radius_condensate)
	contraction_force_per_area = stretched_linkers_mult_cos / (4*np.pi*radius_condensate*radius_condensate)
	
	return  ave_cos, contraction_force, surface_tension, contraction_force_per_area
	
	
	
	