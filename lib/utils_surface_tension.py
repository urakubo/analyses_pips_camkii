
import os, sys, glob, pickle, pprint
import numpy as np
import copy

# import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1


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
	
	
def plot_polar_scatter(angles, distances_to_hub, max_legnth, dir_imgs, prefix, color = 'k', center_direction = True):
	
	fig = plt.figure(figsize=(2, 2))
	ax = fig.add_subplot(111, projection='polar')
	ax.scatter(angles, distances_to_hub, color=color, s=0.5, linewidths=0, marker='o')
	ax.set_title('Polar coordinates',fontsize=18)
	
	ax.set_title('{}'.format(prefix) )
	ax.set_theta_offset(np.pi / 2.0 * 3)
	ax.set_rlabel_position(180+30)
	
	
	if center_direction == True:
		theta_ticks      = np.linspace(0, 315, 8)
		theta_ticklabels = ['0°', '45°', '90°', '135°', '180°', \
			'\N{MINUS SIGN}135°', '\N{MINUS SIGN}90°', '\N{MINUS SIGN}45°']
		ax.set_thetagrids(theta_ticks, labels=theta_ticklabels)
	
	fig.savefig( os.path.join(dir_imgs, '{}_polar_scatter.svg'.format( prefix ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}_polar_scatter.png'.format( prefix ) ) , dpi=150)
	plt.show()
	
	return
	
	
	
	
def plot_polar_histogram(angles, \
		distances_to_hub, \
		max_linker_length, \
		dir_imgs, \
		prefix, \
		cmap="Greys", step_color='k', \
		center_direction = True):
	
	ymax = 100
	
	# Get a histogram for polar graph
	rbins = np.arange(0, max_linker_length+1, 1)
	
	if center_direction == False:
		abins = np.linspace(0,2*np.pi, 33) # must be a odd number
	else:
		abins = np.linspace(-np.pi,np.pi, 33)
	
	hist, _, _ = np.histogram2d(angles, distances_to_hub, bins=(abins, rbins))
	hist = hist + np.flipud(hist) # 0-180 + -0-180
	
	# Sum it for angle or length only
	hist_angle  = copy.deepcopy( hist )
	hist_length = copy.deepcopy( hist )
	
	## Normalization
	## https://en.wikipedia.org/wiki/Spherical_segment
	normal_abins = 2 * np.pi * np.abs( np.cos(abins[1:]) - np.cos(abins[:-1]) )
	normal_abins = normal_abins.reshape(normal_abins.shape[0],-1)
	normal_rbins = 1/3 * ( rbins[1:]**3 - rbins[:-1]**3 )
	hist = hist / normal_abins / normal_rbins / angles.shape[0]
	hist = hist / np.max(hist[:,1:]) * ymax
	A, R = np.meshgrid(abins, rbins)
	
	hist_angle = hist_angle / normal_abins
	hist_angle = np.sum(hist_angle, axis=1)
	hist_angle = hist_angle/np.max(hist_angle)*ymax
	
	hist_length = hist_length / normal_rbins
	hist_length = np.sum(hist_length, axis=0)
	hist_length = hist_length/np.max(hist_length)*ymax
	
	'''
	print('hist_angle.shape ', hist_angle.shape)
	print('normal_abins.shape ', normal_abins.shape)
	print('hist_angle')
	print(hist_angle)
	'''
	
	
	##
	## Polar plot 
	##
	fig = plt.figure(figsize=(7, 2.2))
	fig.suptitle('{}'.format(prefix) )
	ax1 = fig.add_axes((0, 0.1, 0.5, 0.7), projection='polar') # (left, bottom, width, height)
	
	
	pc = ax1.pcolormesh(A, R, hist.T, cmap= cmap, vmin=0, vmax= ymax )  # , vmax = 0.0005, vmax = 0.0005
	ax1.set_rlim(0, max_linker_length)
	plt.colorbar(pc, ax=ax1)

	#ax1.set_theta_offset(np.pi / 2.0 * 3)
	ax1.set_theta_offset(-np.pi / 2.0)
	ax1.set_rlabel_position(180+90)
	
	
	if center_direction == True:
		theta_ticks      = np.linspace(0, 315, 8)
		theta_ticklabels = ['0°', '45°', '90°', '135°', '180°', \
			'\N{MINUS SIGN}135°', '\N{MINUS SIGN}90°', '\N{MINUS SIGN}45°']
		ax1.set_thetagrids(theta_ticks, labels=theta_ticklabels)
		
	##
	## Angle plot
	##
	ax2 = fig.add_axes((0.6, 0.2, 0.22, 0.15)) # (left, bottom, width, height)
	
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.hist( abins[:-1]/np.pi*180, abins/np.pi*180, weights=hist_angle, color=step_color )
	#ax2.step(abins[1:]/np.pi*180, hist_angle, color=step_color )
	ax2.set_xlabel('delta/sigma (deg)')
	ax2.set_ylabel('Probability\n(%max)')
	
	ax2.set_yticks([0,100])
	ax2.set_ylim([0,ymax*1.2])
	
	if center_direction == False:
		ax2.set_xlim([0,360])
		ax2.set_xticks([0,90,180,270,360])
		ax2.invert_xaxis()
	else:
		ax2.set_xlim([-180,180])
		ax2.set_xticks([-180,-90,0,90,180])
	
	
	ax3 = fig.add_axes((0.6, 0.6, 0.07, 0.15)) # (left, bottom, width, height)
	
	#ax3.set_title('{}'.format(prefix) )
	ax3.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.hist( hist_length, rbins, color=step_color )
	ax3.hist( rbins[:-1], rbins, weights=hist_length, color=step_color )
	#ax3.step(rbins[1:], hist_length, color=step_color, where='mid' )
	ax3.set_xlabel('length (l.s.)')
	ax3.set_ylabel('Probability\n(%max)')
	ax3.set_xlim([np.min(rbins),np.max(rbins)])
	ax3.set_ylim([0,ymax*1.2])
	ax3.set_xticks([np.min(rbins),np.max(rbins)])
	ax3.set_yticks([0,100])
	##
	
	fig.savefig( os.path.join(dir_imgs, '{}_polar_hist.svg'.format( prefix ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}_polar_hist.png'.format( prefix ) ) , dpi=150)
	plt.show()
	
	return fig, ax1
	
	
def get_properties_beads_CaMKII_surrogate(g_largest_cluster):
	
	ids_CaMKII = []
	for id, v in g_largest_cluster.nodes.items():
		if v['species'] == 'CaMKII':
			ids_CaMKII.append(id)
	
	import random
	random.shuffle(ids_CaMKII)
	
	CaMKII_binding_site = {}
	for i in range(len(ids_CaMKII) - 1):
		id   = ids_CaMKII[i]
		id_p = ids_CaMKII[i+1]
		v    = g_largest_cluster.nodes[ids_CaMKII[i]]
		v_p  = g_largest_cluster.nodes[ids_CaMKII[i+1]]
		
		'''
		rng = np.random.RandomState(123)
		theta = rng.uniform(-1,1)
		phi = rng.uniform(0,2*np.pi)
		r = rng.uniform(0, 10)
		x=r**(1/3)*(1-theta**2)*np.cos(phi)
		y=r**(1/3)*(1-theta**2)*np.sin(phi)
		z=r**(1/3)*theta
		position_CaMKII_hub_p = np.array([x,y,z])
		'''
		
		# Make CaMKII hub beads in the new graphs
		id_CaMKII_hub = np.nonzero(v['types_bead'] == p.subunits['CaMKII hub']['id'])
		position_CaMKII_hub = np.ravel( v['positions_grid_coord'][id_CaMKII_hub, :] )
		position_CaMKII_hub_p = np.ravel( v_p['positions_grid_coord'][id_CaMKII_hub, :] )
		position_hub_shift = - position_CaMKII_hub + position_CaMKII_hub_p
		#position_hub_shift = - 2* position_CaMKII_hub
		#print(position_hub_shift)
		
		
		#id_bead       = v['ids_bead'][0][id_CaMKII_hub[0][0]]
		# Make CaMKII interaction beads
		id_CaMKII_binding = np.nonzero(v['types_bead'] == p.subunits['CaMKII binding site']['id'])
		position_CaMKII_binding = v['positions_grid_coord'][id_CaMKII_binding, :][0]
		position_CaMKII_binding_surrogate = position_CaMKII_binding + position_hub_shift
		ids_bead       = v['ids_bead'][0][id_CaMKII_binding]
		
		#print('position_CaMKII_binding: ', position_CaMKII_binding, ', ids_bead: ',ids_bead)
		for id_bead_, loc, loc_surrogate in zip(ids_bead, position_CaMKII_binding, position_CaMKII_binding_surrogate):
			vector_to_hub = position_CaMKII_hub - loc
			CaMKII_binding_site[id_bead_] = {}
			CaMKII_binding_site[id_bead_]['loc'] = loc
			CaMKII_binding_site[id_bead_]['connect'] = 0
			CaMKII_binding_site[id_bead_]['id_molecule'] = id
			CaMKII_binding_site[id_bead_]['vector_to_hub'] = vector_to_hub
			CaMKII_binding_site[id_bead_]['distance_to_hub'] = np.linalg.norm( vector_to_hub ) / np.sqrt(3)
			
			CaMKII_binding_site[id_bead_]['angle_from_condensate_center'] = get_angle(-loc_surrogate, vector_to_hub)
			CaMKII_binding_site[id_bead_]['angle_from_hub'] = get_angle(position_CaMKII_hub, vector_to_hub)
			
	return CaMKII_binding_site
	
	
	
def get_properties_beads_CaMKII(g_largest_cluster):
	
	CaMKII_binding_site = {}
	for id, v in g_largest_cluster.nodes.items():
	
		if v['species'] == 'CaMKII':
			# Make CaMKII hub beads in the new graphs
			id_CaMKII_hub = np.nonzero(v['types_bead'] == p.subunits['CaMKII hub']['id'])
			position_CaMKII_hub = np.ravel( v['positions_grid_coord'][id_CaMKII_hub, :] )
			#id_bead       = v['ids_bead'][0][id_CaMKII_hub[0][0]]
			if 1: # np.linalg.norm( position_CaMKII_hub ) / np.sqrt(3) > 9:
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
					
	return CaMKII_binding_site
	
	
def calc_angle_and_distance_to_hub(d, surrogate=True):
	
	# Obtain graph
	multi_graph = d['multi_graph']
	# cond_CaMKII = d['condensate_CaMKII']['condensate_CaMKII_in_grid_mesh']
	
	# Select the largest cluster
	nodes_clusters = sorted(nx.connected_components(multi_graph), key=len, reverse=True)
	g_largest_cluster = nx.MultiGraph( multi_graph.subgraph(nodes_clusters[0]) )
	
	# Get the ids of CaMKII binding beads of the CaMKII condensate.
	if surrogate == True:
		CaMKII_binding_sites = get_properties_beads_CaMKII_surrogate(g_largest_cluster)
	else:
		CaMKII_binding_sites = get_properties_beads_CaMKII(g_largest_cluster)
	
	angles_from_condensate_center = [ v['angle_from_condensate_center'] for v in CaMKII_binding_sites.values() ]
	angles_from_hub = [ v['angle_from_hub'] for v in CaMKII_binding_sites.values() ]
	distance_to_hub = [v['distance_to_hub'] for v in CaMKII_binding_sites.values()]
	
	
	angles_from_condensate_center = np.array( angles_from_condensate_center )
	angles_from_hub = np.array( angles_from_hub )
	distance_to_hub = np.array( distance_to_hub )
	
	return  angles_from_condensate_center, angles_from_hub, distance_to_hub
	
	
def calc_contraction_force(angles, distance_to_hub, max_linker_length, radius_condensate):
	
	streched_linkers = distance_to_hub >= max_linker_length - 2 # - 1 -1/np.sqrt(3)##############
	
	streched_linkers = (distance_to_hub >= max_linker_length - 2) * (distance_to_hub - (max_linker_length - 2)) / 2
	
	stretched_linkers_mult_cos = np.sum( streched_linkers * np.cos(angles))
	
	ave_cos           = np.sum( np.cos(angles) ) / distance_to_hub.shape[0]
	contraction_force = stretched_linkers_mult_cos
	surface_tension   = stretched_linkers_mult_cos / (2*np.pi*radius_condensate)
	contraction_force_per_area = stretched_linkers_mult_cos / (4*np.pi*radius_condensate*radius_condensate)
	
	return  ave_cos, contraction_force, surface_tension, contraction_force_per_area
	
	
	
	