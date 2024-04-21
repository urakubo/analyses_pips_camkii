
import os, sys, glob, pickle, pprint
import numpy as np
import math

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot, patches

import networkx as nx

import utils
import parameters as p
import colormap as c

import pyvista
import trimesh

import itertools


plt.rcParams.update(p.rc_param)


def angle(x, y):
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
	
	
def plot_polar_scatter(angle_interface, distance_to_hub_interface, prefix):
	
	fig = plt.figure(figsize=(4, 4))
	ax = fig.add_subplot(111, projection='polar')
	ax.scatter(angle_interface, distance_to_hub_interface, color='k', s=0.5)
	ax.set_title('Polar coordinates',fontsize=18)
	
	ax.set_title('{}'.format(prefix) )
	ax.set_theta_offset(np.pi / 2.0 * 3)
	ax.set_rlabel_position(180-20)
	
	return fig, ax
	
	
	
def plot_polar_histogram(angle_interface, distance_to_hub_interface, prefix):
	
	rbins = np.arange(0, max(distance_to_hub_interface), 1)
	abins = np.linspace(0,2*np.pi, 60)
	abins_s = np.hstack( [abins[1:], abins[0]] )
	# https://en.wikipedia.org/wiki/Spherical_segment
	normal_abins = 2 * np.pi * np.abs( np.cos(abins[1:]) - np.cos(abins[:-1]) )
	normal_abins = normal_abins.reshape(normal_abins.shape[0],-1)
	normal_rbins = 1/3 * ( rbins[1:]**3 - rbins[:-1]**3 )
	
	hist, _, _ = np.histogram2d(angle_interface, distance_to_hub_interface, bins=(abins, rbins))
	
	hist = hist / normal_abins / normal_rbins / angle_interface.shape[0]
	hist = hist + np.flipud(hist)
	A, R = np.meshgrid(abins, rbins)
	fig, ax = plt.subplots(figsize=(4, 4), subplot_kw=dict(projection="polar"))
	pc = ax.pcolormesh(A, R, hist.T, cmap="Greys", vmin=0)  # 0.025 0.0015
	ax.set_title('{}'.format(prefix) )
	#fig.colorbar(pc)
	ax.set_theta_offset(np.pi / 2.0 * 3)
	ax.set_rlabel_position(180-20)
	return fig, ax
	
	
def get_properties_beads_CaMKII(g_largest_cluster, ids_CaMKII_binding_bead_select = None):
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
				if (ids_CaMKII_binding_bead_select is not None) and (id_bead_ in ids_CaMKII_binding_bead_select):
					vector_to_hub = position_CaMKII_hub - loc
					CaMKII_binding_site[id_bead_] = {}
					CaMKII_binding_site[id_bead_]['loc'] = loc
					CaMKII_binding_site[id_bead_]['connect'] = 0
					CaMKII_binding_site[id_bead_]['id_molecule'] = id
					CaMKII_binding_site[id_bead_]['vector_to_hub'] = vector_to_hub
					CaMKII_binding_site[id_bead_]['distance_to_hub'] = np.linalg.norm( vector_to_hub ) / np.sqrt(3)
					#CaMKII_binding_site[id_bead_]['angle'] = angle(vector_to_hub, -loc)
	
	return CaMKII_binding_site
	
	
if __name__ == '__main__':
	
	

	
	#'''
	dir_target  =  'small_colony'
	filenames = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(3) for id_f in range(10) ]
	filenames = ['00_008'] # 00: len-3, 01: len-9, 02: linear
	filenames = ['01_008'] # 00: len-3, 01: len-9, 02: linear
	length = {'00_008':3, '01_008':9 }
	#'''
	
	
	# Dataset 2:  Valency length
	#'''
	filenames = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
	dir_target  = 'valency_length'
	#filenames = ['12_{}'.format(str(id_f).zfill(3)) for id_f in range(7)] # '12_002', '12_004', '12_005', '12_006'
	filenames = ['12_000', '12_001', '12_002', '12_003', '12_004', '12_005', '12_006']
	filenames = ['04_000', '04_001', '04_002', '04_003', '04_004', '04_005', '04_006']
	filenames = ['06_000', '06_001', '06_002', '06_003', '06_004', '06_005', '06_006']
	length    = {	'04_000':1, '04_001':2, '04_002':3, '04_003':4, '04_004':5, '04_005':6 , '04_006':9, \
					'06_000':1, '06_001':2, '06_002':3, '06_003':4, '06_004':5, '06_005':6 , '06_006':9, \
					'12_000':1, '12_001':2, '12_002':3, '12_003':4, '12_004':5, '12_005':6 , '12_006':9  }
	
	#'''
	
	
	
	type_condensate = 'condensate_all_in_grid_mesh' # 'condensate_all_in_grid_mesh', 'condensate_CaMKII_in_grid_mesh'
	
	# Shared part of initialization
	dir_edited_data  = os.path.join('data3', dir_target)
	os.makedirs(dir_edited_data, exist_ok=True)
	dir_imgs = os.path.join('imgs3', dir_target,'surface_tension')
	os.makedirs(dir_imgs, exist_ok=True)
	
	#angles_interface     = np.zeros( len(filenames), dtype='float' )
	#angles_all           = np.zeros( len(filenames), dtype='float' )
	#cluster_coefficients = np.zeros( len(filenames), dtype='float' )
	
	
	for ii, prefix in enumerate( filenames ):
		print()
		print(prefix)
		
		# Load graph and mesh
		d = utils.load(dir_edited_data, prefix, 'connectivity_graph')
		multi_graph = d['multi_graph']
		cond_CaMKII = d['condensate_CaMKII']['condensate_CaMKII_in_grid_mesh']
		
		
		# Pickup the largest cluster
		clusters = sorted(nx.connected_components(multi_graph), key=len, reverse=True)
		lengths_clusters = [len(c) for c in clusters]
		nodes_cluster = clusters[lengths_clusters.index(lengths_clusters[0])] # 0: max
		g_largest_cluster = nx.MultiGraph( multi_graph.subgraph(nodes_cluster) )
		
		
		##
		## Get the ids of CaMKII binding beads outside the CaMKII condensate.
		##
		
		ids_CaMKII_binding_bead_grid_mesh    = d['condensate_CaMKII']['locs_in_grid_mesh'][ p.subunits['CaMKII binding site']['id'] ]
		cond_in_grid_mesh             = d['condensate_CaMKII'][type_condensate]
		ids_CaMKII_binding_bead_outside_cond = ids_CaMKII_binding_bead_grid_mesh[ cond_in_grid_mesh == 0 ] # 
		
		
		##
		## Get the ids of CaMKII binding beads outside the CaMKII condensate.
		##
		
		CaMKII_binding_sites_surface = get_properties_beads_CaMKII(g_largest_cluster, \
			ids_CaMKII_binding_bead_select = ids_CaMKII_binding_bead_outside_cond)
		
		mesh = utils.generate_mesh(cond_in_grid_mesh, flipz = False)
		t_centers = mesh.triangles_center
		f_normals = mesh.face_normals
		ids_closest_surface = [ np.argmin( np.linalg.norm(v['loc'] - t_centers,axis=1) ) for v in CaMKII_binding_sites_surface.values()]
		face_normals = [ f_normals[i] for i in ids_closest_surface ]
		angles_interface = [ angle(-fn, v['vector_to_hub']) for fn, v in zip( face_normals, CaMKII_binding_sites_surface.values() ) ]
		
		distance_to_hub_interface = [v['distance_to_hub'] for v in CaMKII_binding_sites_surface.values()]
		
		angles_interface = np.array( angles_interface )
		distance_to_hub_interface = np.array( distance_to_hub_interface )
		
		
		
		
		# Plot the polar scatter and save the figure
		#'''
		fig, ax = plot_polar_scatter( angles_interface, distance_to_hub_interface, prefix )
		fig.savefig( os.path.join(dir_imgs, '{}_polar_scatter.svg'.format( prefix ) ) )
		fig.savefig( os.path.join(dir_imgs, '{}_polar_scatter.png'.format( prefix ) ) , dpi=150)
		plt.show()
		#plt.clf()
		#plt.close(fig=fig)
		#'''
		streched_beads = distance_to_hub_interface >= length[prefix] - 1/np.sqrt(3)
		
		print('Total number of beads       : ', distance_to_hub_interface.shape[0] )
		print('Sigma total * cos(theta)    : ', np.sum( np.cos(angles_interface) ) )
		print('Average: ', np.sum( np.cos(angles_interface) ) / distance_to_hub_interface.shape[0] )
		print('Number of stretched beads   : ', np.sum( streched_beads ) )
		print('Sigma stretched * cos(theta): ', np.sum( streched_beads * np.cos(angles_interface) ) )
		print('Average: ',  np.sum( streched_beads * np.cos(angles_interface) ) / np.sum( streched_beads ) )
		
		
		
		# Plot the polar histogram and save the figure
		'''
		fig, ax = plot_polar_histogram(angles_interface, distance_to_hub_interface, prefix)
		fig.savefig( os.path.join(dir_imgs, '{}_polar_hist.svg'.format( prefix ) ) )
		fig.savefig( os.path.join(dir_imgs, '{}_polar_hist.png'.format( prefix ) ) , dpi=150)
		plt.show()
		plt.clf()
		plt.close(fig=fig)
		'''
		
		'''
		distance_to_hub_interface = np.average( distance_to_hub_interface )
		angle_interface = np.average( angle_interface )
		distance_to_hub_all = np.average( distance_to_hub_all )
		angle_all = np.average( angle_all )
		
		## For save
		angles_interface[ii] = angle_interface
		angles_all[ii] = angle_all
		
		ids_molecule_ids_CaMKII_beads_interfaces = [ v['id_molecule'] for id, v in CaMKII_binding_site.items() if id in ids_CaMKII_binding_bead_interface ]
		nums = [ids_molecule_ids_CaMKII_beads_interfaces.count(id) for id in list(set(ids_molecule_ids_CaMKII_beads_interfaces)) ]
		nums.sort()
		
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_title('Num beads participating surface region per holoenzyme, {}'.format(prefix) )
		ax.hist( nums, bins=np.arange(0,13), color = 'b')
		plt.show()
		print('num_CaMKII_molecules_participating_multiple (at interface) ', nums )
		
		
		print(prefix)
		print('num_CaMKII_binding_beads (at interface) ', len(ids_CaMKII_binding_bead_interface) )
		print('unconnected num_CaMKII_binding_beads (at interface)', num_free_interface)
		
		print('num_CaMKII_binding_beads in the total of cluster ', num_CaMKII_binding_beads_all )
		print('unconnected num_CaMKII_binding_beads  in the all of cluster ', num_free_all)
		
		
		#print('Ratio of Free   (surface): ', num_free_interface / len(ids_CaMKII_binding_bead_interface))
		#print('Ratio of Free       (all): ', num_free_all / num_CaMKII_binding_beads_all )
		
		
		print('Distance_to_hub (surface): ', distance_to_hub_interface )
		print('Distance_to_hub     (all): ', distance_to_hub_all)
		
		print('Angle (surface): ', angle_interface)
		print('Angle     (all): ', angle_all)
		'''
		
		
	# Save the edited data
	'''
	dd = {	'cluster_coefficients': cluster_coefficients,\
			'angles_interface': angles_interface, \
			'angles_all': angles_all}
	prefix = 'tension'
	suffix = 'cluster'
	utils.save(dir_edited_data, prefix, suffix, dd)
	'''
	
	
