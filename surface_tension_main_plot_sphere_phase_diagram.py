
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
	
	
def arrange_graph_bar(ax, panel_size_x, panel_size_y):
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax_h, ax_w = ax.bbox.height, ax.bbox.width
	loc = ax.get_position()
	ax.set_position([loc.x0, loc.y0, panel_size_x, panel_size_y])
	ax.set_xlabel('Linker length')
	
def plot_graphs(ave_cos, pull_force, pull_force_per_area):
	
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
	ax.set_title('Cosine similarity')
	ax.bar(*zip(*ave_cos.items()), width=0.6 , color='gray' ) # , color=colormap_conc_bargraph
	ax.set_ylim([0,1.0])
	ax.set_ylabel('E[cos(theta)]')
	arrange_graph_bar(ax, panel_size/4, panel_size)
	
	column = 5
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*1 )
	ax.set_title('Contraction force')
	ax.bar(*zip(*pull_force.items()), width=0.6, color='gray' ) # , color=colormap_conc_bargraph
	arrange_graph_bar(ax, panel_size/4, panel_size)
	
	column = 7
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*1 )
	ax.set_title('Contraction force per area')
	ax.bar(*zip(*pull_force_per_area.items()), width=0.6, color='gray' ) # , color=colormap_conc_bargraph
	arrange_graph_bar(ax, panel_size/4, panel_size)
	
	return fig
	
	
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
	
	fig = plt.figure(figsize=(2, 2))
	ax = fig.add_subplot(111, projection='polar')
	ax.scatter(angle_interface, distance_to_hub_interface, color='k', s=0.5)
	ax.set_title('Polar coordinates',fontsize=18)
	
	ax.set_title('{}'.format(prefix) )
	ax.set_theta_offset(np.pi / 2.0 * 3)
	ax.set_rlabel_position(180-20)
	
	return fig, ax
	
	
def plot_polar_histogram(angle_interface, distance_to_hub_interface, prefix, length):
	
	rbins = np.arange(0, length+1, 1)
	abins = np.linspace(0,2*np.pi, 30)
	abins_s = np.hstack( [abins[1:], abins[0]] )
	# https://en.wikipedia.org/wiki/Spherical_segment
	normal_abins = 2 * np.pi * np.abs( np.cos(abins[1:]) - np.cos(abins[:-1]) )
	normal_abins = normal_abins.reshape(normal_abins.shape[0],-1)
	normal_rbins = 1/3 * ( rbins[1:]**3 - rbins[:-1]**3 )
	
	hist, _, _ = np.histogram2d(angle_interface, distance_to_hub_interface, bins=(abins, rbins))
	
	hist = hist / normal_abins / normal_rbins / angle_interface.shape[0]
	hist = hist + np.flipud(hist)
	A, R = np.meshgrid(abins, rbins)
	fig, ax = plt.subplots(figsize=(2, 2), subplot_kw=dict(projection="polar"))
	pc = ax.pcolormesh(A, R, hist.T, cmap="Greys", vmin=0)  # 0.025 0.0015
	ax.set_rlim(0, length)
	ax.set_title('{}'.format(prefix) )
	#fig.colorbar(pc)
	ax.set_theta_offset(np.pi / 2.0 * 3)
	ax.set_rlabel_position(180-20)
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
			
			# print('position_CaMKII_binding: ', position_CaMKII_binding)
			# print('position_CaMKII_hub    : ',position_CaMKII_hub)
			for id_bead_, loc in zip(ids_bead, position_CaMKII_binding):
				vector_to_hub = position_CaMKII_hub - loc
				CaMKII_binding_site[id_bead_] = {}
				CaMKII_binding_site[id_bead_]['loc'] = loc
				CaMKII_binding_site[id_bead_]['connect'] = 0
				CaMKII_binding_site[id_bead_]['id_molecule'] = id
				CaMKII_binding_site[id_bead_]['vector_to_hub'] = vector_to_hub
				CaMKII_binding_site[id_bead_]['distance_to_hub'] = np.linalg.norm( vector_to_hub ) / np.sqrt(3)
				
				CaMKII_binding_site[id_bead_]['angle_center_of_mass_hub'] = angle(-loc, vector_to_hub)
				
				#CaMKII_binding_site[id_bead_]['angle'] = angle(vector_to_hub, -loc)
	
	return CaMKII_binding_site
	
	
class MatrixValencyLength():
	
	def __init__( self ):
		# Parameters
		self.sigma = 2
		dir_target = 'CG_valency_length'
		self.dir_edited_data = os.path.join('data3', dir_target)
		self.dir_imgs        = os.path.join('imgs3', dir_target,'phase_diagram')
		os.makedirs(self.dir_imgs, exist_ok=True)
		
		self.num_rows		= len( p.valencies )
		self.num_columns	= len( p.lengths )
		
		self.radius = utils.load(self.dir_edited_data, 'Radiuses', 'CaMKII_bead')
		
		
		self.ave_cos             = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		self.pull_force          = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		self.pull_force_per_area = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		
		
	def plot_polar_histogram( angles, distance_to_hub, prefix, length ):
		fig, ax = plot_polar_histogram( angles, distance_to_hub, prefix, length )
		fig.savefig( os.path.join(dir_imgs, '{}_polar_hist.svg'.format( prefix ) ) )
		fig.savefig( os.path.join(dir_imgs, '{}_polar_hist.png'.format( prefix ) ) , dpi=150)
		plt.show()
		#plt.clf()
		#plt.close(fig=fig)
		
		
	def plot_polar_scatter( angles, distance_to_hub, prefix ):
		fig, ax = plot_polar_scatter( angles, distance_to_hub, prefix )
		fig.savefig( os.path.join(dir_imgs, '{}_polar_scatter.svg'.format( prefix ) ) )
		fig.savefig( os.path.join(dir_imgs, '{}_polar_scatter.png'.format( prefix ) ) , dpi=150)
		plt.show()
		#plt.clf()
		#plt.close(fig=fig)
		
	def prepare_plot( self, title ):
		# self.fig  = plt.figure(figsize=(5, 5))
		self.fig  = plt.figure(figsize=(3.5, 3.5))
		self.fig.subplots_adjust(wspace=0.4,  hspace=0.6)
		ax = self.fig.add_subplot( 1, 1, 1 )
		ax.set_title(title)
		ax.set_xlabel('Linker length (l.u.)')
		ax.set_ylabel('Valency')
		return ax
		
	def run_calc( self ):
		
		for i, valency in enumerate(p.valencies):
			for j, length in enumerate(p.lengths):
				
				# Load data
				
				prefix      = p.fnames_valency[valency]+'_'+p.fnames_length[length]
				print('Target file: ', prefix)
				d           = utils.load(self.dir_edited_data, prefix, 'connectivity_graph')
				multi_graph = d['multi_graph']
				cond_CaMKII = d['condensate_CaMKII']['condensate_CaMKII_in_grid_mesh']
				
				
				# Pickup the largest cluster
				
				clusters = sorted(nx.connected_components(multi_graph), key=len, reverse=True)
				lengths_clusters = [len(c) for c in clusters]
				nodes_cluster = clusters[lengths_clusters.index(lengths_clusters[0])] # 0: max
				g_largest_cluster = nx.MultiGraph( multi_graph.subgraph(nodes_cluster) )
				
				
				## Get the ids of CaMKII binding beads outside the CaMKII condensate.
				
				ids_CaMKII_binding_bead_grid_mesh = \
					d['condensate_CaMKII']['locs_in_grid_mesh'][ p.subunits['CaMKII binding site']['id'] ]
				
				
				## Get the ids of CaMKII binding beads outside the CaMKII condensate.
				
				CaMKII_binding_sites = get_properties_beads_CaMKII(g_largest_cluster)
				
				angles = [ v['angle_center_of_mass_hub'] for v in CaMKII_binding_sites.values() ]
				distance_to_hub = [v['distance_to_hub'] for v in CaMKII_binding_sites.values()]
				
				angles = np.array( angles )
				distance_to_hub = np.array( distance_to_hub )
				
				
				# plot_polar_scatter( angles, distance_to_hub, prefix )
				# plot_polar_histogram( angles, distance_to_hub, prefix, length )
				
				streched_linkers1 = distance_to_hub >= length - 1 # - 1 -1/np.sqrt(3)##############
				stretched_linkers_cos = np.sum( streched_linkers1 * np.cos(angles) ) 
				ave_stretched_linkers_cos = stretched_linkers_cos / np.sum( streched_linkers1 )
				
				#'''
				print('Total number of beads       : ', distance_to_hub.shape[0] )
				print('Sigma total * cos(theta)    : ', np.sum( np.cos(angles) ) )
				print('Average: ', np.sum( np.cos(angles) ) / distance_to_hub.shape[0] )
				
				print('Number of stretched_linkers   : ', np.sum( streched_linkers1 ) )
				print('Sigma stretched_linkers * cos(theta): ', stretched_linkers_cos )
				print('Average: ',  ave_stretched_linkers_cos )
				#'''
				
				self.ave_cos[j,i]    = np.sum( np.cos(angles) ) / distance_to_hub.shape[0]
				self.pull_force[j,i] = stretched_linkers_cos
				self.pull_force_per_area[j,i] = stretched_linkers_cos / (4*np.pi*self.radius[prefix]*self.radius[prefix])
				
				
	def save_data( self ):
			
			d = {}
			d['ave_cos']             = self.ave_cos
			d['pull_force']          = self.pull_force
			d['pull_force_per_area'] = self.pull_force_per_area
			
			prefix = 'Surface_tension'
			suffix = 'matrix'
			utils.save(self.dir_edited_data, prefix, suffix, d)
			
	def load_data( self ):
			
			prefix = 'Surface_tension'
			suffix = 'matrix'
			d = utils.load(self.dir_edited_data, prefix, suffix)
			
			self.ave_cos             = d['ave_cos']
			self.pull_force          = d['pull_force']
			self.pull_force_per_area = d['pull_force_per_area']
			
			
	def plot_figures( self ):
		
		title    = 'Cos similarity'
		filename = 'Cos_similarity'
		colormap =  c.cmap_white_green_universal
		levels   = np.linspace(0,1,10)
		ax = self.prepare_plot(title)
		cs, cb = utils.plot_a_panel(ax, self.ave_cos, p.lengths, p.valencies, colormap, levels)
		self.save_a_fig( filename )
		
		title    = 'Contraction force'
		filename = 'Contraction_force'
		colormap = c.cmap_white_red_universal
		levels   = np.linspace(0,3000,10)
		ax = self.prepare_plot(title)
		cs, cb = utils.plot_a_panel(ax, self.pull_force, p.lengths, p.valencies, colormap, levels)
		self.save_a_fig( filename )
		
		title    = 'Contraction force per area'
		filename = 'Contraction_force_per_area'
		colormap = plt.colormaps['Greys']
		levels   = np.linspace(0,0.6,10)
		ax = self.prepare_plot(title)
		cs, cb = utils.plot_a_panel(ax, self.pull_force_per_area, p.lengths, p.valencies, colormap, levels)
		self.save_a_fig( filename )
		
		
	def save_a_fig( self, filename ):
		
		self.fig.savefig( os.path.join(self.dir_imgs, filename + '.svg' ) )
		self.fig.savefig( os.path.join(self.dir_imgs, filename + '.png' ) , dpi=150)
		plt.show()
		#plt.clf()
		#plt.close(fig=self.fig)
		
		
if __name__ == '__main__':
	
	graph = MatrixValencyLength()
	#graph.run_calc()
	#graph.save_data()
	graph.load_data()
	graph.plot_figures()
	
