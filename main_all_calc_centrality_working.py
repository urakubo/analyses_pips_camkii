
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
from networkx.algorithms.community.centrality import girvan_newman

def flatten(sequence):
    result = []

    for item in sequence:
        if isinstance(item, (list, tuple, range, dict, set, frozenset)):
            result.extend(flatten(item))
        else:
            result.append(item)

    return result
	
	
def plot_condensates_pyvista(pl, cond, color='green'): 
	
	flipz = False
	
	mesh = utils.generate_mesh(cond, flipz = flipz)
	
	
	# Add cube
	'''
	cube  = pyvista.Cube(center=(0,0,0), \
		x_length=utils.space[0], y_length=utils.space[1], z_length=utils.space[2])
	pl.add_mesh(cube, color='black', style='wireframe')
	'''
	
	pl.add_mesh(mesh, color=color, show_edges=False,  opacity=0.4)
	pl.set_background('white')
	
	
def draw_adjacency_matrix(G, node_order=None, partitions=[], colors=[]):
	"""
	- G is a netorkx graph
	- node_order (optional) is a list of nodes, where each node in G
	      appears exactly once
	- partitions is a list of node lists, where each node in G appears
	      in exactly one node list
	- colors is a list of strings indicating what color each
	      partition should be
	If partitions is specified, the same number of colors needs to be
	specified.
	"""
	# https://sociograph.blogspot.com/2012/11/visualizing-adjacency-matrices-in-python.html
	
	# comp = nx.community.girvan_newman(G)
	# for communities in itertools.islice(comp, k):
	
	
	#rcm = list(nx.utils.reverse_cuthill_mckee_ordering(G))
	
	def iterative_division(nodes, k):
		if k == 0 or len(nodes) < 3:
			return nodes
		else:
			comp = nx.community.girvan_newman(G.subgraph(nodes))
			segmented_nodes = [iterative_division(list(c), k-1) for c in next(comp)]
			segmented_nodes = sorted(segmented_nodes, key=len)
			return segmented_nodes
	nodes = list( G.nodes )
	rcm = iterative_division(nodes, k = 40)
	rcm = flatten(rcm)
	#print(rcm)
	
	adjacency_matrix = nx.to_numpy_array(G, dtype=np.bool, nodelist=rcm)

	#Plot adjacency matrix in toned-down black and white
	fig = pyplot.figure(figsize=(5, 5)) # in inches
	pyplot.imshow(adjacency_matrix,
	              cmap="Greys",
	              interpolation="none")

	# The rest is just if you have sorted nodes by a partition and want to
	# highlight the module boundaries
	assert len(partitions) == len(colors)
	ax = pyplot.gca()
	for partition, color in zip(partitions, colors):
	    current_idx = 0
	    for module in partition:
	        ax.add_patch(patches.Rectangle((current_idx, current_idx),
	                                      len(module), # Width
	                                      len(module), # Height
	                                      facecolor="none",
	                                      edgecolor=color,
	                                      linewidth="1"))
	        current_idx += len(module)
	plt.show()
	


def draw_network_sentrality(multi_graph, type_centrality):
	
	if type_centrality == 'betweenness':
	        centrality = nx.betweenness_centrality(g_largest_cluster)
	elif type_centrality == 'parcolation':
	        centrality = nx.percolation_centrality(g_largest_cluster)


	color_map = [c.cmap_universal_ratio[ v['species']] for n, v in g_largest_cluster.nodes.items()]
	#labels = {k: '{:.3f}'.format(v) for k, v in centrality.items() if v > 0.04}
	labels = {k: '{:.3f}'.format(v) for k, v in centrality.items() if g_largest_cluster.degree[k] > 2}

	degree = [g_largest_cluster.degree[k] for k, v in centrality.items() if v > 0.04]
	print('Numbers of partners of hub CaMKII ', degree)
	degree = [g_largest_cluster.degree[k] for k, v in centrality.items()]
	print('Numbers of partners of all CaMKII ', degree)

	nx.draw_networkx(	g_largest_cluster, \
						node_color=color_map, \
						with_labels=True, \
						labels=labels, \
						node_size=20, \
						edge_color ='.4',\
						pos = nx.kamada_kawai_layout(g_largest_cluster))
						# nx.kamada_kawai_layout, nx.spring_layout(g_largest_cluster)
	plt.show()
	
		
def draw_network_simple(multi_graph):
	nx.draw_networkx(	multi_graph, \
						with_labels=False, \
						node_size=20, \
						edge_color ='.4', pos = nx.kamada_kawai_layout(g_largest_cluster)) # pos = nx.kamada_kawai_layout(multi_graph)
	plt.suptitle(prefix)
	plt.show()
		
		
def draw_network_of_multi_graph(multi_graph):
	pos = nx.kamada_kawai_layout(multi_graph)
	
	
	parts = list( nx.community.greedy_modularity_communities(multi_graph) )
	values = {n: i for i, ns in enumerate(parts) for n in ns}
	n_color = np.asarray([values[n] for n in multi_graph.nodes()])
	
	
	nx.draw_networkx_nodes(multi_graph, pos, node_color = n_color, node_size = 100, alpha = 1)
	ax = plt.gca()
	for e in multi_graph.edges:
	    ax.annotate("",
	                xy=pos[e[0]], xycoords='data',
	                xytext=pos[e[1]], textcoords='data',
	                arrowprops=dict(arrowstyle="-", color="0.5",
	                                shrinkA=5, shrinkB=5,
	                                patchA=None, patchB=None,
	                                connectionstyle="arc3,rad=rrr".replace('rrr',str(0.3*e[2])
	                                ),
	                                ),
	                )
	plt.axis('off')
	plt.show()
	

def plot_histogram_centrality(multi_graph, type_centrality):
	
	if type_centrality == 'betweenness':
	        centrality = nx.betweenness_centrality(multi_graph)
	elif type_centrality == 'parcolation':
	        centrality = nx.percolation_centrality(multi_graph)
	
	
	fig  = plt.figure(figsize=(4, 4), tight_layout=True)
	ax = fig.add_subplot( 1,1,1 )
	ax.grid()
	ax.set_title(prefix)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_ylabel('(Number)')
	ax.set_ylabel('(Num of molecules)')
	ax.set_ylim([0.1, 1000])
	ax.hist(centrality.values(), range = (0.0001, 0.06), bins= 40, log=True )
	plt.show()

	
def make_new_graphs_CaMKII_connectivity(d, nth_largest = 0):

	# Use the multi graph.
	#g = d['simple_graph_CaMKII_GluN2B']
	g = d['multi_graph']
	
	# Pickup the largest cluster
	clusters = sorted(nx.connected_components(g), key=len, reverse=True)
	lengths_clusters = [len(c) for c in clusters]
	print('lengths_clusters ', lengths_clusters[0:10])
	id_max = lengths_clusters.index(lengths_clusters[nth_largest])
	print('Picked-up cluster ', lengths_clusters[id_max] )
	
	
	max_cluster = clusters[id_max]
	g_largest_cluster = nx.MultiGraph( g.subgraph(clusters[id_max]) )
	#g_largest_simple = nx.Graph( g_largest )
	#g_largest_multi  = nx.MultiGraph( g_largest )
	
	# Remove GluN2Bs that have single connections.
	ids_GluN2B = [n for n, v in g_largest_cluster.nodes.items() if v['species'] == 'GluN2B' ]
	single_con_GluN2B = [id for id in ids_GluN2B if g_largest_cluster.degree[id] == 1 ]
	g_largest_cluster.remove_nodes_from( single_con_GluN2B )
	
	# Check the node ids of GluN2B and CaMKII.
	ids_GluN2B = [n for n, v in g_largest_cluster.nodes.items() if v['species'] == 'GluN2B' ]
	ids_CaMKII = [n for n, v in g_largest_cluster.nodes.items() if v['species'] == 'CaMKII' ]
	
	# Set the new graphs of CaMKIIs
	multi_graph_CaMKII  = nx.MultiGraph()
	simple_graph_CaMKII = nx.Graph()
	# Location of hubs
	locs_hub = []
	# Location info of CaMKII interaction beads
	CaMKII_binding_site = {}
	
	for id, v in g_largest_cluster.nodes.items():
		if v['species'] == 'CaMKII':
			# Make CaMKII hub beads in the new graphs
			id_CaMKII_hub = np.nonzero(v['types_bead'] == p.subunits['CaMKII hub']['id'])
			position_CaMKII_hub = np.ravel( v['positions_grid_coord'][id_CaMKII_hub, :] )
			id_bead       = v['ids_bead'][0][id_CaMKII_hub[0][0]]
			# print('id_CaMKII_hub[0][0] ', id_CaMKII_hub[0][0], ', id_bead ', id_bead)
			multi_graph_CaMKII.add_node(id, id_bead = id_bead, loc_hub = position_CaMKII_hub)
			simple_graph_CaMKII.add_node(id, id_bead = id_bead, loc_hub = position_CaMKII_hub)
			locs_hub.append( position_CaMKII_hub )
			
			# Make CaMKII interaction beads
			id_CaMKII_binding = np.nonzero(v['types_bead'] == p.subunits['CaMKII binding site']['id'])
			position_CaMKII_binding = v['positions_grid_coord'][id_CaMKII_binding, :][0]
			ids_bead       = v['ids_bead'][0][id_CaMKII_binding]
			#print('position_CaMKII_binding: ', position_CaMKII_binding, ', ids_bead: ',ids_bead)
			for id_bead_, loc in zip(ids_bead, position_CaMKII_binding):
				CaMKII_binding_site[id_bead_] = {}
				CaMKII_binding_site[id_bead_]['loc'] = loc
				CaMKII_binding_site[id_bead_]['connect'] = 0
				#print('id2: ', id2, ' loc: ', loc)
	
	# Centering of the hub locations.
	locs_hub = np.array( locs_hub )
	barycenter_hub = np.average( locs_hub, axis = 0)
	for k, v in multi_graph_CaMKII.nodes.items():
		v['loc_hub'] -= barycenter_hub
	locs_hub -= barycenter_hub
	
	
	# Make the connections between CaMKIIs through GluN2B.
	for id in ids_GluN2B:
		neighbors = [n for n in g_largest_cluster.neighbors(id)]
		#print('neighbors ', neighbors)
		#Make the graph connection. The GluN2B may bind to PSD95.
		if (len(neighbors) == 2) and (neighbors[0] in ids_CaMKII) and (neighbors[1] in ids_CaMKII):
			multi_graph_CaMKII.add_edge(  neighbors[0], neighbors[1], id = id)
			simple_graph_CaMKII.add_edge( neighbors[0], neighbors[1], id = id)
		# Loop connection may exist (from-and-to one CaKII)
		#else:
		#	raise ValueError('Something wrong with neighbors', neighbors)
		
		# Get the edges for GluN2B to get the beads of CaMKII.
		if len(neighbors) != 0:
			connected_edges = g_largest_cluster.edges(id, keys=True)
			#print('connected_edges ', connected_edges)
			for connected_edge in connected_edges:
				id_bead1_connected = g_largest_cluster.edges[connected_edge]['id_bead1']
				id_bead2_connected = g_largest_cluster.edges[connected_edge]['id_bead2']
				
				# The other side of GluN2B may bind to PSD95. Such connection was rejected.
				# Loop connections are considered to be connected.
				if id_bead1_connected in CaMKII_binding_site.keys():
					CaMKII_binding_site[id_bead1_connected]['connect'] +=1
				if id_bead2_connected in CaMKII_binding_site.keys():
					CaMKII_binding_site[id_bead2_connected]['connect'] +=1
				#print('connected bead1: ',  g_largest_cluster.edges[connected_edge]['id_bead1'] )
				#print('connected bead2: ',  g_largest_cluster.edges[connected_edge]['id_bead2'] )
				
				# CaMKII_binding_site[neighbors[0]]['connect'] +=1
				# CaMKII_binding_site[neighbors[1]]['connect'] +=1
	
	
	print('Num of GluN2B (edges): ', len(ids_GluN2B))
	print('Num of CaMKII (nodes): ', len(ids_CaMKII))
	
	return multi_graph_CaMKII, simple_graph_CaMKII, locs_hub, CaMKII_binding_site



def plot_3D_pvista_CaMKII_interface_region(d):
	
	pl = pyvista.Plotter(window_size=[400,1000], shape=(3, 1), border=False)
	pl.subplot(0, 0)
	cond = d['condensate_CaMKII']['condensate_CaMKII_in_grid_mesh']
	plot_condensates_pyvista(pl, cond, color='green')
	cond = d['condensate_CaMKII']['region_interface']
	plot_condensates_pyvista(pl, cond, color='red')
	pl.add_mesh(utils.square_yz(), color='black', style='wireframe')
	pl.view_yz()
	pl.camera.roll -= 90
	pl.camera.Zoom(1)
	pl.show(interactive=True, auto_close=True) # off_screen = True
	#filename = os.path.join(dir_imgs, '{}_{}.png'.format(prefix, suffix))
	#pl.screenshot(filename)
	
	
def plot_3d_binding_beads_matplotlib(prefix, CaMKII_binding_site, ids_CaMKII_binding_interface):
	
	fig = plt.figure()
	ax = fig.add_subplot(111, projection="3d")
	ax.set_title(prefix)
	for id, v in CaMKII_binding_site.items():
		if id in ids_CaMKII_binding_interface:
			if v['connect'] > 0:
				ax.scatter(v['loc'][0],v['loc'][1],v['loc'][2], s=100, ec="w", color= 'k')
			else:
				ax.scatter(v['loc'][0],v['loc'][1],v['loc'][2], s=100, ec="w", color= 'red')
		# ax.text(loc[0], loc[1], loc[2], str(num), color = col)
	fig.tight_layout()
	plt.show()
	
	
def plot_3d_hubs_matplotlib(d, prefix, multi_graph_CaMKII, locs_hub):
	
	ids_CaMKII_hub_bead_interface = d['condensate_CaMKII']['interface_beads'][p.subunits['CaMKII hub']['id']].tolist()
	# Plot 3D
	nums_neighbors = np.asarray([v for k, v in multi_graph_CaMKII.degree]) # multi_graph.degree[id] for id in ids
	ids_bead = [v['id_bead'] for k, v in multi_graph_CaMKII.nodes.items()]
	
	fig = plt.figure()
	ax = fig.add_subplot(111, projection="3d")
	ax.set_title( prefix )
	ax.scatter(*locs_hub.T, s=100, ec="w", c = nums_neighbors)
	num_partners_interface_CaMKII_hub = []
	for loc, num, id_bead in zip( locs_hub, nums_neighbors, ids_bead ):
		if id_bead in ids_CaMKII_hub_bead_interface:
			col = 'red'
			num_partners_interface_CaMKII_hub.append(num)
		else:
			col = 'k'
		ax.text(loc[0], loc[1], loc[2], str(num), color = col)
	# edge_xyz = np.array([(multi_graph.nodes[u]['loc_hub'], multi_graph.nodes[v]['loc_hub']) for u, v in multi_graph.edges()])
	#for vizedge in edge_xyz:
	#    ax.plot(*vizedge.T, color="tab:gray")
	
	print('num_partners_interface_CaMKII_hub ', num_partners_interface_CaMKII_hub)
	fig.tight_layout()
	plt.show()
	
	
def plot_hist_num_connections(prefix, multi_graph_CaMKII, num_partners_interface_CaMKII_hub = None):
	nums_neighbors = np.asarray([v for k, v in multi_graph_CaMKII.degree])
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_title('num_partners_interface_CaMKII_hub, {}'.format(prefix) )
	ax.hist( nums_neighbors, bins=np.arange(0,13), color = 'b')
	if num_partners_interface_CaMKII_hub is not None:
		ax.hist( num_partners_interface_CaMKII_hub, bins=np.arange(0,13), color = 'r')
	plt.show()
	
	
if __name__ == '__main__':
	
	# Dataset 2:  Valency length
	'''
	filenames = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
	dir_edited_data  = 'valency_length'
	filenames = ['12_004'] # '12_002', '12_004', '12_005', '12_006'
	'''
        
	# Conc dependence
	'''
	dir_edited_data  =  'conc_dependence'
	filenames 	 = [str(i).zfill(3) for i in range(48) ] # 70
	filenames 	 = [str(9).zfill(3) ] # 70
	filenames 	 = [str(11).zfill(3) ] # 70
	'''	
	
	
	#'''
	dir_edited_data  =  'small_colony'
	filenames = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(3) for id_f in range(10) ]
	
	filenames = ['01_004'] # 00: len-3, 01: len-9, 02: linear
	filenames = ['00_004'] # 00: len-3, 01: len-9, 02: linear
	#filenames = ['02_000'] # 00: len-3, 01: len-9, 02: linear
	#'''
	
	type_centrality = 'betweenness' # 'parcolation', 'betweenness'
	
	# Shared part of initialization
	dir_edited_data  = os.path.join('data3', dir_edited_data)
	os.makedirs(dir_edited_data, exist_ok=True)
	
	#
	for prefix in filenames:
		print(prefix)		
		d = utils.load(dir_edited_data, prefix, 'connectivity_graph')
		
		#plot_3D_pvista_CaMKII_interface_region(d)
		
		# Make new graphs of CaMKII
		multi_graph_CaMKII, simple_graph_CaMKII, locs_hub, CaMKII_binding_site = \
			make_new_graphs_CaMKII_connectivity(d, nth_largest = 0)
		
		
		# plot_histogram_centrality(g_largest_cluster, type_centrality)
		
		### Plot adjacency matrix
		#draw_adjacency_matrix(multi_graph_CaMKII)
		
		
		### Count the number of unbinding beads at the interface region.
		ids_CaMKII_binding_bead_interface = d['condensate_CaMKII']['interface_beads'][p.subunits['CaMKII binding site']['id']].tolist()
		num_free = [ (id in ids_CaMKII_binding_bead_interface) and v['connect'] == 0 for id, v in CaMKII_binding_site.items()]
		num_free = sum( num_free )
		
		num_CaMKII_binding_interface = len(set(ids_CaMKII_binding_bead_interface) & set(CaMKII_binding_site.keys()))
		print(prefix)
		print('num_CaMKII_binding_beads (at interface) ', num_CaMKII_binding_interface)
		print('unconnected num_CaMKII_binding_beads (at interface)', num_free)
		print('Ratio of Free : ', num_free / num_CaMKII_binding_interface)
		
		num_CaMKII_binding_beads_total = len(CaMKII_binding_site)
		num_free = sum( [v['connect'] == 0 for v in CaMKII_binding_site.values() ] )
		print('num_CaMKII_binding_beads in the total of cluster ', num_CaMKII_binding_beads_total )
		print('unconnected num_CaMKII_binding_beads  in the total of cluster ', num_free)
		print('Ratio of Free : ', num_free / num_CaMKII_binding_beads_total )
		
		
		### Plot 3d binding beads in matplotlib
		
		#plot_3d_binding_beads_matplotlib(prefix, CaMKII_binding_site, ids_CaMKII_binding_bead_interface)
		
		
		### Plot 3d hubs in matplotlib
		
		#plot_3d_hubs_matplotlib(d, prefix, multi_graph_CaMKII, locs_hub)
		
		
		plot_hist_num_connections(prefix, multi_graph_CaMKII, num_partners_interface_CaMKII_hub = None)
		
		
		# draw_network_of_multi_graph(multi_graph)
		
		# shortest_path = nx.average_shortest_path_length(multi_graph)
		#print('Average shortest path of {} : {}'.format(prefix, shortest_path))
		
		#communities = nx.community.greedy_modularity_communities(multi_graph)
		'''
		print('average_clustering ', nx.average_clustering(simple_graph_CaMKII) )
		print('average_clustering ', nx.average_clustering(multi_graph_CaMKII) )
		
		'''
		
		
		
		# utils.save(dir_edited_data, prefix, type_centrality, centrality )
		
 		
 		# Obsolete
		#print('transitivity       ', nx.transitivity(simple_graph) )
		#print('average_shortest_path_length ', nx.average_shortest_path_length(simple_graph) )
		
		## Modularity
		#parts = list( nx.community.greedy_modularity_communities(multi_graph) )
		#values = {n: i for i, ns in enumerate(parts) for n in ns}
		#n_color = np.asarray([values[n] for n in multi_graph.nodes()])
