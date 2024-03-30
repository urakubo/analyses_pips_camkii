
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
	
	rcm = list(nx.utils.reverse_cuthill_mckee_ordering(G))
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

	#'''
	type_centrality = 'betweenness' # 'parcolation', 'betweenness'

	# Shared part of initialization
	dir_edited_data  = os.path.join('data3', dir_edited_data)
	os.makedirs(dir_edited_data, exist_ok=True)
	
	#
	for prefix in filenames:
		print(prefix)		
		d = utils.load(dir_edited_data, prefix, 'connectivity_graph')
		
		g = d['simple_graph_CaMKII_GluN2B']
		g = d['multi_graph']
		
		clusters = sorted(nx.connected_components(g), key=len, reverse=True)
		lengths_clusters = [len(c) for c in clusters]
		print('lengths_clusters ', lengths_clusters[0:10])
		id_max = lengths_clusters.index(lengths_clusters[0])
		print('Picked-up cluster ', lengths_clusters[id_max] )
		
		
		max_cluster = clusters[id_max]
		g_largest_cluster = g.subgraph(max_cluster)
		g_largest_cluster = nx.Graph(g_largest_cluster)

		plot_histogram_centrality(g_largest_cluster, type_centrality)


		ids_GluN2B = [n for n, v in g_largest_cluster.nodes.items() if v['species'] == 'GluN2B' ]
		single_con_GluN2B = [id for id in ids_GluN2B if g_largest_cluster.degree[id] == 1 ]
		g_largest_cluster.remove_nodes_from( single_con_GluN2B )


		ids_GluN2B = [n for n, v in g_largest_cluster.nodes.items() if v['species'] == 'GluN2B' ]
		ids_CaMKII = [n for n, v in g_largest_cluster.nodes.items() if v['species'] == 'CaMKII' ]
		multi_graph = nx.MultiGraph()
		
		
		locs_hub = []
		for id, v in g_largest_cluster.nodes.items():
			if v['species'] == 'CaMKII':
				id_CaMKII_hub = np.nonzero(v['types_bead'] == p.subunits['CaMKII hub']['id'])
				
				position_CaMKII_hub = v['positions_grid_coord'][id_CaMKII_hub, :]
				position_CaMKII_hub = np.ravel( position_CaMKII_hub )
				#print('id ', id, ', position_CaMKII_hub ', position_CaMKII_hub)
				multi_graph.add_node(id, species = 'CaMKII', loc_hub = position_CaMKII_hub)
				locs_hub.append( position_CaMKII_hub )
		
		for id in ids_GluN2B:
			neighbors = [n for n in g_largest_cluster.neighbors(id)]
			if len(neighbors) == 2:
				multi_graph.add_edge(neighbors[0], neighbors[1], id = id)
			else:
				raise ValueError('Something wrong with neighbors', neighbors)
		
		
		locs_hub = np.array( locs_hub )
		barycenter_hub = np.average( locs_hub, axis = 0)
		for k, v in multi_graph.nodes.items():
			v['loc_hub'] -= barycenter_hub
		
		locs_hub -= barycenter_hub
		
		print('Num of GluN2B ', len(ids_GluN2B))
		print('Num of CaMKII ', len(ids_CaMKII))
		
		print('Total number of edges ', multi_graph.size())
		
		### Plot adjacency matrix
		#draw_adjacency_matrix(multi_graph)
		
		
		### Plot 3d matplotlib
		
		'''# Modularity
		parts = list( nx.community.greedy_modularity_communities(multi_graph) )
		values = {n: i for i, ns in enumerate(parts) for n in ns}
		n_color = np.asarray([values[n] for n in multi_graph.nodes()])
		'''
		n_color = np.asarray([v for k,v in multi_graph.degree]) # multi_graph.degree[id] for id in ids
		
		
		edge_xyz = np.array([(multi_graph.nodes[u]['loc_hub'], multi_graph.nodes[v]['loc_hub']) for u, v in multi_graph.edges()])
		fig = plt.figure()
		ax = fig.add_subplot(111, projection="3d")
		ax.scatter(*locs_hub.T, s=100, ec="w", c = n_color)
		for loc, num in zip( locs_hub, n_color ):
			ax.text(loc[0], loc[1], loc[2], str(num)) 
		
		#for vizedge in edge_xyz:
		#    ax.plot(*vizedge.T, color="tab:gray")
		fig.tight_layout()
		plt.show()
		
		
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_title(prefix)
		ax.hist( n_color , bins=np.arange(0,13))
		plt.show()
		
		###
		
		#draw_network_of_multi_graph(multi_graph)
		
		shortest_path = nx.average_shortest_path_length(multi_graph)
		print('Average shortest path of {} : {}'.format(prefix, shortest_path))
		
		communities = nx.community.greedy_modularity_communities(multi_graph)
		print('nx.average_clustering ', nx.average_clustering(g_largest_cluster) )
		
		sys.exit(0)
		
		utils.save(dir_edited_data, prefix, type_centrality, centrality )
		
 		
		
               
