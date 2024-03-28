
import os, sys, glob, pickle, pprint
import numpy as np
import math

import matplotlib
import matplotlib.pyplot as plt

import networkx as nx


import utils
import parameters as p
import colormap as c
	
	
	#
	# grid_coord: Molecular coordinate in the grid space (0, 1,..., 119) 
	# (numpy uint) [(x0,y0,z0), (x1,y1,z1), ..., (xn,yn,zn)]
	#
	# real_coord: Molecular coordinate in the real space [-60, 60) 
	# (numpy float) [(x0,y0,z0), (x1,y1,z1), ..., (xn,yn,zn)]
	#
	# grid_mesh : molecular existance in the 3D space
	# (numpy bool in 3D space) (p.space[0], p.space[1], p.space[2])
	#
	
	
	
if __name__ == '__main__':
	
	# Dataset 2:  Valency length
	'''
	filenames = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
	dir_edited_data  = 'valency_length'
	'''
        
	# Conc dependence
	#'''
	dir_edited_data  =  'conc_dependence'
	filenames 	 = [str(i).zfill(3) for i in range(48) ] # 70
	filenames 	 = [str(9).zfill(3) ] # 70
	#'''	
	

	dir_edited_data  =  'small_colony'
	filenames = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(3) for id_f in range(10) ]
	
	#filenames = ['01_003'] # 00: len-3, 01: len-9, 02: linear
	filenames = ['00_004'] # 00: len-3, 01: len-9, 02: linear

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
		# g = d['multi_graph']
		
		clusters = list(nx.connected_components(g))
		lengths_clusters = [len(c) for c in clusters]
		print('max(lengths_clusters) ', max(lengths_clusters))
		id_max = lengths_clusters.index(max(lengths_clusters))
		max_cluster = clusters[id_max]
		g_largest_cluster = g.subgraph(max_cluster)
		g_largest_cluster = nx.Graph(g_largest_cluster)

		ids_GluN2B = [n for n, v in g_largest_cluster.nodes.items() if v['species'] == 'GluN2B' ]
		single_con_GluN2B = [id for id in ids_GluN2B if g_largest_cluster.degree[id] == 1 ]
		g_largest_cluster.remove_nodes_from( single_con_GluN2B )


		ids_GluN2B = [n for n, v in g_largest_cluster.nodes.items() if v['species'] == 'GluN2B' ]
		ids_CaMKII = [n for n, v in g_largest_cluster.nodes.items() if v['species'] == 'CaMKII' ]
		print('Num of GluN2B ', len(ids_GluN2B))
		print('Num of CaMKII ', len(ids_CaMKII))
		
		multi_graph = nx.MultiGraph()
		multi_graph.add_nodes_from(ids_CaMKII, species = 'CaMKII')
		for id in ids_GluN2B:
			neighbors = [n for n in g_largest_cluster.neighbors(id)]
			if len(neighbors) == 2:
				multi_graph.add_edge(neighbors[0], neighbors[1], id = id)
			else:
				raise ValueError('Something wrong with neighbors', neighbors)
		
		
		nx.draw_networkx(	multi_graph, \
							with_labels=False, \
							node_size=20, \
							edge_color ='.4', pos = nx.kamada_kawai_layout(g_largest_cluster)) # pos = nx.kamada_kawai_layout(multi_graph)
		plt.suptitle(prefix)
		plt.show()
		
		'''
		if type_centrality == 'betweenness':
		        centrality = nx.betweenness_centrality(g_largest_cluster)
		elif type_centrality == 'parcolation':
		        centrality = nx.percolation_centrality(g_largest_cluster)


		color_map = [c.cmap_universal_ratio[ v['species']] for n, v in g_largest_cluster.nodes.items()]

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
		
		'''
		
		
		centrality = nx.betweenness_centrality(multi_graph)
		
		
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
		
		
		shortest_path = nx.average_shortest_path_length(multi_graph)
		print('Average shortest path of {} : {}'.format(prefix, shortest_path))
		
		#communities = nx.community.greedy_modularity_communities(multi_graph)
		print('nx.average_clustering(G) ', nx.average_clustering(g_largest_cluster) )
		
		sys.exit(0)
		
		utils.save(dir_edited_data, prefix, type_centrality, centrality )
		
 		
		
               
