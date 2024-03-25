
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
	filenames = ['00_002']

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

		color_map = [c.cmap_universal_ratio[ v['species']] for n, v in g_largest_cluster.nodes.items()]

		nx.draw_networkx(	g_largest_cluster, \
							node_color=color_map, \
							with_labels=False, \
							node_size=20, \
							edge_color ='.4',\
							pos = nx.kamada_kawai_layout(g_largest_cluster))
							# nx.kamada_kawai_layout, nx.spring_layout(g_largest_cluster)
		plt.show()
		
		sys.exit(0)

		if type_centrality == 'betweenness':
		        centrality = nx.betweenness_centrality(g_largest_cluster)
		elif type_centrality == 'parcolation':
		        centrality = nx.percolation_centrality(g_largest_cluster)

		utils.save(dir_edited_data, prefix, type_centrality, centrality )

                
