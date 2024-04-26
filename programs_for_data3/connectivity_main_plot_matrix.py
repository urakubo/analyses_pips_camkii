
# https://stackoverflow.com/questions/59821151/plot-the-dendrogram-of-communities-found-by-networkx-girvan-newman-algorithm
# https://qiita.com/s-wakaba/items/a93f03f27137cff4a26c

import networkx as nx
from itertools import chain, combinations
import matplotlib.pyplot as plt
from matplotlib import patches

from scipy.cluster.hierarchy import dendrogram
import os, sys, pprint, itertools
import numpy as np

import utils
import parameters as p
import colormap as c

plt.rcParams.update(p.rc_param)

import centrality_main_calc_working as calc_central

	
def draw_adjacency_matrix(ax, G, node_order=None, partitions=None, color=None):
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
	adjacency_matrix = nx.to_numpy_array(G, dtype='bool', nodelist=node_order)

	#Plot adjacency matrix in toned-down black and white
	#fig = pyplot.figure(figsize=(5, 5)) # in inches
	ax.imshow(adjacency_matrix,
	              cmap="Greys",
	              interpolation="none")

	# The rest is just if you have sorted nodes by a partition and want to
	# highlight the module boundaries
	#ax = pyplot.gca()

	if partitions is not None:
		for p in partitions:
			p_start = p[0]
			p_end   = p[1]
			p_size  = p_end - p_start
			ax.add_patch(patches.Rectangle((p_start-0.5, p_start-0.5),
			                              p_size, # Width
			                              p_size, # Height
			                              facecolor=None,
			                              fill=False,
			                              edgecolor=color,
			                              linewidth=1))

	# girvan_newman1
def girvan_newman1(G, num_div = 4):
	def iterative_division(nodes, k):
		if k == 0 or len(nodes) < 3:
			return nodes
		else:
			comp = nx.community.girvan_newman(G.subgraph(nodes))
			segmented_nodes = [iterative_division(list(c), k-1) for c in next(comp)]
			segmented_nodes = sorted(segmented_nodes, key=len)
			return segmented_nodes
	
	nodes = list( G.nodes )
	rcm = iterative_division(nodes, k = num_div)
	print('rcm ', rcm)
	rcm_flat = utils.flatten(rcm)
	return rcm, rcm_flat
	
	

if __name__ == '__main__':

	## Init
	dir_target  =  'small_colony'
	prefix = '01_008'
	#prefix = '01_004'

	# lengths_clusters  [1503, 881, 699, 447, 274, 1, 1, 1, 1, 1]
	prefix = '01_004'
	nth_largest = 1

	# lengths_clusters  [845, 838, 793, 443, 372, 368, 1, 1, 1, 1]
	prefix = '00_004'
	nth_largest = 0 #0


	#'''
	dir_target  =  'CG_valency_length'
	prefix = '12_002'
	nth_largest = 0
	#'''

	# dir_edited_data  = 'valency_length'
	# prefix = '04_002'


	fig_title = '{}_{}'.format(nth_largest, prefix)
	print(fig_title)
	dir_edited_data  = os.path.join('data3', dir_target)
	dir_imgs         = os.path.join('imgs3', dir_target, 'connectivity_matrix_dendrogram')
	os.makedirs(dir_imgs, exist_ok=True)
	
	# Load graph.
	d = utils.load(dir_edited_data, prefix, 'connectivity_graph')
	
	# Make new graphs of CaMKII
	multi_graph_CaMKII, simple_graph_CaMKII, locs_hub, CaMKII_binding_site = \
		calc_central.make_new_graphs_CaMKII_connectivity(d, nth_largest =  nth_largest)
	G_ = multi_graph_CaMKII
	
	

	print('girvan_newman')
	# rcm, rcm_flat = girvan_newman1(G, num_div = 4)
	k = 3
	comp = nx.community.girvan_newman(G_)
	step = 0
	for communities in itertools.islice(comp, k):
		rcm = tuple(sorted(c) for c in communities)
		print('step ', step)
		step += 1
	print(rcm)
	
	a_rcm = list(itertools.accumulate([len(r) for r in rcm], initial=0))
	partitions = [[a_rcm[i], a_rcm[i+1]] for i in range(len(a_rcm)-1) ]
	rcm_flat = utils.flatten(rcm)
	#
	
	# Set color of dendrogram
	from cycler import cycler
	plt.rcParams["axes.prop_cycle"] = cycler( color=c.cols_list_ratio )
	
	# Plot and save figure
	fig = plt.figure(figsize=(4,4))
	
	ax2 = fig.add_axes([0.3,0.1,0.6,0.6])
	ax2.set_title(fig_title)
	draw_adjacency_matrix(ax2, G_, node_order = rcm_flat, partitions= partitions, color= 'r')
	#'''
	fig.savefig( os.path.join(dir_imgs, '{}_connect_matrix.svg'.format( fig_title ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}_connect_matrix.png'.format( fig_title ) ) , dpi=150)
	plt.show()
	#'''
	plt.close(fig=fig)
	
	
	
	data = {}
	data['multi_graph_CaMKII'] = multi_graph_CaMKII
	data['rcm']        = rcm
	data['rcm_flat']   = rcm_flat
	data['partitions']  = partitions
	utils.save(dir_edited_data, prefix, 'cluster', data)
	
