
# https://stackoverflow.com/questions/59821151/plot-the-dendrogram-of-communities-found-by-networkx-girvan-newman-algorithm
# https://qiita.com/s-wakaba/items/a93f03f27137cff4a26c

import os, sys, glob, pprint,itertools
import numpy as np
import itertools

from itertools import chain, combinations
import matplotlib.pyplot as plt

import networkx as nx
from scipy.cluster.hierarchy import dendrogram, set_link_color_palette
from networkx.algorithms.community.centrality import girvan_newman

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph

from specification_datasets import SpecDatasets


plt.rcParams.update(p.rc_param)
from cycler import cycler
plt.rcParams["axes.prop_cycle"] = cycler( color=c.cols_list_ratio )



def calc_dendrogram(communities):
	
	# building initial dict of node_id to each possible subset:
	node_id = 0
	init_node2community_dict = {node_id: communities[0][0].union(communities[0][1])}
	for comm in communities:
	    for subset in list(comm):
	        if subset not in init_node2community_dict.values():
	            node_id += 1
	            init_node2community_dict[node_id] = subset
	
	# turning this dictionary to the desired format in @mdml's answer
	node_id_to_children = {e: [] for e in init_node2community_dict.keys()}
	for node_id1, node_id2 in combinations(init_node2community_dict.keys(), 2):
	    for node_id_parent, group in init_node2community_dict.items():
	        if len(init_node2community_dict[node_id1].intersection(init_node2community_dict[node_id2])) == 0 and group == init_node2community_dict[node_id1].union(init_node2community_dict[node_id2]):
	            node_id_to_children[node_id_parent].append(node_id1)
	            node_id_to_children[node_id_parent].append(node_id2)

	# also recording node_labels dict for the correct label for dendrogram leaves
	node_labels = dict()
	for node_id, group in init_node2community_dict.items():
	    if len(group) == 1:
	        node_labels[node_id] = list(group)[0]
	    else:
	        node_labels[node_id] = ''

	# also needing a subset to rank dict to later know within all k-length merges which came first
	subset_rank_dict = dict()
	rank = 0
	for e in communities[::-1]:
	    for p in list(e):
	        if tuple(p) not in subset_rank_dict:
	            subset_rank_dict[tuple(sorted(p))] = rank
	            rank += 1
	subset_rank_dict[tuple(sorted(chain.from_iterable(communities[-1])))] = rank

	# ここを、切断するedgeの数の逆数にするなどして、浦久保に分かるようにする？？
	# my function to get a merge height so that it is unique (probably not that efficient)
	def get_merge_height(sub):
	    sub_tuple = tuple(sorted([node_labels[i] for i in sub]))
	    n = len(sub_tuple)
	    other_same_len_merges = {k: v for k, v in subset_rank_dict.items() if len(k) == n}
	    min_rank, max_rank = min(other_same_len_merges.values()), max(other_same_len_merges.values())
	    range = (max_rank-min_rank) if max_rank > min_rank else 1
	    return float(len(sub)) + 0.8 * (subset_rank_dict[sub_tuple] - min_rank) / range

	# finally using @mdml's magic, slightly modified:
	G           = nx.DiGraph(node_id_to_children)
	nodes       = G.nodes()
	leaves      = set( n for n in nodes if G.out_degree(n) == 0 )
	inner_nodes = [ n for n in nodes if G.out_degree(n) > 0 ]

	# Compute the size of each subtree
	subtree = dict( (n, [n]) for n in leaves )
	for u in inner_nodes:
	    children = set()
	    node_list = list(node_id_to_children[u])
	    while len(node_list) > 0:
	        v = node_list.pop(0)
	        children.add( v )
	        node_list += node_id_to_children[v]
	    subtree[u] = sorted(children & leaves)

	inner_nodes.sort(key=lambda n: len(subtree[n])) # <-- order inner nodes ascending by subtree size, root is last

	# Construct the linkage matrix
	leaves = sorted(leaves)
	index  = dict( (tuple([n]), i) for i, n in enumerate(leaves) )
	Z = []
	k = len(leaves)
	for i, n in enumerate(inner_nodes):
	    children = node_id_to_children[n]
	    x = children[0]
	    for y in children[1:]:
	        z = tuple(sorted(subtree[x] + subtree[y]))
	        i, j = index[tuple(sorted(subtree[x]))], index[tuple(sorted(subtree[y]))]
	        Z.append([i, j, get_merge_height(subtree[n]), len(z)]) # <-- float is required by the dendrogram function
	        index[z] = k
	        subtree[z] = list(z)
	        x = z
	        k += 1

	# dendrogram
	labels=[node_labels[node_id] for node_id in leaves]
	
	return Z, labels
	
	## 	print(init_node2community_dict)
	# {0: {0, 1, 2, 3, 4, 5, 6, ..., 431}, 
	# 1: {0, 2, 4, 5, 10, 11, 12, ... , 429, 430},
	# 2: {1, 3, 6, 7, 8, 9, 15,  ..., 426, 427, 428, 431},
	
	## Z
	# [[424, 425, 2.000090631018466, 2], [192, 193, 2.8000000000000003, 2], [104, 105, 2.7990936898153396, 2], ...,
	#  [858, 853, 194.0, 194], [860, 856, 297.0, 297], [859, 861, 432.0, 432]]
	
	## labels
	# [253, 430, 336, 11, 85, 360, 278, 372, 222, 128, ..., 381, 176, 380]
	
	#  [list(v)[0] for k, v, in init_node2community_dict.items() if len(v) == 1] # 元のnodesの並び順？
	#  [k for k, v, in init_node2community_dict.items() if len(v) == 1] 
	
	# pprint.pprint(node_labels) ここにシッカリ対応表がある。
	#  [node_labels[node_id] for node_id in leaves] <= この順に並べなおす。
	# leavesが末端のノードの様だ。
	# 
	
	# Save するべきは、node_id_to_children, subtree, node_labels, leaves, init_node2community_dict, G など？
	# もっと small scale で試してからの方が良い？



def get_clusters_from_dendrogram(R):

	node_order =  R['ivl']
	leaves_color_list = R['leaves_color_list']

	# Clustering for graph
	leaves_colors = list(set(leaves_color_list))
	ref = {col: i+1 for i, col in enumerate(leaves_colors)}
	ref['k'] = 0
	partition     = {n: ref[c] for n, c in zip(node_order, leaves_color_list)}
	
	maxid_color = max(partition.values()) + 1
	for k, v in partition.items():
		if v == 0:
			partition[k] = maxid_color
			maxid_color += 1
	#print('leaves_color_list ', leaves_color_list)
	#print('partition         ', partition)


	# Blocks for matrix
	blocks = []
	ref_col = ''
	ref_i   = 0
	for i, col in enumerate(leaves_color_list):
		if col == 'k' and ref_col != 'k':
			blocks.append([ref_i, i])
			ref_i   = i
			ref_col = col
		elif col != ref_col and ref_col != 'k':
			blocks.append([ref_i, i])
			ref_col = col
			ref_i   = i
		elif col != ref_col and ref_col == 'k':
			ref_col = col
			ref_i   = i
		else:
			pass
	if leaves_color_list[-1] != 'k':
		blocks.append([ref_i, len(leaves_color_list)])
	#print('blocks ', blocks)

	return node_order, blocks, partition, leaves_color_list



def plot_dendrogram(ax, Z, labels, G_, color_threshold):
	
	#set_link_color_palette(['k','b','r','y','c','g'])
	#from itertools import cycle
	#cols_list_uint = ['k','b','r','y','c','g']*10000
	# link_color_func=lambda x: 'black' link_color_func=lambda k: c.cols[k]
	
	R=dendrogram(Z, \
		labels=labels, \
		color_threshold = color_threshold, \
		above_threshold_color='k',\
		ax=ax)
	ax.set_yticks([0,200,400])
	ax.set_xticks([])
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	plt.setp(ax.collections,linewidth=1.0)
	
	'''
	# Show root
	xk = R['icoord'][-1]
	yk = R['dcoord'][-1]
	x = np.mean(xk[1:3])
	y1 = yk[1]
	y2 = yk[1] + 20 # y-coordinate of endpoint
	ax.plot([x, x], [y1, y2], color=R['color_list'][-1],linewidth=1.0)
	'''
	
	node_order, blocks, partition, leaves_color_list = get_clusters_from_dendrogram(R)
	print('len(leaves_color_list) ', len(leaves_color_list))
	leaves_color_list_squeezed = [k for k, g in itertools.groupby(leaves_color_list)]
	print()
	
	#x_labels = ax1.get_xmajorticklabels()
	
	return node_order, blocks, partition, leaves_color_list_squeezed



class PlotConnectivityMatrixDendrogram(SpecDatasets):
		
	def __init__( self ):
		
		pass
		
		
	def inspect( self ):
		
		for i, (f_lammpstrj, f_edited) in enumerate( zip(self.filenames_lammpstrj, self.filenames_edited) ):
			print('ID {}: {}, {}'.format(i, f_lammpstrj, f_edited))
		print()
		
		
	def plot_connectivity_save_graph( self, i ):
		
		filename_lammpstrj = self.filenames_lammpstrj[i]
		filename_edited    = self.filenames_edited[i]
		
		print('ID {}: {}, {}'.format(i, filename_lammpstrj, filename_edited))
		
		print('\n File setting')
		# Input
		prefix_load = filename_edited
		suffix_load = 'connectivity_graph'
		time_frame  = 0
		# Output
		self.savename_prefix     = '{}_{}'.format(time_frame, prefix_load)
		self.savename_suffix_img = 'connect_matrix_dendrogram'
		
		
		print('\n Load graph')
		d = utils.load(self.dir_edited_data, prefix_load, suffix_load )
		
		'''
		print('\n Make new graphs of CaMKII')
		nth_largest = 0
		multi_graph_CaMKII, _, _, _ = \
			utils_graph.make_new_graphs_CaMKII_connectivity(d, nth_largest =  nth_largest)
		
		
		print('\n Calc dendrogram')
		communities = list(nx.community.girvan_newman(multi_graph_CaMKII))
		Z, labels = calc_dendrogram( communities )
		
		
		print('\n Save intermidiate file')
		data = {}
		data['multi_graph_CaMKII'] = multi_graph_CaMKII
		data['communities']= communities
		data['Z']          = Z
		data['labels']     = labels
		
		utils.save(self.dir_edited_data, filename_edited, 'tmp_dendrogram', data)
		'''
		
		self.data_ = utils.load(self.dir_edited_data, filename_edited, 'tmp_dendrogram')
		multi_graph_CaMKII = self.data_['multi_graph_CaMKII']
		communities = self.data_['communities']
		Z      = self.data_['Z']
		labels = self.data_['labels']
		
		
		print('\n Plot figure and save it')
		
		color_threshold = 150
		fig = plt.figure(figsize=(4,4))
		ax1 = fig.add_axes([0.3,0.71,0.6,0.2])
		ax2 = fig.add_axes([0.3,0.1,0.6,0.6])
		
		node_order, blocks, partition, color_list = plot_dendrogram(ax1, Z, labels, multi_graph_CaMKII, color_threshold)
		utils_graph.draw_adjacency_matrix(ax2, multi_graph_CaMKII,  node_order = node_order, partitions= blocks, color= color_list)
		
		ax1.set_title(self.savename_prefix)
		self.show_save_img( fig, ax1, ax2 )
		
		# Decode color
		cols = plt.rcParams["axes.prop_cycle"].by_key()['color']
		print('cols ', cols)
		print('color_list ', color_list)
		color_list = [cols[int(id[1:])-1] for id in color_list]
		print('color_list ' , color_list)
		#utils.get_ratio_code(
		
		self.data = {}
		self.data['multi_graph_CaMKII'] = multi_graph_CaMKII
		self.data['Z']          = Z
		self.data['labels']     = labels
		self.data['node_order'] = node_order
		self.data['blocks']     = blocks
		self.data['partition']  = partition
		self.data['color_list'] = color_list
		
		utils.save(self.dir_edited_data, self.savename_prefix, 'cluster_dendrogram', self.data)
		
		
	def show_save_img( self, fig, ax1, ax2 ):
		
		dir_imgs = os.path.join(self.dir_imgs_root, 'connectivity_matrix_dendrogram')
		os.makedirs(dir_imgs, exist_ok=True)
		basename = self.savename_prefix + '_' + self.savename_suffix_img
		fig.savefig( os.path.join(dir_imgs, basename + '.svg' ) )
		fig.savefig( os.path.join(dir_imgs, basename + '.png' ), dpi=150 )
		plt.show()
		plt.clf()
		plt.close(fig=fig)
		
if __name__ == '__main__':
	
	pass
	