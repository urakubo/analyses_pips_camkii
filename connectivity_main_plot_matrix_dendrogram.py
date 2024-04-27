
# https://stackoverflow.com/questions/59821151/plot-the-dendrogram-of-communities-found-by-networkx-girvan-newman-algorithm
# https://qiita.com/s-wakaba/items/a93f03f27137cff4a26c

import os, sys, glob, pprint,itertools
import numpy as np


from itertools import chain, combinations
import matplotlib.pyplot as plt

import networkx as nx
from scipy.cluster.hierarchy import dendrogram
from networkx.algorithms.community.centrality import girvan_newman

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph

plt.rcParams.update(p.rc_param)



def calc_dendrogram(G_):
	communities = list(nx.community.girvan_newman(G_))

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



def draw_network_simple(ax, G, pos, cols):

	# pos = nx.kamada_kawai_layout(G)
	
	# https://stackoverflow.com/questions/43541376/how-to-draw-communities-with-networkx
	# https://qiita.com/yuru/items/4208e13a399773c2b9f8
	
	#draw(G, pos=None, ax=None, **kwds)
	
	nx.draw_networkx_nodes(	G, \
				ax = ax, \
				pos = pos, \
				node_size=25, \
				node_color = cols, \
				edgecolors = 'k', \
				linewidths = 1.0)
	
	nx.draw_networkx_edges(	G, \
				pos = pos, \
				ax = ax, \
				edge_color =c.cmap_universal_ratio['CaMKII'],\
				width = 0.5)
	
	ax.set_aspect('equal')
	
	xy=np.array([v.tolist() for v in pos.values()] )
	x_min_max = [np.min(xy[:,0]), np.max(xy[:,0])/10 + np.min(xy[:,0])* 9/10 ]
	ymax = np.max(xy[:,1])
	
	ax.plot(x_min_max, [ymax, ymax], 'k-')
	ax.plot(x_min_max, [ymax, ymax], 'k-')
	
	ax.set_title('{:.3g}'.format( ( np.max(xy[:,0])-np.min(xy[:,0]))/10 ) )
	plt.box(False)
	#plt.suptitle(prefix)
	return


if __name__ == '__main__':


	# Input file
	i = 6
	prefix      = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(3) for id_f in range(7)][i]
	dir_target  = 'small_colony'
	nth_largest = 0


	# Shared init
	fig_title = '{}_{}'.format(nth_largest, prefix)
	print(fig_title)
	dir_edited_data  = os.path.join('data4', dir_target)
	dir_imgs         = os.path.join('imgs4', dir_target, 'connectivity_matrix_dendrogram')
	os.makedirs(dir_imgs, exist_ok=True)

	# Load graph.
	d = utils.load(dir_edited_data, prefix, 'connectivity_graph')

	# Make new graphs of CaMKII
	multi_graph_CaMKII, simple_graph_CaMKII, locs_hub, CaMKII_binding_site = \
		utils_graph.make_new_graphs_CaMKII_connectivity(d, nth_largest =  nth_largest)
	G_ = multi_graph_CaMKII
	
	# Calc dendrogram
	Z, labels = calc_dendrogram(G_)
	
	# Set color of dendrogram
	from cycler import cycler
	plt.rcParams["axes.prop_cycle"] = cycler( color=c.cols_list_ratio )
	
	# Plot and save figure
	fig = plt.figure(figsize=(4,4))
	ax1 = fig.add_axes([0.3,0.71,0.6,0.2])
	R=dendrogram(Z, \
		labels=labels, \
		color_threshold = 40, \
		above_threshold_color='k' \
		) # link_color_func=lambda x: 'black' link_color_func=lambda k: c.cols[k]
	plt.setp(ax1.collections,linewidth=1.0)
	ax1.set_title(fig_title)
	
	x_labels = ax1.get_xmajorticklabels()
	
	node_order, blocks, partition, leaves_color_list = get_clusters_from_dendrogram(R)
	ax2 = fig.add_axes([0.3,0.1,0.6,0.6]) # fig.add_axes([0.3,0.1,0.6,0.6])	
	utils_graph.draw_adjacency_matrix(ax2, G_,  node_order = node_order, partitions= blocks, color= 'r')
	#'''
	fig.savefig( os.path.join(dir_imgs, '{}_connect_matrix_dendrogram.svg'.format( fig_title ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}_connect_matrix_dendrogram.png'.format( fig_title ) ) , dpi=150)
	plt.show()
	#'''
	plt.close(fig=fig)
	
	data = {}
	data['multi_graph_CaMKII'] = multi_graph_CaMKII
	data['Z']          = Z
	data['labels']     = labels
	data['node_order'] = node_order
	data['blocks']     = blocks
	data['partition']  = partition
	data['leaves_color_list'] = leaves_color_list
	
	
	utils.save(dir_edited_data, prefix, 'cluster_dendrogram', data)
	
	### graph visulization
	
	'''
	#pos = community_layout(G_, partition)
	cols = [leaves_color_list[i] for i in np.argsort(node_order)]
	
	num_rows    = 1
	num_columns = 3
	fig = plt.figure(figsize=(12,4))
	
	# Rotation
	from scipy.spatial.transform import Rotation as Rot
	r1 = Rot.from_euler('x', 70, degrees=True)
	r2 = Rot.from_euler('z', 20, degrees=True)
	rot = r2 * r1
	pos  = {i: rot.apply(v['loc_hub'])  for i,v in G_.nodes.items()}
	# print('pos ', pos)
	for i, d in enumerate( [[0,1],[1,2],[2,0]] ):
		ax = fig.add_subplot( num_rows, num_columns, i+1 )
		pos_  = {i: np.array( [ v[d[0]], v[d[1]] ] ) for i,v in pos.items() }
		draw_network_simple(ax, G_, pos = pos_, cols = cols)
    
	
	plt.suptitle(fig_title)
	plt.savefig( os.path.join(dir_imgs, '{}_connectivity.svg'.format( fig_title ) ) )
	plt.savefig( os.path.join(dir_imgs, '{}_connectivity.png'.format( fig_title ) ) , dpi=150)
	plt.show()
	'''
	