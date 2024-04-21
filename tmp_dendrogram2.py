
# https://stackoverflow.com/questions/59821151/plot-the-dendrogram-of-communities-found-by-networkx-girvan-newman-algorithm
# https://qiita.com/s-wakaba/items/a93f03f27137cff4a26c

import networkx as nx
from itertools import chain, combinations
import matplotlib.pyplot as plt
from matplotlib import patches

from scipy.cluster.hierarchy import dendrogram
import os
import utils
import colormap as c
import main_all_calc_centrality_working as calc_central


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


def draw_network_simple(multi_graph, prefix):
	# https://qiita.com/yuru/items/4208e13a399773c2b9f8
	plt.figure(figsize=(5,5))
	nx.draw_networkx(	multi_graph, \
						with_labels=False, \
						node_size=25, \
						node_color = 'k', \
						edge_color =c.cmap_universal_ratio['CaMKII'], pos = nx.kamada_kawai_layout(multi_graph)) # pos = nx.kamada_kawai_layout(multi_graph)
	plt.suptitle(prefix)
	plt.show()

if __name__ == '__main__':


	dir_edited_data  =  'small_colony'
	prefix = '01_008'
	prefix = '00_004'
	#prefix = '01_004'
	nth_largest = 0

#	dir_edited_data  = 'valency_length'
#	prefix = '04_002'
	fig_title = '{} th in {}'.format(nth_largest, prefix)
	print(fig_title)
	dir_edited_data  = os.path.join('data3', dir_edited_data)
	d = utils.load(dir_edited_data, prefix, 'connectivity_graph')
	#plot_3D_pvista_CaMKII_interface_region(d)

	# Make new graphs of CaMKII
	multi_graph_CaMKII, simple_graph_CaMKII, locs_hub, CaMKII_binding_site = \
		calc_central.make_new_graphs_CaMKII_connectivity(d, nth_largest =  nth_largest)

	G_ = multi_graph_CaMKII
	#org_nodes = list(G_.nodes)
	#mapping = {o: str(o)+'_' for o in org_nodes}
	#G__ = nx.relabel_nodes(G_, mapping)


	draw_network_simple(G_, fig_title)


	Z, labels = calc_dendrogram(G_)


	#  [list(v)[0] for k, v, in init_node2community_dict.items() if len(v) == 1] # 元のnodesの並び順？
	#  [k for k, v, in init_node2community_dict.items() if len(v) == 1] 

	# pprint.pprint(node_labels) ここにシッカリ対応表がある。
	#  [node_labels[node_id] for node_id in leaves] <= この順に並べなおす。
	# leavesが末端のノードの様だ。
	# 

	# Save するべきは、node_id_to_children, subtree, node_labels, leaves, init_node2community_dict, G など？
	# もっと small scale で試してからの方が良い？

	# dendrogram
	#plt.figure()
	fig = plt.figure(figsize=(6,6))
	
	ax1 = fig.add_axes([0.3,0.71,0.6,0.2])
	#ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
	
	# link_color_func=lambda x: 'black'
	R=dendrogram(Z, labels=labels, color_threshold = 40, above_threshold_color='k')
	plt.title(fig_title)
	
	ax2 = fig.add_axes([0.3,0.1,0.6,0.6]) # fig.add_axes([0.3,0.1,0.6,0.6])
	
	#node_order=[node_labels[node_id] for node_id in leaves]
	#print('node oredr: ',  R['ivl'])
	
	leaves_color_list = R['leaves_color_list']
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
	
	
	# blocks = [l.index(i) for i in sorted(list(set(l)))]+[len(l)]
	
	print('blocks ', blocks)
	
	draw_adjacency_matrix(ax2, G_,  node_order =  R['ivl'], partitions= blocks, color= 'r')
	
	plt.savefig('dendrogram.png')
	plt.show()

