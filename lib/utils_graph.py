import os, sys, glob, pprint,itertools
import numpy as np


import networkx as nx
from networkx.algorithms.community.centrality import girvan_newman

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c


import matplotlib.pyplot as plt
from matplotlib import patches
plt.rcParams.update(p.rc_param)


def get_angle(x, y):

	dot_xy = np.dot(x, y)
	norm_x = np.linalg.norm(x)
	norm_y = np.linalg.norm(y)
	cos = dot_xy / (norm_x*norm_y)
	rad = np.arccos(cos)
	theta = rad * 180 / np.pi

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


def _position_communities(g, partition, **kwargs):

    # create a weighted graph, in which each node corresponds to a community,
    # and each edge weight to the number of edges between communities
    between_community_edges = _find_between_community_edges(g, partition)

    communities = set(partition.values())
    hypergraph = nx.DiGraph()
    hypergraph.add_nodes_from(communities)
    for (ci, cj), edges in between_community_edges.items():
        hypergraph.add_edge(ci, cj, weight=len(edges))

    # find layout for communities
    pos_communities = nx.spring_layout(hypergraph, **kwargs)

    # set node positions to position of community
    pos = dict()
    for node, community in partition.items():
        pos[node] = pos_communities[community]

    return pos

def _find_between_community_edges(g, partition):

    edges = dict()

    for (ni, nj) in g.edges():
        ci = partition[ni]
        cj = partition[nj]

        if ci != cj:
            try:
                edges[(ci, cj)] += [(ni, nj)]
            except KeyError:
                edges[(ci, cj)] = [(ni, nj)]

    return edges

def _position_nodes(g, partition, **kwargs):
    """
    Positions nodes within communities.
    """

    communities = dict()
    for node, community in partition.items():
        try:
            communities[community] += [node]
        except KeyError:
            communities[community] = [node]

    pos = dict()
    for ci, nodes in communities.items():
        subgraph = g.subgraph(nodes)
        pos_subgraph = nx.spring_layout(subgraph, **kwargs)
        pos.update(pos_subgraph)

    return pos



# https://stackoverflow.com/questions/43541376/how-to-draw-communities-with-networkx
def community_layout(g, partition):
	"""
	Compute the layout for a modular graph.


	Arguments:
	----------
	g -- networkx.Graph or networkx.DiGraph instance
	    graph to plot

	partition -- dict mapping int node -> int community
	    graph partitions


	Returns:
	--------
	pos -- dict mapping int node -> (float x, float y)
	    node positions

	"""
	pos_communities = _position_communities(g, partition, scale=3.)
	pos_nodes = _position_nodes(g, partition, scale=1.)

	# combine positions
	pos = dict()
	for node in g.nodes():
	    pos[node] = pos_communities[node] + pos_nodes[node]
	
	return pos
	

def make_new_graphs_CaMKII_connectivity(d, nth_largest = 0):

	# Use the multi graph.
	#g = d['simple_graph_CaMKII_GluN2B']
	g = d['multi_graph']
	
	# Pickup the largest cluster
	clusters = sorted(nx.connected_components(g), key=len, reverse=True)
	lengths_clusters = [len(c) for c in clusters]
	#'''
	print('lengths_clusters ', lengths_clusters[0:10])
	print('Picked-up cluster ', lengths_clusters[nth_largest] )
	#'''
	
	g_largest_cluster = nx.MultiGraph( g.subgraph(clusters[nth_largest]) )
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
			multi_graph_CaMKII.add_node(id, id_bead = id_bead, id_bead_all=v['ids_bead'][0], loc_hub = position_CaMKII_hub)
			simple_graph_CaMKII.add_node(id, id_bead = id_bead, id_bead_all=v['ids_bead'][0], loc_hub = position_CaMKII_hub)
			locs_hub.append( position_CaMKII_hub )
			
			# Make CaMKII interaction beads
			id_CaMKII_binding = np.nonzero(v['types_bead'] == p.subunits['CaMKII binding site']['id'])
			position_CaMKII_binding = v['positions_grid_coord'][id_CaMKII_binding, :][0]
			ids_bead       = v['ids_bead'][0][id_CaMKII_binding]
			
			#print('position_CaMKII_binding: ', position_CaMKII_binding, ', ids_bead: ',ids_bead)
			for id_bead_, loc in zip(ids_bead, position_CaMKII_binding):
				vector_to_hub = position_CaMKII_hub - loc
				
				CaMKII_binding_site[id_bead_] = {}
				CaMKII_binding_site[id_bead_]['loc'] = loc
				CaMKII_binding_site[id_bead_]['connect'] = 0
				CaMKII_binding_site[id_bead_]['id_molecule'] = id
				CaMKII_binding_site[id_bead_]['vector_to_hub'] = vector_to_hub
				CaMKII_binding_site[id_bead_]['distance_to_hub'] = np.linalg.norm( vector_to_hub )
				
				CaMKII_binding_site[id_bead_]['angle'] = get_angle(vector_to_hub, -loc)
				
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
	
	#'''
	print('Num of GluN2B (edges): ', len(ids_GluN2B))
	print('Num of CaMKII (nodes): ', len(ids_CaMKII))
	#'''
	return multi_graph_CaMKII, simple_graph_CaMKII, locs_hub, CaMKII_binding_site

	# girvan_newman1
def girvan_newman_by_hand(G, num_div = 4):
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
	

