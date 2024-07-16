import os, sys, glob, pprint,itertools
import numpy as np


import networkx as nx
from networkx.algorithms.community.centrality import girvan_newman

import skimage
from scipy import ndimage
from scipy.optimize import curve_fit

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c


import matplotlib.pyplot as plt
from matplotlib import patches
from cycler import cycler
plt.rcParams.update(p.rc_param)



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


def _draw_a_patch(ax, p, c):
	p_start = p[0]
	p_end   = p[1]
	p_size  = p_end - p_start
	ax.add_patch(patches.Rectangle((p_start-0.5, p_start-0.5),
	                              p_size, # Width
	                              p_size, # Height
	                              facecolor=None,
	                              fill=False,
	                              edgecolor=c,
	                              linewidth=1))
	
	
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
	print('Patches drawing')
	if partitions is not None and not isinstance(color, (list, tuple)):
		for p in partitions:
			_draw_a_patch(ax, p, color)
	elif partitions is not None and isinstance(color, (list, tuple)):
		if (partitions[0][0] == 0) and (partitions[0][1] == 0):
			partitions = partitions[1:]
		for p, c in zip(partitions, color) :
			## print('patch {}, {}, color: '.format(p[0], p[1], c) )
			_draw_a_patch(ax, p, c)

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
	


def get_graphs(ids_molecule, types, bp, positions_grid_coord):
	
	unique_ids_molecule = np.unique( ids_molecule )
	multi_graph = nx.MultiGraph()
	simple_graph_CaMKII_GluN2B = nx.Graph()
	# CaMKII_or_GluN2B = np.zeros( len(ids_molecule), dtype = 'int')
	
	# Create nodes
	for id in unique_ids_molecule:
		flags = (id == ids_molecule)
		ids_bead         = np.nonzero(flags)
		types_bead       = types[flags]
		positions_bead   = positions_grid_coord[flags,:]
		
		ids_partner_bead = bp[flags]
		species = [k for k, i in p.molecules_without_all.items() if types_bead[0] in i['id']][0]
		
		
		multi_graph.add_node(id,\
			species    = species, \
			ids_bead   = ids_bead, \
			types_bead = types_bead, \
			num_beads  = np.sum(flags), \
			positions_grid_coord = positions_bead, \
			ids_partner_bead =ids_partner_bead)
		
		if species in ['CaMKII','GluN2B']:
			simple_graph_CaMKII_GluN2B.add_node(id, \
				species    = species, \
				ids_bead   = ids_bead, \
				types_bead = types_bead, \
				num_beads  = np.sum(flags), \
				positions_grid_coord = positions_bead, \
				ids_partner_bead =ids_partner_bead)
		
	# Make connection
	list_ids = list(multi_graph.nodes.keys())
	ids_bead_already_connected = np.zeros_like( bp, dtype='int' )
	for i in range(ids_molecule.shape[0]):
		if 	(bp[i] >= 0) and \
			(ids_bead_already_connected[i] == 0) and \
			(ids_bead_already_connected[bp[i]] == 0):
			
			ids_bead_already_connected[i]     = 1
			ids_bead_already_connected[bp[i]] = 1
			
			id_molecule = ids_molecule[i]
			id_molecule_partner = ids_molecule[bp[i]]
			
			species = multi_graph.nodes[id_molecule]['species']
			species_partner = multi_graph.nodes[id_molecule_partner]['species']
			connecting_species = [species, species_partner]
			if ('GluN2B' in connecting_species) and ('CaMKII' in connecting_species):
				type_connection = 'GluN2B_CaMKII'
			elif ('GluN2B' in connecting_species) and ('PSD95' in connecting_species):
				type_connection = 'GluN2B_PSD95'
			elif ('STG' in connecting_species) and ('PSD95' in connecting_species):
				type_connection = 'STG_PSD95'
			else:
				raise ValueError("Erronous connection: {}", connecting_species)
			multi_graph.add_edge(id_molecule, id_molecule_partner, type_connection = type_connection, id_bead1 = i ,id_bead2 = bp[i])
			
			if type_connection == 'GluN2B_CaMKII':
				simple_graph_CaMKII_GluN2B.add_edge(id_molecule, id_molecule_partner, type_connection = type_connection, id_bead1 = i ,id_bead2 = bp[i])
	
	return multi_graph, simple_graph_CaMKII_GluN2B



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
	

def get_a_profile( data, i_split = 0 ):
	
	org_num_clusters  = len(data[0]['rcm'])
	# print('org_num_clusters ', org_num_clusters )
	
	rcm_t0_0     = set(data[0]['rcm'][i_split])
	num_rcm_t0_0 = len(rcm_t0_0)
	
	tot_num_CaMKII =  sum( [len(r) for r in data[0]['rcm']] )
	ave_ratio = num_rcm_t0_0 / tot_num_CaMKII
	
	prof = []
	# print('len(rcm_t0_0) ', num_rcm_t0_0)
	for d in data.values():
		num_clusters_t = len(d['rcm'])
		# print('num_clusters_t ', num_clusters_t)
		rcm_0 = [set(d['rcm'][c]) for c in range(num_clusters_t)]
		rcm_0 = [len(rcm_t0_0 & r)/len(r) for r in rcm_0]
		#print('rcm_0 ', rcm_0)
		prof.append(rcm_0)
	
	return prof, ave_ratio


def get_time( data ):
	mc_step_0   = data[0]['mc_step']
	time_points = [d['mc_step'] - mc_step_0 for d in data.values()]
	time_points = np.array(time_points)
	return time_points


def prep_fig(ax, time_points):
	xmin = 0
	xmax = np.max(time_points)
	ax.set_xlabel('Time (/10^9 MC steps)')
	ax.set_ylabel('Standard deviation of cluster-wise mixture ratios')
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_xlim([xmin, xmax])
	ax.set_ylim([-0.1, 0.5])
	
	return xmin, xmax
	
	
def func_exponential(x, tau, a, b):
    return a*(np.exp(-x/tau) + b)

def get_profiles(data):
	# Num of cluster
	partitions  = data[0]['partitions']
	num_cluster = len(partitions)
	
	# Decode time
	time_points = get_time( data ) / 1e9
	
	# Set colors
	cols   = cycler( color=c.cols_list_ratio )
	colors = [col['color'] for p, col in zip(partitions, cols())]
	
	profs  = []
	taus   = []
	params = []
	for k, col in enumerate( colors ):
		prof, ave_ratio = get_a_profile( data, k )
		prof = np.array([np.std(p) for p in prof])
		# exponential fit: func_exponential(x, tau, a, b):
		min_tau, max_tau = 0, 20
		min_a  , max_a   = 0, 0.3
		min_b  , max_b   = 0, 1.0
		param, cov = curve_fit(func_exponential, time_points[1:], prof[1:],\
			p0=[1, 0.2, 0.0 ],\
			bounds = ((min_tau, min_a, min_b), (max_tau, max_a, max_b)),\
			maxfev=5000)
		profs.append(prof)
		params.append(param)
		taus.append(param[0])
		
	return time_points, profs, taus, params, colors
	
	
	
def get_multi_graph(ids_molecule, types, bp, positions_grid_coord):
	
	unique_ids_molecule = np.unique( ids_molecule )
	multi_graph = nx.MultiGraph()
	# CaMKII_or_GluN2B = np.zeros( len(ids_molecule), dtype = 'int')
	
	# Create nodes
	for id in unique_ids_molecule:
		flags = (id == ids_molecule)
		ids_bead         = np.nonzero(flags)
		types_bead       = types[flags]
		positions_bead   = positions_grid_coord[flags,:]
		
		ids_partner_bead = bp[flags]
		species = [k for k, i in p.molecules_without_all.items() if types_bead[0] in i['id']][0]
		multi_graph.add_node(id,\
			species    = species, \
			ids_bead   = ids_bead, \
			types_bead = types_bead, \
			num_beads  = np.sum(flags), \
			positions_grid_coord = positions_bead, \
			ids_partner_bead =ids_partner_bead)
	
		
	# Make connection
	list_ids = list(multi_graph.nodes.keys())
	ids_bead_already_connected = np.zeros_like( bp, dtype='int' )
	for i in range(ids_molecule.shape[0]):
		if 	(bp[i] >= 0) and \
			(ids_bead_already_connected[i] == 0) and \
			(ids_bead_already_connected[bp[i]] == 0):
			
			ids_bead_already_connected[i]     = 1
			ids_bead_already_connected[bp[i]] = 1
			
			id_molecule = ids_molecule[i]
			id_molecule_partner = ids_molecule[bp[i]]
			
			species = multi_graph.nodes[id_molecule]['species']
			species_partner = multi_graph.nodes[id_molecule_partner]['species']
			connecting_species = [species, species_partner]
			if ('GluN2B' in connecting_species) and ('CaMKII' in connecting_species):
				type_connection = 'GluN2B_CaMKII'
			elif ('GluN2B' in connecting_species) and ('PSD95' in connecting_species):
				type_connection = 'GluN2B_PSD95'
			elif ('STG' in connecting_species) and ('PSD95' in connecting_species):
				type_connection = 'STG_PSD95'
			else:
				raise ValueError("Erronous connection: {}", connecting_species)
			multi_graph.add_edge(id_molecule, id_molecule_partner, type_connection = type_connection, id_bead1 = i ,id_bead2 = bp[i])
			
	
	return multi_graph
	
	
def get_connection_statistics_condensate(multi_graph, species, type_analysis):
	
	clusters = sorted(nx.connected_components(multi_graph), key=len, reverse=True)
	g_largest_cluster = multi_graph.subgraph(clusters[0])
	
	return get_connection_statistics(g_largest_cluster, species, type_analysis)
	
	
	
def get_connection_statistics(multi_graph, species, type_analysis):
	
	ids = [ i for i, attr in multi_graph.nodes('species') if attr == species ]
	numbers_of_connections = [ multi_graph.degree[id] for id in ids ]
	
	edges_from_molecules = [ multi_graph.edges(id, keys=True, data=True) for id in ids ]
	connections_from_one_molecules = [ [e[3]['type_connection'] for e in es] for es in edges_from_molecules ]
	
	
	if species == 'GluN2B' and type_analysis in ['CaMKII','PSD95']:
		types_connection = [[len(v) for k, v in multi_graph[id].items() if multi_graph.nodes[k]['species'] == type_analysis] for id in ids]
		reference_types_connection = {\
			'No {}'.format(type_analysis) :[],\
			'One {} (single)'.format(type_analysis) :[1],\
			'One {} (double)'.format(type_analysis) :[2],\
			'Two {}s'.format(type_analysis) :[1, 1]\
			}
		dist = {k: types_connection.count(v) for k, v in reference_types_connection.items()}
		
	elif species in ['CaMKII','STG'] and type_analysis == 'distribution':
		partner      = {'CaMKII':'GluN2B', 'STG': 'PSD95'}
		binding_num  = {'CaMKII':13, 'STG':5}
		dist = {'{:d} {}'.format(i, partner[species]): numbers_of_connections.count(i) for i in range(binding_num[species])}
		
	elif species in ['CaMKII','STG'] and type_analysis == 'average':
		partner = {'CaMKII':'GluN2B', 'STG': 'PSD95'}
		ave_numbers_of_connections = np.average( numbers_of_connections )
		dist = {partner[species]: ave_numbers_of_connections}
		
	elif species == 'GluN2B' and type_analysis == 'average':
		num_GluN2B_CaMKII = np.average( [ c.count('GluN2B_CaMKII') for c in connections_from_one_molecules ] )
		num_GluN2B_PSD95  = np.average( [ c.count('GluN2B_PSD95') for c in connections_from_one_molecules ] )
		dist = {'GluN2B_CaMKII': num_GluN2B_CaMKII,
				'GluN2B_PSD95' : num_GluN2B_PSD95 }
		
	elif species == 'GluN2B' and type_analysis == 'distribution':
		reference_types_connection = [\
			[], \
			['GluN2B_CaMKII'], \
			['GluN2B_CaMKII', 'GluN2B_CaMKII'], \
			['GluN2B_PSD95'], \
			['GluN2B_PSD95', 'GluN2B_PSD95'], \
			['GluN2B_CaMKII', 'GluN2B_PSD95']]
		titles = [\
			'None', \
			'CaMKII, None',\
			'CaMKII, CaMKII',\
			'PSD95, None', \
			'PSD95, PSD95', \
			'CaMKII, PSD95' ]
		dist = {t: 0 for t in titles}
		for connections_from_a_molecule in connections_from_one_molecules:
			for j, ref in enumerate( reference_types_connection ):
				if utils.equal_list(connections_from_a_molecule, ref):
					dist[titles[j]] += 1
		
	elif species == 'PSD95' and type_analysis == 'average':
		num_STG_PSD95    = np.average( [c.count('STG_PSD95')  for c in connections_from_one_molecules ]   )
		num_GluN2B_PSD95 = np.average( [c.count('GluN2B_PSD95') for c in connections_from_one_molecules ] )
		dist = {'STG_PSD95': num_STG_PSD95, \
				'GluN2B_PSD95': num_GluN2B_PSD95 }
		
	elif species == 'PSD95' and type_analysis == 'distribution':
		reference_types_connection = [\
			[], \
			['STG_PSD95'],['STG_PSD95']*2,['STG_PSD95']*3, \
			['GluN2B_PSD95'], ['GluN2B_PSD95']*2, ['GluN2B_PSD95']*3,\
			['STG_PSD95', 'GluN2B_PSD95'],\
			['STG_PSD95']*2 + ['GluN2B_PSD95'],\
			['STG_PSD95'] + ['GluN2B_PSD95']*2]
		titles = [\
			'None', \
			'1 STG', '2 STG', '3 STG', \
			'1 GluN2B', '2 GluN2B', '3 GluN2B', \
			'1 STG, 1 GluN2B', '2 STG, 1 GluN2B', '1 STG, 2 GluN2B'
			 ]
		dist = {t: 0 for t in titles}
		for connections_from_a_molecule in connections_from_one_molecules:
			for j, ref in enumerate( reference_types_connection ):
				if utils.equal_list(connections_from_a_molecule, ref):
					dist[titles[j]] += 1
					
	elif species == 'PSD95' and type_analysis == 'ratio':
		types_connection = []
		#print('connections_from_one_molecules')
		#print(connections_from_one_molecules)
		for c in connections_from_one_molecules:
			if ('STG_PSD95' not in c) and ('GluN2B_PSD95' not in c):
				types_connection.append('None')
			elif ('STG_PSD95' in c) and ('GluN2B_PSD95' not in c):
				types_connection.append('STG only')
			elif ('STG_PSD95' not in c) and ('GluN2B_PSD95' in c):
				types_connection.append('GluN2B only')
			elif ('STG_PSD95' in c) and ('GluN2B_PSD95' in c):
				types_connection.append('Both')
		dist = {t: types_connection.count(t) for t in [ 'None', 'STG only', 'GluN2B only', 'Both' ]}
	
	return dist
	
	
	
def get_surface_condensate_CaMKII(types, positions, ids_molecule, sigma=2):
	
	# Parameters
	
	positions_in_grid_coord = np.floor(positions + p.center_np).astype('int')
	locs_in_grid_mesh = \
		{v['id']:  np.ones([len(p.edge0),len(p.edge1),len(p.edge2)], dtype='int')*(-1) for v in p.subunits.values()}
	
	for i, (t, pos) in enumerate( zip(types, positions_in_grid_coord) ):
		locs_in_grid_mesh[t][pos[0],pos[1],pos[2]] = i
	
	# Define the region periphery.
	region_periphery = utils.get_periphery_in_grid_mesh()
	
	# Get the CaMKII condenate regions.
	conc_CaMKII_in_grid_mesh = [ locs_in_grid_mesh[id] >= 0 for id in p.molecules_without_all['CaMKII']['id'] ]
	conc_CaMKII_in_grid_mesh = sum(conc_CaMKII_in_grid_mesh).astype('float')
	conc_CaMKII_in_grid_mesh = ndimage.gaussian_filter(conc_CaMKII_in_grid_mesh, sigma = sigma)
	condensate_CaMKII_in_grid_mesh = utils.get_high(conc_CaMKII_in_grid_mesh)


	# Get the condensate of all.
	conc_all_in_grid_mesh = [ v >= 0 for v in locs_in_grid_mesh.values()]
	conc_all_in_grid_mesh = sum( conc_all_in_grid_mesh ).astype('float')
	
	# print('conc_all_in_grid_mesh.shape : ', conc_all_in_grid_mesh.shape ) 
	conc_all_in_grid_mesh = ndimage.gaussian_filter(conc_all_in_grid_mesh, sigma = sigma)
	conc_all_periphery    = np.sum(conc_all_in_grid_mesh * region_periphery ) / np.sum( region_periphery )
	condensate_all_in_grid_mesh = utils.get_high(conc_all_in_grid_mesh-conc_all_periphery)


	# Get the interface region
	footprint_inside  = skimage.morphology.ball(1) #### "2": r = 3, "1": r = 2
	footprint_outside = skimage.morphology.ball(15) #### "2": r = 3, "1": r = 2
	dil1 = skimage.morphology.binary_dilation( condensate_CaMKII_in_grid_mesh, footprint=footprint_outside )
	ero1 = skimage.morphology.binary_erosion( condensate_CaMKII_in_grid_mesh, footprint=footprint_inside )
	region_interface_CaMKII = np.logical_xor(dil1, ero1)
	
	dil2 = skimage.morphology.binary_dilation( condensate_all_in_grid_mesh, footprint=footprint_outside  )
	ero2 = skimage.morphology.binary_erosion( condensate_all_in_grid_mesh , footprint=footprint_inside )
	region_interface_all = np.logical_xor(dil2, ero2)
	
	beads_interface_CaMKII = {v['id']: locs_in_grid_mesh[v['id']][region_interface_CaMKII] for v in p.subunits.values() } 
	beads_interface_CaMKII = {k: v[v != -1] for k, v in beads_interface_CaMKII.items() } 
	beads_interface_all    = {v['id']: locs_in_grid_mesh[v['id']][region_interface_all] for v in p.subunits.values() } 
	beads_interface_all    = {k: v[v != -1] for k, v in beads_interface_all.items() } 
	
	# Summary
	d = {
			'locs_in_grid_mesh'				: locs_in_grid_mesh,
			'conc_CaMKII_in_grid_mesh'		: conc_CaMKII_in_grid_mesh,
			'condensate_CaMKII_in_grid_mesh': condensate_CaMKII_in_grid_mesh,
			'conc_all_in_grid_mesh'			: conc_all_in_grid_mesh,
			'condensate_all_in_grid_mesh'	: condensate_all_in_grid_mesh,
			'region_interface_CaMKII'		: region_interface_CaMKII,
			'region_interface_all'			: region_interface_all,
			'beads_interface_CaMKII'		: beads_interface_CaMKII,
			'beads_interface_all'			: beads_interface_all,
		}
	return d
	
	
	
	
def get_graphs(ids_molecule, types, bp, positions_grid_coord):
	
	unique_ids_molecule = np.unique( ids_molecule )
	multi_graph = nx.MultiGraph()
	simple_graph_CaMKII_GluN2B = nx.Graph()
	# CaMKII_or_GluN2B = np.zeros( len(ids_molecule), dtype = 'int')
	
	# Create nodes
	for id in unique_ids_molecule:
		flags = (id == ids_molecule)
		ids_bead         = np.nonzero(flags)
		types_bead       = types[flags]
		positions_bead   = positions_grid_coord[flags,:]
		
		ids_partner_bead = bp[flags]
		species = [k for k, i in p.molecules_without_all.items() if types_bead[0] in i['id']][0]
		
		
		multi_graph.add_node(id,\
			species    = species, \
			ids_bead   = ids_bead, \
			types_bead = types_bead, \
			num_beads  = np.sum(flags), \
			positions_grid_coord = positions_bead, \
			ids_partner_bead =ids_partner_bead)
		
		if species in ['CaMKII','GluN2B']:
			simple_graph_CaMKII_GluN2B.add_node(id, \
				species    = species, \
				ids_bead   = ids_bead, \
				types_bead = types_bead, \
				num_beads  = np.sum(flags), \
				positions_grid_coord = positions_bead, \
				ids_partner_bead =ids_partner_bead)
		
	# Make connection
	list_ids = list(multi_graph.nodes.keys())
	ids_bead_already_connected = np.zeros_like( bp, dtype='int' )
	for i in range(ids_molecule.shape[0]):
		if 	(bp[i] >= 0) and \
			(ids_bead_already_connected[i] == 0) and \
			(ids_bead_already_connected[bp[i]] == 0):
			
			ids_bead_already_connected[i]     = 1
			ids_bead_already_connected[bp[i]] = 1
			
			id_molecule = ids_molecule[i]
			id_molecule_partner = ids_molecule[bp[i]]
			
			species = multi_graph.nodes[id_molecule]['species']
			species_partner = multi_graph.nodes[id_molecule_partner]['species']
			connecting_species = [species, species_partner]
			if ('GluN2B' in connecting_species) and ('CaMKII' in connecting_species):
				type_connection = 'GluN2B_CaMKII'
			elif ('GluN2B' in connecting_species) and ('PSD95' in connecting_species):
				type_connection = 'GluN2B_PSD95'
			elif ('STG' in connecting_species) and ('PSD95' in connecting_species):
				type_connection = 'STG_PSD95'
			else:
				raise ValueError("Erronous connection: {}", connecting_species)
			multi_graph.add_edge(id_molecule, id_molecule_partner, type_connection = type_connection, id_bead1 = i ,id_bead2 = bp[i])
			
			if type_connection == 'GluN2B_CaMKII':
				simple_graph_CaMKII_GluN2B.add_edge(id_molecule, id_molecule_partner, type_connection = type_connection, id_bead1 = i ,id_bead2 = bp[i])
	
	return multi_graph, simple_graph_CaMKII_GluN2B
	
	
	
	