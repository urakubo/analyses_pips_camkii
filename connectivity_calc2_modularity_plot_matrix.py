
# https://stackoverflow.com/questions/59821151/plot-the-dendrogram-of-communities-found-by-networkx-girvan-newman-algorithm
# https://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-on-top-of-a-matrix-of-data
# https://qiita.com/s-wakaba/items/a93f03f27137cff4a26c

import os, sys, glob, pprint,itertools
import numpy as np


import networkx as nx
from networkx.algorithms.community import modularity


import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph
from specification_datasets import SpecDatasets


import matplotlib.pyplot as plt
plt.rcParams.update(p.rc_param)


def girvan_newman_auto(G, num_div = 10):
	comp = nx.community.girvan_newman(G)
	step = 0
	rcms = []
	for communities in itertools.islice(comp, num_div):
		rcm = tuple(sorted(c) for c in communities)
		rcms.append(rcm)
		print('div step ', step)
		step += 1
	print(rcm)
	
	return rcms, rcm


def get_flat_partitions(rcm):
	a_rcm = list(itertools.accumulate([len(r) for r in rcm], initial=0))
	partitions = [[a_rcm[i], a_rcm[i+1]] for i in range(len(a_rcm)-1) ]
	rcm_flat = utils.flatten(rcm)
	return rcm_flat, partitions
	
	
def process_and_save( d, dir_edited_data, prefix_save ):
	
	
	# Make new graphs of CaMKII
	multi_graph_CaMKII, simple_graph_CaMKII, locs_hub, CaMKII_binding_site = \
		utils_graph.make_new_graphs_CaMKII_connectivity(d, nth_largest = 0)
	G = multi_graph_CaMKII
	
	
	'''
	from networkx.algorithms.community import greedy_modularity_communities
	rcm  = list(greedy_modularity_communities(G))
	rcms = None
	suffix_save='greedy_modularity'
	'''

	#'''
	from networkx.algorithms.community import louvain_communities
	rcm  = list(louvain_communities(G))
	rcms = None
	suffix_save='louvain'
	#'''
	
	print('suffix_save ', suffix_save)
	print('lengths of cluster: ', [len(c) for c in rcm])
	
	rcm_flat, partitions = get_flat_partitions(rcm)
	#plot_matrix_and_save_it(dir_imgs, prefix_save, G, rcm, rcm_flat, partitions, suffix_save)
	
	data = {}
	data['mc_step'] = d['mc_step']
	data['sampling_frame'] = d['sampling_frame']
	data['multi_graph_CaMKII'] = multi_graph_CaMKII
	data['rcms']       = rcms
	data['rcm']        = rcm
	data['rcm_flat']   = rcm_flat
	data['partitions']  = partitions	
	utils.save(dir_edited_data, prefix_save, suffix_save, data)
	
	#sys.exit(0)
	
	
def plot_matrix_and_save_it(dir_imgs, fig_title, G, rcm, rcm_flat, partitions, color= 'r', suffix='connect_matrix'):
	fig = plt.figure(figsize=(4,4))
	ax2 = fig.add_axes([0.3,0.1,0.6,0.6])
	
	modl = modularity(G, rcm)
	
	ax2.set_title(fig_title+', modularity: {:.3f}'.format(modl))
	utils_graph.draw_adjacency_matrix(ax2, G, node_order = rcm_flat, partitions= partitions, color= color)
	fig.savefig( os.path.join(dir_imgs, '{}_{}.svg'.format( fig_title, suffix ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}_{}.png'.format( fig_title, suffix ) ) , dpi=150)
	plt.show()
	#plt.close(fig=fig)
	
	
def plot_matrix_histogram_and_save_it(dir_imgs, fig_title, G, rcm, rcm_flat, partitions, color= 'r', suffix='connect_matrix'):
	fig = plt.figure(figsize=(4,4))
	ax1 = fig.add_axes([0.2,0.1,0.6,0.6])
	
	modl = modularity(G, rcm)
	ax1.set_title(fig_title+', modularity: {:.3f}'.format(modl))
	utils_graph.draw_adjacency_matrix(ax1, G, node_order = rcm_flat, partitions= partitions, color= color)
	
	nums_edges = [G.degree[i] for i in rcm_flat]
	indices = list(range(len(nums_edges)))
	#print('nums_edges ', nums_edges )
	ax2 = fig.add_axes([0.84, 0.1, 0.15, 0.6])
	ax2.barh(indices, nums_edges, nums_edges, color=(0.3,0.3,0.3) )
	ax2.set_ylim([min(indices), max(indices)])
	ax2.set_xlim([0, 20])
	ax2.invert_yaxis()
	ax2.set_xticks([0, 12])
	ax2.set_yticks([])
	fig.savefig( os.path.join(dir_imgs, '{}_{}.svg'.format( fig_title, suffix ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}_{}.png'.format( fig_title, suffix ) ) , dpi=150)
	plt.show()
	
	
	#plt.close(fig=fig)
	
	
def load_and_plot_a_matrix( dir_imgs, dir_edited_data, prefix_load_save, suffix_load_save ):
	
	
	data = utils.load(dir_edited_data, prefix_load_save, suffix_load_save)
	G          = data['multi_graph_CaMKII']
	rcm        = data['rcm']
	rcm_flat   = data['rcm_flat']
	partitions = data['partitions']
	
	from cycler import cycler
	cols   = cycler( color=c.cols_list_ratio )
	colors = [c['color'] for p, c in zip(partitions, cols())]
	#print('colors ', colors)
	#print('partitions ', partitions)
	'''
	plot_matrix_and_save_it(dir_imgs, prefix_load_save, G, rcm, rcm_flat, partitions, \
		color = colors, \
		suffix=suffix_load_save)
	'''
	plot_matrix_histogram_and_save_it(dir_imgs, prefix_load_save, G, rcm, rcm_flat, partitions, \
		color = colors, \
		suffix=suffix_load_save)
	
	
def repeat_for_time_development(dir_edited_data, prefix_load):

	time_frames  = list(range(80)) # 0-5
	
	for time_frame in time_frames:
		
		suffix_load = 'connectivity_graph_{}'.format(time_frame)
		prefix_save = '{}_{}'.format(time_frame, prefix_load)
		print('prefix_save ', prefix_save)
		# Load graph.
		d = utils.load(dir_edited_data, prefix_load, suffix_load)
		process_and_save( d, dir_edited_data, prefix_save )
	
def repeat_for_valency_length(dir_edited_data, prefix_load):
	for prefix_load in prefixes_load:
		suffix_load = 'connectivity_graph'
		prefix_save = prefix_load
		print('prefix_save ', prefix_save)
		# Load graph.
		d = utils.load(dir_edited_data, prefix_load, suffix_load)
		process_and_save( d, dir_edited_data, prefix_save )
	
	
def load_calc_modularity_density(dir_edited_data, prefix_load):
	suffix_load ='louvain'
	data = utils.load(dir_edited_data, prefix_load, suffix_load)
	G          = data['multi_graph_CaMKII']
	rcm        = data['rcm']
	rcm_flat   = data['rcm_flat']
	partitions = data['partitions']
	clust_coeff = nx.average_clustering(nx.Graph(G)) # np.mean(list(nx.clustering(G).values()))
	modul = modularity(G, rcm)
	adj_matrices = [nx.adjacency_matrix(G, nodelist=r).todense() for r in rcm]
	densities    = [np.sum(m) for m in adj_matrices]
	#areas        = [m.shape[0]*m.shape[1] for m in adj_matrices]
	density      = np.sum(densities) / G.number_of_nodes()
	tot_density  = np.sum(nx.adjacency_matrix(G).todense()) / G.number_of_nodes()
	
	return modul, tot_density, density, clust_coeff
	
	
def load_calc_modularities_densities_save_them(dir_edited_data, prefixes_load):
	
	modulatiries = {}
	densities    = {}
	clust_coeffs = {}
	
	for prefix_load in prefixes_load:
		modul, _, density, clust_coeff = \
			load_calc_modularity_density(dir_edited_data, prefix_load)
			
		modulatiries[prefix_load] = modul
		densities[prefix_load]    = density
		clust_coeffs[prefix_load] = clust_coeff
		print('prefix_load ', prefix_load)
		# Load graph.
	
	prefix = 'modularities'
	suffix = 'matrix'
	utils.save(dir_edited_data, prefix, suffix, modulatiries)
	
	prefix = 'densities'
	suffix = 'matrix'
	utils.save(dir_edited_data, prefix, suffix, densities)
	
	prefix = 'average_clustering_coefficient'
	suffix = 'matrix'
	utils.save(dir_edited_data, prefix, suffix, clust_coeffs)


if __name__ == '__main__':
	
	# Input file
	t = SpecDatasets()
	t.valency_length_small_colony2()
	
	dir_imgs = os.path.join(t.dir_imgs_root, 'connectivity_matrix_dendrogram')
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	#load_calc_modularities_densities_save_them(dir_edited_data, prefixes_load)
	#repeat_for_valency_length(dir_edited_data, prefix_load)
	
	
	id_length_valency = 7*4+6 # val_12\R2_006
	id_length_valency = 7*4+2 # val_12\R2_002
	id_length_valency = 7*5+5 # val_12\R2_002
	
	
	prefix_load = t.filenames_edited[id_length_valency]
	suffix_load ='louvain'
	#suffix_load ='greedy_modularity'
	
	
	load_and_plot_a_matrix( dir_imgs, t.dir_edited_data, prefix_load, suffix_load )
	
	
	#plot_matrix_histogram_and_save_it( dir_imgs, dir_edited_data, prefix_load, suffix_load )
	
	
	'''
	modul, tot_density, density = load_calc_modularity_density(t.dir_edited_data, prefix_load)
	print('prefix_load ', prefix_load)
	print('tot_density ', tot_density)
	print('density     ', density)
	
	
	#repeat_for_time_development(dir_edited_data, prefix_load)
	
	
	for prefix_load in prefixes_load:
		print('prefix_load ', prefix_load)
		repeat_for_time_development(t.dir_edited_data, prefix_load)
	'''
	
	
