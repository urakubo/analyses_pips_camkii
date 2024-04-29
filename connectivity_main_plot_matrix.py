
# https://stackoverflow.com/questions/59821151/plot-the-dendrogram-of-communities-found-by-networkx-girvan-newman-algorithm
# https://qiita.com/s-wakaba/items/a93f03f27137cff4a26c

import os, sys, glob, pprint,itertools
import numpy as np


import networkx as nx


import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph


import matplotlib.pyplot as plt
plt.rcParams.update(p.rc_param)


def plot_figure_and_save_it(dir_imgs, fig_title, G_, rcm_flat, partitions):
	fig = plt.figure(figsize=(4,4))
	ax2 = fig.add_axes([0.3,0.1,0.6,0.6])
	ax2.set_title(fig_title)
	utils_graph.draw_adjacency_matrix(ax2, G_, node_order = rcm_flat, partitions= partitions, color= 'r')
	fig.savefig( os.path.join(dir_imgs, '{}_connect_matrix.svg'.format( fig_title ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}_connect_matrix.png'.format( fig_title ) ) , dpi=150)
	plt.show()
	#plt.close(fig=fig)

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
	
	a_rcm = list(itertools.accumulate([len(r) for r in rcm], initial=0))
	partitions = [[a_rcm[i], a_rcm[i+1]] for i in range(len(a_rcm)-1) ]
	rcm_flat = utils.flatten(rcm)
	
	return rcms, rcm, rcm_flat, partitions
	
	
def save_data(dir_edited_data, prefix_save, d, multi_graph_CaMKII, rcms, rcm, rcm_flat, partitions):
	data = {}
	data['mc_step'] = d['mc_step']
	data['sampling_frame'] = d['sampling_frame']
	data['multi_graph_CaMKII'] = multi_graph_CaMKII
	data['rcms']       = rcms
	data['rcm']        = rcm
	data['rcm_flat']   = rcm_flat
	data['partitions']  = partitions
	utils.save(dir_edited_data, prefix_save, 'cluster', data)


	
if __name__ == '__main__':

	# Input file
	dir_target  = 'small_colony'
	prefixes    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(3) for id_f in range(7)]
	prefix_load = prefixes[2]
	suffix_load = 'connectivity_graph'

	# Shared init
	dir_edited_data  = os.path.join('data4', dir_target)
	dir_imgs         = os.path.join('imgs4', dir_target, 'connectivity_matrix_dendrogram')
	os.makedirs(dir_imgs, exist_ok=True)
	# 
	time_frames  = [0,1,2,3,4,5] # 0-5
	
	for time_frame in time_frames:
		
		suffix_load = 'connectivity_graph_{}'.format(time_frame)
		prefix_save = '{}_{}'.format(time_frame, prefix_load)
		print('prefix_save ', prefix_save)
		print('suffix_load ', suffix_load)
		
		# Load graph.
		d = utils.load(dir_edited_data, prefix_load, suffix_load)
		
		# Make new graphs of CaMKII
		multi_graph_CaMKII, simple_graph_CaMKII, locs_hub, CaMKII_binding_site = \
			utils_graph.make_new_graphs_CaMKII_connectivity(d, nth_largest = 0)
		G = multi_graph_CaMKII
		
		
		# rcm, rcm_flat = utils_graph.girvan_newman_by_hand(G, num_div = 4)
		
		print('girvan_newman_auto')
		rcms, rcm, rcm_flat, partitions = girvan_newman_auto(G, num_div = 10)
		
		plot_figure_and_save_it(dir_imgs, prefix_save, G, rcm_flat, partitions)
		
		save_data(dir_edited_data, prefix_save, d, multi_graph_CaMKII, rcms, rcm, rcm_flat, partitions)
		
