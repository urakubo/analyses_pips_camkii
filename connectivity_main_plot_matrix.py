
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


if __name__ == '__main__':

	# Input file
	i = 2
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
	
	

	print('girvan_newman')
	# rcm, rcm_flat = utils_graph.girvan_newman_by_hand(G, num_div = 4)
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
	
	
	# Plot and save figure
	fig = plt.figure(figsize=(4,4))
	
	ax2 = fig.add_axes([0.3,0.1,0.6,0.6])
	ax2.set_title(fig_title)
	utils_graph.draw_adjacency_matrix(ax2, G_, node_order = rcm_flat, partitions= partitions, color= 'r')
	fig.savefig( os.path.join(dir_imgs, '{}_connect_matrix.svg'.format( fig_title ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}_connect_matrix.png'.format( fig_title ) ) , dpi=150)
	plt.show()
	plt.close(fig=fig)
	
	
	data = {}
	data['multi_graph_CaMKII'] = multi_graph_CaMKII
	data['rcm']        = rcm
	data['rcm_flat']   = rcm_flat
	data['partitions']  = partitions
	utils.save(dir_edited_data, prefix, 'cluster', data)
	
