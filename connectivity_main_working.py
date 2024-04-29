
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



if __name__ == '__main__':


	# Input file
	
	i = 2
	dir_target  = 'small_colony'
	prefixes    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(3) for id_f in range(7)]
	prefix_load = prefixes[i]
	suffix_load = 'connectivity_graph'
	time_frame  = 0
	
	# 
	time_frames = 6
	suffix_load = 'connectivity_graph_{}'.format(time_frame)
	
	
	# Shared init
	dir_edited_data  = os.path.join('data4', dir_target)
	'''
	print(prefix_save)
	dir_imgs         = os.path.join('imgs4', dir_target, 'connectivity_matrix_dendrogram')
	os.makedirs(dir_imgs, exist_ok=True)
	'''
	
	data = {}
	for time_frame in range(time_frames):
		
		# Load graph.
		prefix_load_ = '{}_{}'.format(time_frame, prefix_load)
		data[time_frame] = utils.load(dir_edited_data, prefix_load_, 'cluster')
		
		'''
		multi_graph_CaMKII = data['multi_graph_CaMKII']
		Z          = data['Z']
		labels     = data['labels']
		node_order = data['node_order']
		blocks     = data['blocks']
		partition  = data['partition']
		leaves_color_list = data['leaves_color_list']
		'''
		
		'''
		multi_graph_CaMKII = data['multi_graph_CaMKII']
		rcm        = data['rcm']
		rcm_flat   = data['rcm_flat']
		partitions = data['partitions']
		'''
