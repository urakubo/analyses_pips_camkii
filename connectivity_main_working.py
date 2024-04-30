
# https://stackoverflow.com/questions/59821151/plot-the-dendrogram-of-communities-found-by-networkx-girvan-newman-algorithm
# https://qiita.com/s-wakaba/items/a93f03f27137cff4a26c

import os, sys, glob, pprint,itertools
import numpy as np


# from itertools import chain, combinations
import matplotlib.pyplot as plt

#import networkx as nx
#from scipy.cluster.hierarchy import dendrogram
#from networkx.algorithms.community.centrality import girvan_newman

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
	num_time_frames = 40
	suffix_load = 'connectivity_graph_{}'.format(time_frame)
	
	
	# Shared init
	dir_edited_data  = os.path.join('data4', dir_target)
	'''
	print(prefix_save)
	dir_imgs         = os.path.join('imgs4', dir_target, 'connectivity_matrix_dendrogram')
	os.makedirs(dir_imgs, exist_ok=True)
	'''
	
	data = {}
	for time_frame in range(num_time_frames):
		
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
	
	
	total_num_CaMKII = len(utils.flatten(data[0]['rcms'][0]))
	print('total_num_CaMKII ', total_num_CaMKII )
	num_splits = 4 # 2-10
	
	profs = []
	for i_split in range(num_splits):
		prof = []
		rcm_t0_0 = set(data[0]['rcms'][num_splits-2][i_split])
		num_rcm_t0_0 = len(rcm_t0_0)
		print('len(rcm_t0_0) ', num_rcm_t0_0)
		for time_frame in range(num_time_frames):
			rcm_0 = [set(data[time_frame]['rcms'][num_splits-2][j_split]) for j_split in range(num_splits)]
			rcm_0 = [len(rcm_t0_0 & r)/num_rcm_t0_0 for r in rcm_0]
			#print('rcm_0 ', rcm_0)
			prof.append(max(rcm_0))
		profs.append(prof)
	
	time = list(range(num_time_frames))
	num_rows    = 2
	num_columns = 2

	fig  = plt.figure(figsize=(8, 8))
	
	ax = fig.add_subplot( num_rows, num_columns, 1 )
	# ax.set_title(title)

	ax.set_xlabel('Time (/10 frame)')
	ax.set_ylabel('Num overalpped voxels')
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_xlim([np.min(time), np.max(time)])
	
	for prof in profs:
		ax.plot(time, prof, '-')

	ax.set_ylim(bottom=0)
	#fig.savefig( os.path.join(self.dir_imgs, self.basename + '.svg' ) )
	#fig.savefig( os.path.join(self.dir_imgs, self.basename + '.png' ), dpi=150 )
	plt.show()

