
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

def get_profiles( data ):
	
	org_num_clusters  = len(data[0]['rcm'])
	print('org_num_clusters ', org_num_clusters )
	profs = []
	for i_split in range(org_num_clusters):
		prof = []
		rcm_t0_0 = set(data[0]['rcm'][i_split])
		num_rcm_t0_0 = len(rcm_t0_0)
		# print('len(rcm_t0_0) ', num_rcm_t0_0)
		for time_frame in range(num_time_frames):
			num_clusters_t = len(data[time_frame]['rcm'])
			# print('num_clusters_t ', num_clusters_t)
			rcm_0 = [set(data[time_frame]['rcm'][c]) for c in range(num_clusters_t)]
			rcm_0 = [len(rcm_t0_0 & r) for r in rcm_0]
			#print('rcm_0 ', rcm_0)
			prof.append(max(rcm_0)/num_rcm_t0_0)
		profs.append(prof)
		
	time_points = [data[time_frame]['mc_step'] for time_frame in range(num_time_frames)]
	print('time_points: ', time_points)
	profs = np.array(profs)
	return profs


def load_data(dir_edited_data, prefix_load_, num_time_frames):
	data = {}
	for time_frame in range(num_time_frames):
		# Load graph.
		prefix_load = '{}_{}'.format(time_frame, prefix_load_)
		data[time_frame] = utils.load(dir_edited_data, prefix_load, suffix_load)
	return data

if __name__ == '__main__':


	# Input file
	
	'''
	dir_target  = 'small_colony'
	prefixes    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(3) for id_f in range(7)]
	suffix_load  = 'connectivity_graph'
	time_frame  = 0
	ids_target_prefix  = [1,2,5,6]
	'''
	
	#'''
	dir_target  = 'small_colony2'
	prefixes = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7)]	
	suffix_load  = 'connectivity_graph'
	time_frame  = 0
	ids_target_prefix  = [7*5+2 ,7*4+2, 7*3+2 ,7*2+2]
	#'''

	
	#
	num_time_frames = 40
	#suffix_load = 'greedy_modularity'
	suffix_load = 'louvain'
	
	# Shared init
	dir_edited_data  = os.path.join('data4', dir_target)
	dir_imgs         = os.path.join('imgs4', dir_target, 'connectivity_matrix_dendrogram')
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	
	time = list(range(num_time_frames))
	num_rows    = 2
	num_columns = 2
	fig  = plt.figure(figsize=(8, 8))
	
	#for i, i_targ in enumerate([1,2,5,6]):
	for i, i_targ in enumerate(ids_target_prefix):

		prefix_load_ = prefixes[i_targ]
		data = load_data(dir_edited_data, prefix_load_, num_time_frames)
		profs = get_profiles( data )


		from networkx.algorithms.community import modularity
		
		modl = []
		for d in data.values():
			G   = d['multi_graph_CaMKII']
			rcm = d['rcm']
			modl.append( modularity(G, rcm) )
		# print('modl ', modl)
		modl_mean = np.mean(modl)
		modl_std  = np.std(modl)
		ax = fig.add_subplot( num_rows, num_columns, i+1 )
		ax.set_title(prefix_load_+', modularity: {:.3f} ± {:.3f}'.format(modl_mean, modl_std))
		
		ax.set_xlabel('Time (/10 frame)')
		ax.set_ylabel('Overalpped voxels')
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.set_xlim([np.min(time), np.max(time)])
		
		for prof in profs:
			ax.plot(time, prof, '-')
		profs
		ax.plot(time, np.average(profs, axis=0), 'k-', linewidth=3.0)

		ax.set_ylim(bottom=0)
		
	fig.savefig( os.path.join(dir_imgs, suffix_load + '.svg' ) )
	fig.savefig( os.path.join(dir_imgs, suffix_load + '.png' ), dpi=150 )
	plt.show()
	
	#total_num_CaMKII = len(utils.flatten(data[0]['rcm']))
	#print('total_num_CaMKII ', total_num_CaMKII )
	
