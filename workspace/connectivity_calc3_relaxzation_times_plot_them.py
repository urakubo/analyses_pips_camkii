
# https://stackoverflow.com/questions/59821151/plot-the-dendrogram-of-communities-found-by-networkx-girvan-newman-algorithm
# https://qiita.com/s-wakaba/items/a93f03f27137cff4a26c

import os, sys, glob, pprint,itertools
import numpy as np


# from itertools import chain, combinations
import matplotlib.pyplot as plt

#import networkx as nx
#from scipy.cluster.hierarchy import dendrogram
#from networkx.algorithms.community.centrality import girvan_newman
from scipy.optimize import curve_fit

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph

plt.rcParams.update(p.rc_param)




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
	ids_target_prefix  = [0,1,2,3,4,5,6]
	nums_time_frames = [40]*len(ids_target_prefix)
	nums_time_frames[-1] = 80
	'''
	
	#'''
	dir_target  = 'small_colony2'
	prefixes = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7)]	
	suffix_load  = 'connectivity_graph'
	time_frame  = 0
	#ids_target_prefix  = [7*5+0, 7*5+1, 7*5+2, 7*5+3, 7*5+4, 7*5+5] # Valency 12
	#ids_target_prefix  = [7*4, 7*4+1, 7*4+2, 7*4+3, 7*4+4, 7*4+5] # Valency 10
	ids_target_prefix  = [7*3, 7*3+1, 7*3+2, 7*3+3, 7*3+4, 7*3+5] # Valency 8
	#ids_target_prefix  = [7*2, 7*2+1, 7*2+2, 7*2+3, 7*2+4, 7*2+5] # Valency 6
	#ids_target_prefix  = [7*1, 7*1+1, 7*1+2, 7*1+3, 7*1+4, 7*1+5] # Valency 4
	nums_time_frames = [80]*len(ids_target_prefix)
	#'''
	
	
	#
	#suffix_load = 'greedy_modularity'
	suffix_load = 'louvain'
	
	# Shared init
	dir_edited_data  = os.path.join('data4', dir_target)
	dir_imgs         = os.path.join('imgs4', dir_target, 'connectivity_relaxzation_time')
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	num_rows    = 2
	num_columns = len( ids_target_prefix )
	fig  = plt.figure(figsize=(18, 6),tight_layout=True)
	
	for i, (i_targ, num_time_frames) in enumerate(zip(ids_target_prefix, nums_time_frames)):
	
		# Load data
		prefix_load_ = prefixes[i_targ]
		data = load_data(dir_edited_data, prefix_load_, num_time_frames)
		
		time_points, profs, taus, params, colors = utils_graph.get_profiles(data)
		
		
		# Plot figure
		ax = fig.add_subplot( num_rows, num_columns, i+1 )
		utils_graph.prep_fig(ax, time_points)
		for col, prof, param in zip( colors, profs, params ):
			ax.plot(time_points, prof, 'o', markersize=3, color= col )
			ax.plot(time_points, utils_graph.func_exponential(time_points,param[0],param[1],param[2]), '-', color= col )
		ax.set_title(prefix_load_+', Tau: {:.4f}'.format(np.median(taus)))
		
		ax = fig.add_subplot( num_rows, num_columns, num_columns+i+1 )
		ax.set_title('Tau (median) {:.4f}'.format(np.median(taus)))
		ax.plot([0]*len(taus), taus, 'o', markersize=3, color= col )

	fig.savefig( os.path.join(dir_imgs, suffix_load + '.svg' ) )
	fig.savefig( os.path.join(dir_imgs, suffix_load + '.png' ), dpi=150 )
	plt.show()
	
	
	

