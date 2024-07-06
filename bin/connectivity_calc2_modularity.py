
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

from specification_datasets import SpecDatasets


def process( d ):
	
	
	
	# Make new graphs of CaMKII
	multi_graph_CaMKII, _, _, _ = \
		utils_graph.make_new_graphs_CaMKII_connectivity(d, nth_largest = 0)
	
	G = multi_graph_CaMKII
	
	'''
	from networkx.algorithms.community import greedy_modularity_communities
	rcm  = list(greedy_modularity_communities( G ))
	suffix_save='greedy_modularity'
	'''
	
	#'''
	from networkx.algorithms.community import louvain_communities
	rcm = list(louvain_communities( G ))
	suffix_save='louvain'
	#'''
	
	adj_matrices = [nx.adjacency_matrix(G, nodelist=r).todense() for r in rcm]
	densities = [np.sum(m) for m in adj_matrices]
	data = {}
	
	data['multi_graph_CaMKII']  = multi_graph_CaMKII
	data['clust_coeff']  = nx.average_clustering(nx.Graph(G)) # np.mean(list(nx.clustering(G).values()))
	data['modularity']   = modularity(G, rcm)
	data['adj_matrices'] = adj_matrices
	data['densities']    = densities
	data['areas']        = [m.shape[0]*m.shape[1] for m in adj_matrices]
	data['density']      = np.sum(densities) / G.number_of_nodes()
	data['tot_density']  = np.sum(nx.adjacency_matrix(G).todense()) / G.number_of_nodes()
	data['rcm']          = rcm
	data['suffix_save']  = suffix_save
	
	return data
	
	
	
class CalcModularityClustringC(SpecDatasets):
	def __init__( self ):
		
		self.suffix = 'connectivity_graph'
		
	def run( self ):
		rows	= self.valencies
		columns	= self.lengths
		
		self.num_rows		= len( rows )
		self.num_columns	= len( columns )
		print('num_column, num_row :', self.num_columns, self.num_rows)
		
		self.data = {}
		self.modulatiry        = np.zeros([self.num_rows, self.num_columns])
		self.density           = np.zeros([self.num_rows, self.num_columns])
		self.clust_coefficient = np.zeros([self.num_rows, self.num_columns])
		
		for i, column in enumerate( columns ):
			for j, row in enumerate( rows ):
				# Load data
				filename_edited_prefix = self.filename_edited_matrix(row, column)
				print('Target file: ', filename_edited_prefix, ', column: ', column, ', row: ', row)
				d = utils.load(self.dir_edited_data, \
							filename_edited_prefix, \
							self.suffix)
				dd = process( d )
				self.data[filename_edited_prefix] = dd
				self.modulatiry[j, i]        = dd['modularity']
				self.density[j, i]           = dd['density']
				self.clust_coefficient[j, i] = dd['clust_coeff']
				
	def save( self ):
		
		prefix = 'connectivity_CaMKII'
		suffix = 'matrix'
		utils.save(self.dir_edited_data, prefix, suffix, self.data)
		
		prefix = 'modularities'
		suffix = 'matrix'
		utils.save(self.dir_edited_data, prefix, suffix, self.modulatiry)
		
		prefix = 'densities'
		suffix = 'matrix'
		utils.save(self.dir_edited_data, prefix, suffix, self.density)
		
		prefix = 'average_clustering_coefficient'
		suffix = 'matrix'
		utils.save(self.dir_edited_data, prefix, suffix, self.clust_coefficient)
		
		
if __name__ == '__main__':
	
	pass
	
	
