
import os, glob, pickle, pprint
import numpy as np
import matplotlib.pyplot as plt


import lib.utils as utils
import lib.parameters as p
import lib.colormap as c

plt.rcParams.update(p.rc_param)

from specification_datasets import SpecDatasets
from lib.utils_matrix import Matrix
	
	
class PlotConnectivity():
	def __init__( self, species, type_analysis ):
		
		if species not in p.molecules_without_all.keys():
			raise ValueError("Erronous argument 'species': ", species)
		elif type_analysis not in ['average','distribution', 'CaMKII', 'PSD95', 'ratio']:
			raise ValueError("Erronous argument 'type_analysis': ", type_analysis)
		else:
			self.species = species
			self.type_analysis = type_analysis
		
		self.suffix = 'connectivity_graph'
		self.ymax = {'STG'   : {'average': 4, 'distribution':4000},\
				'GluN2B': {'average': 2, 'distribution':8000}, \
				'PSD95' : {'average': 3, 'distribution':3000},
				'CaMKII': {'average':12, 'distribution':1000}}
		self.ymax['GluN2B']['CaMKII'] =  10000
		self.ymax['GluN2B']['PSD95']  =  10000
		self.ymax['PSD95']['ratio']   =  10000
		
		self.basename = '{}_{}'.format( self.species, self.type_analysis )
		
	def plot_a_graph( self, row, column, d, title, legend = None ):
		dist = d[self.species][self.type_analysis]
		ax = self.fig.add_subplot( self.num_rows, self.num_columns, row*self.num_columns+column )
		ax.set_title(title)
		ax.bar(*zip(*dist.items()), width=0.5)
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.set_ylim(0,self.ymax[self.species][self.type_analysis])
		ax.set_ylabel('(Number)')
		if  self.type_analysis in ['distribution', 'CaMKII', 'PSD95'] :
			ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
		return True, ax
		
		
class PlotConc():
	def __init__( self, target ):
		
		if target not in ['region_condensates', \
			'conc_CaMKII', 'conc_STG', 'conc_GluN2B', 'conc_PSD95', 'rdf', 'rdf_CG', 'rdf_PSD95', \
			'concs_in_CaMKII', 'concs_in_STG',\
			'unshared_PSD95', 'shared_PSD95','conc_unrotated_CaMKII']:
			raise ValueError("Erronous argument 'target': ", target)
		else:
			self.target = target
		sigma = 2
		self.suffix = 'sigma_{}'.format(sigma)
		self.basename = 'matrix_{}'.format( self.target )
		
	def plot_a_graph( self, row, column, d, title, legend = None ):
		val = utils.select_plot(self.target, self.fig, self.num_rows, self.num_columns, row, column, d, title)
		return val, None
	
	
class PlotConcMatrix(PlotConc, Matrix, SpecDatasets):
	def __init__( self, target ):
		PlotConc.__init__(self, target)
		Matrix.__init__(self)
	
	
class PlotConnectivityMatrix(PlotConnectivity, Matrix, SpecDatasets):
	def __init__( self, species, type_analysis ):
		PlotConnectivity.__init__(self, species, type_analysis )
		Matrix.__init__(self)
	
	
if __name__ == '__main__':
	
	
	# 'region_condensates', 'conc_CaMKII', 'conc_PSD95', 'conc_STG', 'conc_GluN2B', 'rdf',  'rdf_PSD95'
	# 'concs_in_CaMKII', 'concs_in_STG',
	# 'shared_PSD95', 'unshared_PSD95', 'conc_unrotated_CaMKII'
	#'''
	target = 'conc_unrotated_CaMKII' # 'conc_CaMKII'
	plot_concs = PlotConcMatrix(target) # PlotConcMatrixConcDependence, PlotConcMatrixValencyLength
	values = plot_concs.run()
	plot_concs.save()
	#'''

	'''
	species       = 'CaMKII' # 'STG','GluN2B', 'PSD95','CaMKII'
	type_analysis = 'distribution'
	# 'average' and 'distribution' for all,
	# species: 'GluN2B', type_analysis 'CaMKII' or 'PSD95'
	# species: 'PSD95' , type_analysis 'ratio'

	connectivity = PlotConnectivityMatrix(species, type_analysis)
	values = connectivity.run()
	connectivity.save()	
	'''
	

