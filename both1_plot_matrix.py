
import os, glob, pickle, pprint
import numpy as np
import matplotlib.pyplot as plt


import lib.utils as utils
import lib.parameters as p
import lib.colormap as c

plt.rcParams.update(p.rc_param)

class MatrixValencyLength():
	def __init__( self, valencies, lengths, dir_target ):
		plt.rcParams.update( {'font.size': 6} )
		
		# Input files
		self.dir_edited_data = os.path.join('data4', dir_target)
		self.dir_imgs        = os.path.join('imgs4', dir_target, 'matrix')
		os.makedirs(self.dir_imgs, exist_ok=True)
		
		self.valencies		= valencies
		self.lengths		= lengths
		
		self.num_rows		= len( valencies )
		self.num_columns	= len( lengths   )
		
		
	def run( self ):
		
		vals = {}
		# vals = np.zeros([self.num_rows, self.num_columns])
		self.fig  = plt.figure(figsize=(8, 8), tight_layout=True)
		#fig.subplots_adjust(wspace=0.4,  hspace=0.6)
		
		for i, v in enumerate( self.valencies ):
			for j, l in enumerate( self.lengths ):
				# Load data
				prefix = str(v).zfill(2)+'_'+str(l).zfill(3)
				row    = self.num_rows-i-1
				column = j+1
				print('Target file: ', prefix, ', column: ', column, ', row: ', row)
				d      = utils.load(self.dir_edited_data, prefix, self.suffix)
				
				legend = (row == 1) and (column == 1)
				title = prefix # prefix, None
				vv, _ = self.plot_a_graph(row, column, d, title, legend = legend)
				vals[prefix] = vv
		return vals
		
	def save( self ):
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.svg' ) )
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.png' ), dpi=150 )
		plt.show()
		plt.clf()
		plt.close(fig=self.fig)
		


class MatrixConcDependence():
	def __init__( self ):
		plt.rcParams.update( {'font.size': 6} )
		
		# Input files
		self.dir_edited_data 		= os.path.join('data4', 'conc_dependence_merged')
		self.filenames_edited_data 	= [str(i).zfill(3) for i in range(81) ]
		
		# Output files
		self.dir_imgs = os.path.join('imgs4', 'conc_dependence_merged','matrix')
		os.makedirs(self.dir_imgs, exist_ok=True)
		
		self.num_rows		= len( p.GluN2Bs )
		self.num_columns	= len( p.STGs )
		
		
	def run( self ):
		vals = np.zeros([self.num_rows, self.num_columns])
		self.fig  = plt.figure(figsize=(10, 10), tight_layout=True)
		#fig.subplots_adjust(wspace=0.4,  hspace=0.6)
		for i, stg in enumerate(p.STGs):
			for j, glun in enumerate(p.GluN2Bs):
				# Load data
				prefix = str(i).zfill(2)+'_'+str(j).zfill(2)
				
				title = prefix # prefix, None
				row    = self.num_rows-j-1
				column = i+1
				print('Target file: ', prefix, ', column: ', column, ', row: ', row)
				d = utils.load(self.dir_edited_data, prefix, self.suffix)
				vv, _ = self.plot_a_graph(row, column, d, title)
				vals[row, column-1] = vv
		return vals
		
		
	def save( self ):
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.svg' ) )
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.png' ), dpi=150 )
		plt.show()
		plt.clf()
		plt.close(fig=self.fig)
	
	
class PlotConnectivity():
	def __init__( self, species, type_analysis ):
		
		if species not in p.molecules_without_all.keys():
			raise ValueError("Erronous augment 'species': ", species)
		elif type_analysis not in ['average','distribution', 'CaMKII', 'PSD95', 'ratio']:
			raise ValueError("Erronous augment 'type_analysis': ", type_analysis)
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
		if  type_analysis in ['distribution', 'CaMKII', 'PSD95'] :
			ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
		return True, ax
		
		
class PlotConc():
	def __init__( self, target ):
		
		if target not in ['region_condensates', \
			'conc_CaMKII', 'conc_STG', 'conc_GluN2B', 'conc_PSD95', 'rdf', 'rdf_PSD95', 'concs_in_CaMKII', 'concs_in_STG',\
			'unshared_PSD95', 'shared_PSD95']:
			raise ValueError("Erronous augment 'target': ", target)
		else:
			self.target = target
		sigma = 2
		self.suffix = 'sigma_{}'.format(sigma)
		self.basename = 'matrix_{}'.format( self.target )
		
	def plot_a_graph( self, row, column, d, title, legend = None ):
		val = utils.select_plot(self.target, self.fig, self.num_rows, self.num_columns, row, column, d, title)
		return val, None
		
		
		
class PlotConcMatrixConcDependence(PlotConc, MatrixConcDependence):
	def __init__( self, target ):
		PlotConc.__init__(self, target)
		MatrixConcDependence.__init__(self)
	
class PlotConcMatrixValencyLength(PlotConc, MatrixValencyLength):
	def __init__( self, target ):
		
		valencies = range(2,14,2)
		lengths   = range(7)
		dir_target= 'valency_length'
		
		PlotConc.__init__(self, target )
		MatrixValencyLength.__init__(self, valencies, lengths, dir_target)
	
	
class PlotConnectivityMatrixConcDependence(PlotConnectivity, MatrixConcDependence):
	def __init__( self, species, type_analysis ):
		PlotConnectivity.__init__(self, species, type_analysis )
		MatrixConcDependence.__init__(self)
		
class PlotConnectivityMatrixValencyLength(PlotConnectivity, MatrixValencyLength):
	def __init__( self, species, type_analysis ):
		
		valencies = range(2,14,2)
		lengths   = range(7)
		dir_target= 'valency_length'
		
		PlotConnectivity.__init__(self, species, type_analysis )
		MatrixValencyLength.__init__(self, valencies, lengths, dir_target)
	
	
if __name__ == '__main__':
	
	
	# 'region_condensates', 'conc_CaMKII', 'conc_PSD95', 'conc_STG', 'conc_GluN2B', 'rdf',  'rdf_PSD95'
	# 'concs_in_CaMKII', 'concs_in_STG',
	# 'shared_PSD95', 'unshared_PSD95'
	#'''
	target = 'conc_CaMKII'
	plot_concs = PlotConcMatrixConcDependence(target) # PlotConcMatrixConcDependence, PlotConcMatrixValencyLength
	values = plot_concs.run()
	plot_concs.save()
	#'''

	'''
	species       = 'CaMKII' # 'STG','GluN2B', 'PSD95','CaMKII'
	type_analysis = 'distribution'
	# 'average' and 'distribution' for all,
	# species: 'GluN2B', type_analysis 'CaMKII' or 'PSD95'
	# species: 'PSD95' , type_analysis 'ratio'

	connectivity = PlotConnectivityMatrixValencyLength(species, type_analysis)
	values = connectivity.run()
	connectivity.save()	
	'''
	

