
import os, glob, pickle, pprint
import numpy as np
import utils
import parameters as p
import matplotlib.pyplot as plt

plt.rcParams.update(p.rc_param)
plt.rcParams.update( {'font.size': 6} )

class MatrixValencyLength():
	def __init__( self ):
		
		# Input files
		dir_target      = 'valency_length'
		self.dir_edited_data = os.path.join('data3', dir_target)
		self.dir_imgs        = os.path.join('imgs3', dir_target, 'matrix')
		os.makedirs(self.dir_imgs, exist_ok=True)
		
		self.valency = list(range(2,14,2)) 
		self.length  = [1, 2, 3, 4, 5, 6, 9]
		
		self.fnames_valency= { v: str(v).zfill(2) for v in self.valency }
		self.fnames_length = {ll: str(i).zfill(3) for i,ll in enumerate(self.length) }
		
		self.num_rows		= len( self.valency )
		self.num_columns	= len( self.length )
		
		
	def run( self ):
		
		vals = np.zeros([self.num_rows, self.num_columns])
		self.fig  = plt.figure(figsize=(8, 8), tight_layout=True)
		#fig.subplots_adjust(wspace=0.4,  hspace=0.6)
		
		for i, v in enumerate(self.valency):
			for j, ll in enumerate(self.length):
				# Load data
				prefix = self.fnames_valency[v]+'_'+self.fnames_length[ll]
				row    = self.num_rows-i-1
				column = j+1
				print('Target file: ', prefix, ', column: ', column, ', row: ', row)
				d      = utils.load(self.dir_edited_data, prefix, self.suffix)
				
				
				title = prefix # prefix, None
				vals[row, column-1] = self.plot_a_graph(row, column, d, title)
		
		
	def save( self ):
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.svg' ) )
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.png' ), dpi=150 )
		plt.show()
		plt.clf()
		plt.close(fig=self.fig)
		


class MatrixConcDependence():
	def __init__( self ):
		
		# Input files
		self.dir_edited_data 		= os.path.join('data3', 'conc_dependence')
		self.filenames_edited_data 	= [str(i).zfill(3) for i in range(48) ]
		
		# Output files
		self.dir_imgs = os.path.join('imgs3', 'conc_dependence','matrix')
		os.makedirs(self.dir_imgs, exist_ok=True)
		
		STG    = [540, 1620, 2160, 2700, 3240, 4520] 
		GluN2B = [570, 1080, 4320, 6480, 8640, 10800, 12960, 17280]
		volume = np.prod(p.space_np)
		self.STG    = [ s / volume for s in STG    ]
		self.GluN2B = [ n / volume for n in GluN2B ]
		
		self.num_rows		= len( GluN2B )
		self.num_columns	= len( STG )
		
		
	def run( self ):
		vals = np.zeros([self.num_rows, self.num_columns])
		self.fig  = plt.figure(figsize=(10, 10), tight_layout=True)
		#fig.subplots_adjust(wspace=0.4,  hspace=0.6)
		for i, stg in enumerate(self.STG):
			for j, glun in enumerate(self.GluN2B):
				# Load data
				id = i + j * len(self.STG)
				prefix = str(id).zfill(3)
				
				print('Target: {}'.format(prefix))
				
				title = prefix # prefix, None
				row    = self.num_rows-j-1
				column = i+1
				print('Target file: ', prefix, ', column: ', column, ', row: ', row)
				d = utils.load(self.dir_edited_data, prefix, self.suffix)
				vals[row, column-1] = self.plot_a_graph(row, column, d, title)
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
		
	def plot_a_graph( self, row, column, d, title ):
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
		return True
		
		
class PlotConc():
	def __init__( self, target ):
		
		if target not in ['region_condensates', 'conc_CaMKII', 'conc_STG', 'rdf', 'concs_in_CaMKII', 'concs_in_STG']:
			raise ValueError("Erronous augment 'target': ", target)
		else:
			self.target = target
		sigma = 2
		self.suffix = 'sigma_{}'.format(sigma)
		self.basename = 'matrix_{}'.format( self.target )
		
	def plot_a_graph( self, row, column, d, title ):
		val = utils.select_plot(self.target, self.fig, self.num_rows, self.num_columns, row, column, d, title)
		return val
		
		
		
class PlotConnectivityMatrixConcDependence(PlotConnectivity, MatrixConcDependence):
	def __init__( self, species, type_analysis ):
		PlotConnectivity.__init__(self, species, type_analysis )
		MatrixConcDependence.__init__(self)

class PlotConcMatrixConcDependence(PlotConc, MatrixConcDependence):
	def __init__( self, target ):
		PlotConc.__init__(self, target)
		MatrixConcDependence.__init__(self)
	
class PlotConnectivityMatrixValencyLength(PlotConnectivity, MatrixValencyLength):
	def __init__( self, species, type_analysis ):
		PlotConnectivity.__init__(self, species, type_analysis )
		MatrixValencyLength.__init__(self)

class PlotConcMatrixValencyLength(PlotConc, MatrixValencyLength):
	def __init__( self, target ):
		PlotConc.__init__(self, target )
		MatrixValencyLength.__init__(self)
	
	
if __name__ == '__main__':
	
	
	# 'region_condensates', 'conc_CaMKII', 'conc_STG', 'rdf', 'concs_in_CaMKII', 'concs_in_STG'
	'''
	target = 'conc_STG'
	plot_concs = PlotConcMatrixConcDependence(target)
	values = plot_concs.run()
	plot_concs.save()
	'''

	species       = 'GluN2B' # 'STG','GluN2B', 'PSD95','CaMKII'
	type_analysis = 'CaMKII'
	# 'average' and 'distribution' for all,
	# species: 'GluN2B', type_analysis 'CaMKII' or 'PSD95'
	# species: 'PSD95' , type_analysis 'ratio'
	#'''
	connectivity = PlotConnectivityMatrixValencyLength(species, type_analysis)
	values = connectivity.run()
	connectivity.save()	
	#'''
	

