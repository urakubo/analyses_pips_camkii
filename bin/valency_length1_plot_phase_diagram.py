
import os, glob, pickle, pprint, copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import mpl_toolkits.axes_grid1

from scipy.interpolate import griddata, RegularGridInterpolator

import lib.utils as utils
import lib.utils_panel as utils_panel
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph
from specification_datasets import SpecDatasets

plt.rcParams.update( p.rc_param )

def prepare_plot(title):
	#self.fig  = plt.figure(figsize=(5, 5))
	#self.fig  = plt.figure(figsize=(4, 4))
	fig  = plt.figure(figsize=(3.5, 3.5))
	fig.subplots_adjust(wspace=0.4,  hspace=0.6)
	ax = fig.add_subplot( 1, 1, 1 )
	ax.set_title(title)
	ax.set_xlabel('Linker length (l.u.)')
	ax.set_ylabel('Valency')
	return fig, ax


class PlotValencyLength( SpecDatasets ):
	def __init__( self ):
		
		pass
		
	def run( self ):
		
		self.num_rows		= len( self.valencies )
		self.num_columns	= len( self.lengths )
		
		data = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		for i, v in enumerate(self.valencies):
			for j, l in enumerate(self.lengths):
				# Load data
				prefix    = self.filename_edited_matrix(v, l)
				print('Target file: ', prefix)
				d         = utils.load(self.dir_edited_data, prefix, self.suffix)
				data[j,i] = self._modify_data(d)
		
		print('data ', data)
		self.fig, ax = prepare_plot(self.title)
		cs, cb = utils_panel.plot_a_panel(ax, data, self.real_lengths, self.valencies, self.colormap, self.levels, ticks = self.ticks, mx_min=2.0, my_min=4.0)
		# , mx_min=2.0, my_min=2.0
	def save( self ):
		
		dir_imgs = os.path.join(self.dir_imgs_root, 'phase_diagram')
		os.makedirs(dir_imgs, exist_ok=True)
		self.fig.savefig( os.path.join( dir_imgs, self.basename + '.svg' ) )
		self.fig.savefig( os.path.join( dir_imgs, self.basename + '.png' ) , dpi=150)
		plt.show()
		plt.clf()
		plt.close(fig=self.fig)
		

		
class PhaseDiagram():
	def __init__( self ):
		
		# self.dir_target = 'valency_length'
		
		self.p_lengths   = [1, 2, 3, 4, 5, 6, 9]
		self.p_valencies = [12, 10, 8, 6, 4, 2, -0] # [4, 6, 8, 10, 12] 
		self.basename = 'phase_diagram_valency_length'
		self.title    = 'Phase diagram'
		# 1: PIPS
		# 2: Partial engulfment
		# 3: iPIPS
		# 4: Homogeneous LLPS
		
		# 1, 2, 3, 4, 5, 6, 9
		self.phase_diagram = [\
			[ 1, 1, 1, 2,  2, 2, 2], # 12
			[ 1, 1, 1, 2,  2, 2, 2], # 10
			[ 1, 1, 1, 1,  2, 2, 2], # 8
			[ 1, 1, 1, 1,  2, 2, 2], # 6
			[ 1, 1, 1, 1,  2, 2, 2], # 4
			[ 1, 1, 1, 1,  2, 2, 2], # 2
			[ 1, 1, 1, 1,  2, 2, 2], # 0
			]
		
		self.STG_only = [\
			[ 0, 0, 0, 0,  0, 0, 0], # 12
			[ 0, 0, 0, 0,  0, 0, 0], # 10
			[ 0, 0, 0, 0,  0, 0, 0], # 8
			[ 0, 0, 0, 0,  0, 0, 0], # 6
			[ 1, 1, 1, 1,  1, 1, 1], # 4
			[ 2, 2, 2, 2,  2, 2, 2], # 2
			[ 2, 2, 2, 2,  2, 2, 2], # 0
			]
		
		self.phase_diagram = np.array(self.phase_diagram).T
		self.STG_only	= np.array(self.STG_only).T
		
		#self.phase_diagram = np.flipud(self.phase_diagram)
		#self.STG_only      = np.flipud(self.STG_only)
		
		
	def plot( self ):
		
		levels1		= np.array([0.5, 1.5, 2.5, 3.5, 4.5])
		colormap1	= c.cmap_phase_diagram2
		
		levels2		= [-0.5, 0.5, 1.5, 2.5]
		colormap2	= c.cmap_phase_diagram3
		
		max_conc      = np.max(self.phase_diagram)
		
		print('phase_diagram.shape ', self.phase_diagram.shape)
		print('STG_only.shape ',  self.STG_only.shape)
		print('self.lengths ', self.p_lengths)
		print('self.valencies ', self.p_valencies)
		
		
		self.fig, ax = prepare_plot(self.title)
		cs, cb = utils_panel.plot_a_panel(ax, self.phase_diagram, self.p_lengths, self.p_valencies, colormap1, levels1, draw_border = True, mx_min = 2.0, my_min = 2.0)
		utils_panel.plot_a_panel_overlay(ax, self.STG_only, self.p_lengths, self.p_valencies, colormap2, levels2, mx_min = 2.0, my_min = 2.0)
		
		
class Connectivity():
	def __init__( self, species, type_analysis ):
		
		#self.dir_target = 'valency_length'
		#self.lengths    = p.lengths
		#self.valencies  = p.valencies
		
		if species == 'CaMKII' and type_analysis == 'average':
			self.title    = 'Number of GluN2B bound to one CaMKII'
			self.basename = 'num_GluN2B_bound_to_one_CaMKII'
			self.colormap =  c.cmap_white_green_universal
			self.levels   = np.linspace(0,12,7)
			self.ticks = [0, 4, 8, 12]
			self.suffix = 'connectivity_graph'
		elif species == 'PSD95' and type_analysis == 'average':
			self.title    = 'Number of STG bound to one PSD95'
			self.basename = 'num_STG_bound_to_one_PSD95'
			self.colormap = c.cmap_white_red_universal
			self.levels   = np.linspace(0,3,6)
			self.ticks = [0, 1, 2, 3]
			self.suffix = 'connectivity_graph'
		elif species == 'PSD95' and type_analysis == 'ratio':
			self.title    = 'Ratio of PSD95 bound to both GluN2B and STG'
			self.basename = 'PSD95_bound_to_both_GluN2B_STG'
			self.colormap =  plt.colormaps['Greys']
			self.levels   = np.linspace(0,1.0,11)
			self.ticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
			self.suffix = 'connectivity_graph'
		elif species == 'CaMKII' and type_analysis == 'conc_in_CaMKII_condensate':
			self.title    = 'conc CaMKII in CaMKII condensate'
			self.basename = 'conc_CaMKII_in_CaMKII_condensate'
			self.colormap =  c.cmap_white_green_universal
			self.levels   = np.linspace(0,1.0,11)
			self.ticks = [0, 0.5, 1.0]
			self.suffix = 'sigma_2'
		elif species == 'GluN2B' and type_analysis == 'conc_in_CaMKII_condensate':
			self.title    = 'conc GluN2B in CaMKII condensate'
			self.basename = 'conc_GluN2B_in_CaMKII_condensate'
			self.colormap =  c.cmap_white_purple_universal
			self.levels   = np.linspace(0,1.0,11)
			self.ticks = [0, 0.5, 1.0]
			self.suffix = 'sigma_2'
		elif species == 'All' and type_analysis == 'conc_in_CaMKII_condensate':
			self.title    = 'conc Total in CaMKII condensate'
			self.basename = 'conc_total_in_CaMKII_condensate'
			self.colormap =  plt.colormaps['Greys']
			self.levels   = np.linspace(0,1.0,11)
			self.ticks = [0, 0.5, 1.0]
			self.suffix = 'sigma_2'
		
		
		else:
			raise ValueError("Not implemented, species: ", species, ", type_analysis: ", type_analysis)
		
		self.species = species
		self.type_analysis = type_analysis
		self.basename = '{}_{}'.format( self.species, self.type_analysis )
	
	
	def _modify_data(self, d):
		#print('d ')
		#print(d[self.species][self.type_analysis])
		if self.species == 'CaMKII' and self.type_analysis == 'average':
			data      = d[self.species][self.type_analysis]['GluN2B']
		elif self.species == 'PSD95' and self.type_analysis == 'average':
			data      = d[self.species][self.type_analysis]['STG_PSD95']
		elif self.species == 'PSD95' and self.type_analysis == 'ratio':
			num_total = sum( d[self.species][self.type_analysis].values() )
			data      = d[self.species][self.type_analysis]['Both'] / num_total
		elif self.species == 'CaMKII' and self.type_analysis == 'conc_in_CaMKII_condensate':
			data      = d['conc_condensate']['CaMKII']['CaMKII']
		elif self.type_analysis == 'conc_in_CaMKII_condensate':
			data      = d['conc_condensate']['CaMKII'][self.species]
			
		return data


		
class RelaxzationTime():
	def __init__( self ):
		
		#self.dir_target = 'small_colony2'
		#self.lengths    = p.lengths2
		#self.valencies  = p.valencies2
		
		self.title    = 'Relaxzation time'
		self.basename = 'relaxzation_time_for_mixture'
		self.colormap = plt.colormaps['Greys']
		self.levels   = np.linspace(-1,1,30)
		self.ticks_level = [-1,0,1]
		
		self.prefix_loadname = 'relaxzation_time'
		
	def plot_mixture( self , no_log = False):
		#
		
		valencies    = self.valencies
		lengths      = self.lengths
		real_lengths = self.real_lengths
		num_columns = len( lengths )
		num_rows    = len( valencies )
		
		suffix = 'matrix'
		self.data    = utils.load(self.dir_edited_data, self.prefix_loadname, suffix).T
		if no_log == False:
			self.data    = np.log10( self.data )
		#print('OK!')
		self.fig, ax = prepare_plot(self.title)
		cs, cb = utils_panel.plot_a_panel_log(ax, self.data, real_lengths, valencies, self.colormap, self.levels,\
			ticks_level = self.ticks_level)
		
		print(self.data)
		
		'''
		data = np.zeros([num_columns, num_rows], dtype = 'float')
		for i, v in enumerate( valencies ):
			for j, l in enumerate( lengths ):
				# Load data
				prefix    = self.filename_edited_matrix(v,l)
				print('Target file: ', prefix)
				data[j,i] = np.log10( np.median(d[prefix]) )
		print('data ', data)
		self.fig, ax = prepare_plot(self.title)
		cs, cb = utils_panel.plot_a_panel_log(ax, data, real_lengths, valencies, self.colormap, self.levels)
		'''
		
		
class ModularityDensityClustering():
	def __init__( self, property ):
		
		#self.dir_target     = 'small_colony2'
		#self.lengths   = p.lengths2
		#self.valencies = p.valencies2
		self.log = False
		if property == 'modularity':
			self.title    = 'Modularity'
			self.basename = 'modularity'
			self.prefix = 'modularities'
			self.colormap = plt.colormaps['Greys']
			self.levels   = np.linspace(0,1,10)
			self.ticks = [0, 0.2, 0.4, 0.6, 0.8, 1]
		elif property == 'density':
			self.title    = 'Density'
			self.prefix = 'densities'
			self.basename = 'density'
			self.colormap = plt.colormaps['Greys']
			self.levels   = np.linspace(0,8,30)		
		elif property == 'clustering':
			self.title    = 'Clustering coefficient'
			self.prefix   = 'average_clustering_coefficient'
			self.basename = 'clustering'
			self.colormap = plt.colormaps['Greys']
			self.levels   = np.linspace(0,0.20,10)
			self.ticks    = [0, 0.05, 0.1, 0.15, 0.2]
		elif property == 'clustering_log':
			self.title    = 'Clustering coefficient'
			self.prefix   = 'average_clustering_coefficient'
			self.basename = 'clustering_log'
			self.colormap = plt.colormaps['Greys']
			self.levels   = np.linspace(-3,0,10)
			self.ticks    = [-3, -2, -1, 0]
			self.log      = True
		elif property == 'FRAP':
			self.title    = 'CaMKII FRAP'
			self.prefix   = 'FRAP_merged'
			self.basename = 'FRAP'
			self.colormap = plt.colormaps['Greys']
			self.levels   = np.linspace(-2,2.5,10)
			self.ticks = [-2,-1, 0, 1, 2]
		elif property == 'FRAP_GluN2B':
			# self.dir_target     = 'small_colony3'
			self.title    = 'GluN2B FRAP'
			self.prefix   = 'FRAP_GluN2B'
			self.basename = 'FRAP_GluN2B'
			self.colormap = plt.colormaps['Greys']
			self.levels   = np.linspace(0, 0.02, 10)
			self.ticks = [0, 0.01, 0.02]
		else:
			sys.exit('Target unspecified.')
		#
	def plot( self ):
		#
		
		valencies    = self.valencies
		lengths      = self.lengths
		real_lengths = self.real_lengths
		num_columns = len( lengths )
		num_rows    = len( valencies )
		
		#
		suffix = 'matrix'
		self.d = utils.load(self.dir_edited_data, self.prefix, suffix)
		#
		data = np.zeros([num_columns, num_rows], dtype = 'float')
		for i, v in enumerate(valencies):
			for j, l in enumerate(lengths):
				# Load data
				prefix    = self.filename_edited_matrix(v,l)
				print('Target file: ', prefix)
				data[j,i] = self.d[prefix]
		
		print('data ', data)
		self.fig, ax = prepare_plot( self.title )
		cs, cb = utils_panel.plot_a_panel_log(ax, data, \
			lengths, \
			valencies, \
			self.colormap, \
			self.levels, \
			ticks_level = self.ticks)
		#
	def plot2( self ):
		#
		valencies    = self.valencies
		lengths      = self.lengths
		real_lengths = self.real_lengths
		#
		suffix = 'matrix'
		self.data = utils.load(self.dir_edited_data, self.prefix, suffix)
		self.data = self.data.T
		print('data ', self.data)
		if self.log == True:
			self.data = np.log10(self.data)
		print('data ', self.data)
		self.fig, ax = prepare_plot( self.title )
		cs, cb = utils_panel.plot_a_panel_log(ax, self.data, \
			real_lengths, \
			valencies, \
			self.colormap, \
			self.levels, \
			ticks_level = self.ticks)		
		
		
class PlotConnectivityValencyLength(Connectivity, PlotValencyLength):
	def __init__( self, species, type_analysis ):
		Connectivity.__init__(self, species, type_analysis)
		PlotValencyLength.__init__(self )
		
class PlotPhaseDiagramValencyLength(PhaseDiagram, PlotValencyLength):
	def __init__( self ):
		PhaseDiagram.__init__(self)
		PlotValencyLength.__init__(self )
		
		
class PlotRelaxzationTimeValencyLength(RelaxzationTime, PlotValencyLength):
	def __init__( self ):
		RelaxzationTime.__init__(self)
		PlotValencyLength.__init__(self )
		

class PlotPropertiesValencyLength(ModularityDensityClustering, PlotValencyLength):
	def __init__( self, property ):
		ModularityDensityClustering.__init__(self, property)
		PlotValencyLength.__init__(self )


if __name__ == '__main__':
	
	
	'''
	pl = PlotPhaseDiagramValencyLength()
	pl.valency_length() # Only for the save directory?
	pl.plot()
	pl.save()
	
	
	#species, type_analysis = 'CaMKII', 'average'
	#species, type_analysis = 'PSD95' , 'average'
	species, type_analysis = 'PSD95' , 'ratio'
	pl = PlotConnectivityValencyLength(species, type_analysis)
	pl.valency_length()
	pl.run()
	pl.save()
	
	pl = PlotRelaxzationTimeValencyLength()
	pl.valency_length_small_colony2()
	pl.run_mixture()
	pl.save()
	'''	
	
	#'''	
	property = 'modularity' # 'density', 'modularity', 'clustering', 'FRAP'
	pl = PlotPropertiesValencyLength(property)
	pl.valency_length_small_colony2()
	pl.plot()
	pl.save()
	#'''
	
	'''
	property = 'FRAP_GluN2B' # 'density', 'modularity', 'clustering', 'FRAP', 'FRAP_GluN2B'
	pl = PlotPropertiesValencyLength(property)
	pl.valency_length_small_colony3()
	pl.plot()
	pl.save()
	'''
	