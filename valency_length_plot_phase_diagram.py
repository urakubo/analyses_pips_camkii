
import os, glob, pickle, pprint, copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import mpl_toolkits.axes_grid1

from scipy.interpolate import griddata, RegularGridInterpolator

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph

plt.rcParams.update( p.rc_param )


class PlotValencyLength():
	def __init__( self ):
		
		# Parameters
		self.sigma = 2
		
		self.dir_edited_data = os.path.join('data4', self.dir_target)
		self.dir_imgs        = os.path.join('imgs4', self.dir_target, 'phase_diagram')
		os.makedirs(self.dir_imgs, exist_ok=True)
		self.num_rows		= len( self.valencies )
		self.num_columns	= len( self.lengths )
		
		
	def run( self ):
		#
		data = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		for i, v in enumerate(self.valencies):
			for j, l in enumerate(self.lengths):
				
				# Load data
				prefix    = p.fnames_valency[v]+'_'+p.fnames_length[l]
				print('Target file: ', prefix)
				d         = utils.load(self.dir_edited_data, prefix, self.suffix)
				data[j,i] = self._modify_data(d)
		
		print('data ', data)
		ax = self.prepare_plot()
		cs, cb = utils.plot_a_panel(ax, data, p.lengths, p.valencies, self.colormap, self.levels)
		
	def save( self ):
		
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.svg' ) )
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.png' ) , dpi=150)
		plt.show()
		plt.clf()
		plt.close(fig=self.fig)
		
	def prepare_plot( self ):
		#self.fig  = plt.figure(figsize=(5, 5))
		#self.fig  = plt.figure(figsize=(4, 4))
		self.fig  = plt.figure(figsize=(3.5, 3.5))
		self.fig.subplots_adjust(wspace=0.4,  hspace=0.6)
		ax = self.fig.add_subplot( 1, 1, 1 )
		ax.set_title(self.title)
		ax.set_xlabel('Linker length (l.u.)')
		ax.set_ylabel('Valency')
		return ax
		
class PhaseDiagram():
	def __init__( self ):
		
		self.dir_target = 'valency_length'
		
		self.valencies   = [12, 10, 8, 6, 4, 2, -2] # [4, 6, 8, 10, 12] 
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
			[ 1, 1, 1, 2,  2, 2, 2], # 8
			[ 1, 1, 1, 2,  2, 2, 2], # 6
			[ 1, 1, 1, 2,  2, 2, 2], # 4
			[ 1, 1, 1, 2,  2, 2, 2], # 2
			[ 1, 1, 1, 2,  2, 2, 2], # -2
			]
		
		self.STG_only = [\
			[ 0, 0, 0, 0,  0, 0, 0], # 12
			[ 0, 0, 0, 0,  0, 0, 0], # 10
			[ 0, 0, 0, 0,  0, 0, 0], # 8
			[ 0, 0, 0, 0,  0, 0, 0], # 6
			[ 0, 0, 0, 0,  0, 0, 0], # 4
			[ 1, 1, 1, 1,  1, 1, 1], # 2
			[ 1, 1, 1, 1,  1, 1, 1], # -2
			]
		
	def plot( self ):
		
		phase_diagram = np.array(self.phase_diagram).T
		levels1		= np.array([0.5, 1.5, 2.5, 3.5, 4.5])
		colormap1	= c.cmap_phase_diagram2
		
		STG_only	= np.array(self.STG_only).T
		levels2		= [-0.5, 0.5, 1.5]
		colormap2	= c.cmap_phase_diagram3
		
		max_conc      = np.max(phase_diagram)
		
		ax = self.prepare_plot()
		
		cs, cb = utils.plot_a_panel(ax, phase_diagram, p.lengths, self.valencies, colormap1, levels1, draw_border = True)
		utils.plot_a_panel_overlay(ax, STG_only, p.lengths, self.valencies, colormap2, levels2)
		
		
class Connectivity():
	def __init__( self, species, type_analysis ):
		
		self.dir_target = 'valency_length'
		self.lengths   = p.lengths
		self.valencies = p.valencies
		
		if species == 'CaMKII' and type_analysis == 'average':
			self.title    = 'Number of GluN2B bound to one CaMKII'
			self.basename = 'num_GluN2B_bound_to_one_CaMKII'
			self.colormap =  c.cmap_white_green_universal
			self.levels   = np.linspace(0,12,7)
		elif species == 'PSD95' and type_analysis == 'average':
			self.title    = 'Number of STG bound to one PSD95'
			self.basename = 'num_STG_bound_to_one_PSD95'
			self.colormap = c.cmap_white_red_universal
			self.levels   = np.linspace(0,3,6)
		elif species == 'PSD95' and type_analysis == 'ratio':
			self.title    = 'Ratio of PSD95 bound to both GluN2B and STG'
			self.basename = 'PSD95_bound_to_both_GluN2B_STG'
			self.colormap =  plt.colormaps['Greys']
			self.levels   = np.linspace(0,1.0,11)		
		else:
			raise ValueError("Not implemented, species: ", species, ", type_analysis: ", type_analysis)
		
		self.species = species
		self.type_analysis = type_analysis
		self.basename = '{}_{}'.format( self.species, self.type_analysis )
		self.suffix = 'connectivity_graph'
	
	
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
		return data


		
class RelaxzationTimeForMixture():
	def __init__( self ):
		
		self.dir_target     = 'small_colony2'
		self.lengths   = p.lengths2
		self.valencies = p.valencies2
		
		self.title    = 'Relaxzation time'
		self.basename = 'relaxzation_time_for_mixture'
		self.colormap = plt.colormaps['Greys']
		self.levels   = np.linspace(-1,1,30)
		
		
	def run_mixture( self ):
		#
		prefix = 'relaxzation_times'
		suffix = 'matrix'
		d    = utils.load(self.dir_edited_data, prefix, suffix)
		data = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		for i, v in enumerate(self.valencies):
			for j, l in enumerate(self.lengths):
				# Load data
				prefix    = p.fnames_valency2[v]+'_'+p.fnames_length2[l]
				print('Target file: ', prefix)
				data[j,i] = np.log10( np.median(d[prefix]) )
				if prefix in ['12_000', '10_000', '04_000', '04_001', '04_006']:
					data[j,i] = np.nan
		
		print('data ', data)
		ax = self.prepare_plot()
		cs, cb = utils.plot_a_panel_log(ax, data, self.lengths, self.valencies, self.colormap, self.levels)
		
		
class ModularityDensityClustering():
	def __init__( self, property ):
		
		self.dir_target     = 'small_colony2'
		self.lengths   = p.lengths2
		self.valencies = p.valencies2
		
		if property == 'modularity':
			self.title    = 'Modularity'
			self.basename = 'modularity'
			self.prefix = 'modularities'
			self.colormap = plt.colormaps['Greys']
			self.levels   = np.linspace(0,1,10)
		elif property == 'density':
			self.title    = 'Density'
			self.prefix = 'densities'
			self.basename = 'density'
			self.colormap = plt.colormaps['Greys']
			self.levels   = np.linspace(0,8,30)		
		elif property == 'clustering':
			self.title    = 'Clustering coefficient'
			self.prefix = 'average_clustering_coefficient'
			self.basename = 'clustering'
			self.colormap = plt.colormaps['Greys']
			self.levels   = np.linspace(0,0.30,10)
		elif property == 'FRAP':
			self.title    = 'CaMKII FRAP'
			self.prefix   = 'FRAP'
			self.basename = 'FRAP'
			self.colormap = plt.colormaps['Greys']
			self.levels   = np.linspace(-1,2.5,10)
			
		
	def plot( self ):
		#
		suffix = 'matrix'
		d    = utils.load(self.dir_edited_data, self.prefix, suffix)
		data = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		for i, v in enumerate(self.valencies):
			for j, l in enumerate(self.lengths):
				# Load data
				prefix    = p.fnames_valency[v]+'_'+p.fnames_length[l]
				print('Target file: ', prefix)
				data[j,i] = d[prefix]
				if prefix in ['04_000', '04_001', '04_006']:
					data[j,i] = np.nan
		
		print('data ', data)
		ax = self.prepare_plot()
		cs, cb = utils.plot_a_panel_log(ax, data, self.lengths, self.valencies, self.colormap, self.levels)
		
		
class PlotConnectivityValencyLength(Connectivity, PlotValencyLength):
	def __init__( self, species, type_analysis ):
		Connectivity.__init__(self, species, type_analysis)
		PlotValencyLength.__init__(self )
		
class PlotPhaseDiagramValencyLength(PhaseDiagram, PlotValencyLength):
	def __init__( self ):
		PhaseDiagram.__init__(self)
		PlotValencyLength.__init__(self )
		
		
class PlotRelaxzationTimeForMixtureValencyLength(RelaxzationTimeForMixture, PlotValencyLength):
	def __init__( self ):
		RelaxzationTimeForMixture.__init__(self)
		PlotValencyLength.__init__(self )
		

class PlotPropertiesValencyLength(ModularityDensityClustering, PlotValencyLength):
	def __init__( self, property ):
		ModularityDensityClustering.__init__(self, property)
		PlotValencyLength.__init__(self )


if __name__ == '__main__':
	
	'''
	pl = PlotPhaseDiagramValencyLength()
	pl.plot()
	pl.save()
	'''
	
	'''
	species, type_analysis = 'CaMKII', 'average'
	#species, type_analysis = 'PSD95' , 'average'
	#species, type_analysis = 'PSD95' , 'ratio'
	
	pl = PlotConnectivityValencyLength(species, type_analysis)
	pl.run()
	pl.save()
	'''
	'''	
	pl = PlotRelaxzationTimeForMixtureValencyLength()
	pl.run_mixture()
	pl.save()
	'''	
	

	property =  'FRAP' # 'density', 'modularity', 'clustering', 'FRAP'
	pl = PlotPropertiesValencyLength(property)
	pl.plot()
	pl.save()


