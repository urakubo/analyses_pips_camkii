
import os, glob, pickle, pprint, copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from skimage.measure import label

import lib.utils as utils
import lib.utils_panel as utils_panel
import lib.parameters as p
import lib.colormap as c
from specification_datasets import SpecDatasets

plt.rcParams.update( p.rc_param )


class PhaseDiagramConcDependence():
	def __init__( self ):
		
		pass
		
		
	def prepare_plot( self ):
		#self.fig  = plt.figure(figsize=(5, 5))
		self.fig  = plt.figure(figsize=(4, 4))
		self.fig.subplots_adjust(wspace=0.4,  hspace=0.6)
		ax = self.fig.add_subplot( 1, 1, 1 )
		ax.set_title( self.title )
		ax.set_xlabel('STG (beads / voxel) x 10-3')	
		ax.set_ylabel('GluN2B (beads / voxel) x 10-3')
		return ax
		
		
	def edit_data_save_them( self ):
		
		num_rows	= len( self.concs_glun2b )
		num_columns	= len( self.concs_stg )
		
		self.data = np.zeros([num_columns, num_rows], dtype = 'float')
		
		for i, stg in enumerate(self.stgs):
			for j, glun2b in enumerate(self.glun2bs):
				# Load data
				prefix = self.filename_edited_matrix(glun2b, stg)
				
				print('Target file: ', prefix)
				d              = utils.load(self.dir_edited_data, prefix, self.suffix)
				self.data[i,j] = self._modify_data(d)
		
		print('data ', self.data)
		prefix = self.basename
		suffix = 'data'
		utils.save( self.dir_edited_data, prefix, suffix, self.data )
		
		
	def load_data( self ):
		prefix = self.basename
		suffix = 'data'
		self.data = utils.load( self.dir_edited_data, prefix, suffix  )
		
		'''
		prefix = 'max_vol_All'
		suffix = 'data'
		volume = utils.load( self.dir_edited_data, prefix, suffix  )
		self.data = self.data / volume * 4e4
		print(self.data)
		'''
		
	def plot_data( self ):
		mx_min = 0
		my_min = 0
		mx_max = np.max(self.concs_stg[:-1]) * 1.05# None # 2.6
		my_max = np.max(self.concs_glun2b)   * 1.05 # None # 10.1
		
		ax = self.prepare_plot()
		cs, cb = utils_panel.plot_a_panel(ax, self.data, self.concs_stg, self.concs_glun2b, self.colormap, self.levels,\
			mx_min=mx_min, my_min=my_min, \
			mx_max=mx_max, my_max=my_max, \
			margin = 0
			)
		
		
	def save_plots( self ):
		
		dir_imgs = os.path.join(self.dir_imgs_root, 'phase_diagram')
		os.makedirs(dir_imgs, exist_ok=True)
		self.fig.savefig( os.path.join(dir_imgs, self.basename + '.svg' ) )
		self.fig.savefig( os.path.join(dir_imgs, self.basename + '.png' ) , dpi=150)
		plt.show()
		plt.clf()
		plt.close(fig=self.fig)
		
		
class PlotConnectivity():
	def __init__( self, species, type_analysis ):
		
		if species == 'CaMKII' and type_analysis == 'average':
			self.title    = 'Number of GluN2B bound to one CaMKII'
			self.basename = 'num_GluN2B_bound_to_one_CaMKII'
			self.colormap =  c.cmap_white_green_universal
			self.levels   = np.linspace(0,12,7)
		elif species == 'PSD95' and type_analysis == 'average_STG':
			self.title    = 'Number of STG bound to one PSD95'
			self.basename = 'num_STG_bound_to_one_PSD95'
			self.colormap = c.cmap_white_red_universal
			self.levels   = np.linspace(0,3,6)
		elif species == 'PSD95' and type_analysis == 'average_GluN2B':
			self.title    = 'Number of GluN2B bound to one PSD95'
			self.basename = 'num_GluN2B_bound_to_one_PSD95'
			self.colormap = c.cmap_white_purple_universal
			self.levels   = np.linspace(0,3,6)
		elif species == 'PSD95' and type_analysis == 'average_ratio':
			self.title    = 'Ratio of STG /(STG+GluN2B) bound to PSD95'
			self.basename = 'PSD95_bound_to_both_GluN2B_STG'
			self.colormap = plt.get_cmap('Greys') # plt.colormaps['Greys']
			self.levels   = np.linspace(0,1.0,11)
		elif species == 'PSD95' and type_analysis == 'ratio':
			self.title    = 'Ratio of PSD95 bound to both GluN2B and STG'
			self.basename = 'PSD95_bound_to_both_GluN2B_STG'
			self.colormap = plt.get_cmap('Greys') # plt.colormaps['Greys']
			self.levels   = np.linspace(0,1.0,11)
		elif species == 'PSD95' and type_analysis == 'ratio_condensate':
			self.title    = 'Ratio of PSD95 bound to both GluN2B and STG'
			self.basename = 'PSD95_bound_to_both_GluN2B_STG_condensate'
			self.colormap = plt.get_cmap('Greys') # plt.colormaps['Greys']
			self.levels   = np.linspace(0,1.0,11)
			#self.levels   = np.linspace(0,2.0,11)
		elif species == 'PSD95' and type_analysis == 'distribution':
			self.title    = 'PSD95 occupied with neigther GluN2B nor STG'
			self.basename = 'PSD95_occupied_with_neigther_GluN2B_nor_STG'
			self.colormap = plt.get_cmap('Greys') # plt.colormaps['Greys']
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
		elif self.species == 'PSD95' and self.type_analysis == 'average_STG':
			data      = d[self.species]['average']['STG_PSD95']
		elif self.species == 'PSD95' and self.type_analysis == 'average_GluN2B':
			data      = d[self.species]['average']['GluN2B_PSD95']
		elif self.species == 'PSD95' and self.type_analysis == 'average_ratio':
			average_STG    = d[self.species]['average']['STG_PSD95']
			average_GluN2B = d[self.species]['average']['GluN2B_PSD95']
			data = average_STG / (average_STG+average_GluN2B)
			
		elif self.species == 'PSD95' and 'ratio' in self.type_analysis:
			both   = d[self.species][self.type_analysis]['Both']
			STG    = d[self.species][self.type_analysis]['STG only']
			GluN2B = d[self.species][self.type_analysis]['GluN2B only']
			No_bound= d[self.species][self.type_analysis]['None']
			
			num_total = both + No_bound + STG + GluN2B
			data      = both / num_total
			
		elif self.species == 'PSD95' and 'distribution' in self.type_analysis:
			both0  = d[self.species]['distribution']['1 STG, 1 GluN2B']
			both1  = d[self.species]['distribution']['2 STG, 1 GluN2B']
			both2  = d[self.species]['distribution']['1 STG, 2 GluN2B']
			No_bound= d[self.species]['distribution']['None']
			STG1   = d[self.species]['distribution']['1 STG']
			STG2   = d[self.species]['distribution']['2 STG']
			STG3   = d[self.species]['distribution']['3 STG']
			GluN2B1= d[self.species]['distribution']['1 GluN2B']
			GluN2B2= d[self.species]['distribution']['2 GluN2B']
			GluN2B3= d[self.species]['distribution']['3 GluN2B']
			
			#num_total = both0 + both1 + both2 + No_bound + STG1 + STG2 + STG3 + GluN2B1 + GluN2B2 + GluN2B3
			#unoccupied = both0 + both1 + both2 + STG1 + STG2 + GluN2B1 + GluN2B2 
			#data      = unoccupied / num_total
			data = (both1+both2) / (both1 + both2 + STG3 + GluN2B3 )
		return data



class PlotCondVolume():
	def __init__( self, species ):
		
		if species == 'CaMKII':
			self.title    = 'Largest condensate volume of {}'.format(species)
			self.basename = 'volume_largest_cond_CaMKII'
			self.colormap =  c.cmap_white_green_universal
			self.levels   = np.linspace(0,7e4,8)
		elif species == 'STG':
			self.title    = 'Largest condensate volume of {}'.format(species)
			self.basename = 'volume_largest_cond_STG'
			self.colormap = c.cmap_white_red_universal
			self.levels   = np.linspace(0,3e4,7)
		elif species == 'All':
			self.title    = 'Largest condensate volume of {}'.format(species)
			self.basename = 'volume_largest_cond_All'
			self.colormap = c.cmap_white_red_universal
			self.levels   = np.linspace(0,10e4,7)
		else:
			raise ValueError("Not implemented, species: ", species)
		
		self.species  = species
		self.basename = 'max_vol_{}'.format( self.species )
		self.suffix   = 'sigma_2'
		
	def _modify_data(self, d):
		#
		data = d['region_condensate_in_grid_mesh'][self.species]
		labels, num_labels = label( data, return_num = True )
		vols_label = [np.sum( labels == i ) for i in range(1, num_labels+1)]
		data = np.max( vols_label )
		print( data )
		return data
		
		
class PlotPhaseDiagramConcDependence( PhaseDiagramConcDependence, SpecDatasets ):
	def __init__( self ):
		
		super().__init__()
		
		GluN2Bs = np.array([270,540,1080,2160,4320,6480,8640,12960,17280])
		GluN2Bs = GluN2Bs / np.prod(p.space_np)
		GluN2Bs = GluN2Bs * 1000
		GluN2Bs = np.flip(GluN2Bs)
		#GluN2Bs = np.log(GluN2Bs)
		self.GluN2Bs = GluN2Bs
		
		
		STGs = np.array([108,216,432,576,864,1728,2592,3456,4320,5184])
		STGs = STGs / np.prod(p.space_np)
		STGs = STGs * 1000
		#STGs = np.log(STGs)
		self.STGs = STGs
		
		
		self.num_rows = len( GluN2Bs )
		self.title    = 'Phase diagram'
		self.basename = 'phase_diagram_conc_dependence'
		# -1: Unclear
		# 1: Homogeneous LLPS (CaMKII)
		# 2: PIPS
		# 3: Homogeneous LLPS (STG)
		
		phase_diagram = [\
			[ 1, 1, 1, 1,   1, 1, 1, 1,   1, 1],
			[ 1, 1, 1, 1,   1, 1, 1, 1,   1, 1],
			[ 1, 1, 1, 1,   1, 1, 1, 1,   1, 1],
			[ 1, 1, 1, 1,   1, 1, 1, 1,   1, 1],
			[ 1, 1, 1, 1,   1, 1, 1, 1,   1, 1],
			[ 3, 3, 3, 3,   3, 3, 3, 3,   3, 3],
			[ 3, 3, 3, 3,   3, 3, 3, 3,   3, 3],
			[ 3, 3, 3, 3,   3, 3, 3, 3,   3, 3],
			[ 3, 3, 3, 3,   3, 3, 3, 3,   3, 3]]
			
		# 1: Partial engulfment
		phase_diagram_partial_engulfment = [\
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0],
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0],
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0],
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0],
			[ 0, 0, 0, 0,   0, 0, 1, 1,   1, 1],
			[ 0, 0, 0, 0,   0, 0, 1, 1,   1, 1],
			[ 0, 0, 0, 0,   0, 0, 0, 1,   1, 1],
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0],
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0]]
		phase_diagram_two_condensates = [\
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0],
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0],
			[ 0, 0, 0, 0,   1, 1, 1, 1,   1, 1],
			[ 0, 0, 1, 1,   1, 1, 1, 1,   1, 1],
			[ 0, 0, 1, 1,   1, 1, 1, 1,   1, 1],
			[ 0, 0, 0, 0,   0, 1, 1, 1,   1, 1],
			[ 0, 0, 0, 0,   0, 0, 0, 1,   1, 1],
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0],
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0]]
			
			
		self.phase_diagram = np.array(phase_diagram).T
		self.phase_diagram_partial_engulfment = np.array(phase_diagram_partial_engulfment).T
		self.phase_diagram_two_condensates = np.array(phase_diagram_two_condensates).T
		
	def plot( self ):
		
		''' Range for log plot
		mx_min = np.min(self.STGs) +np.log(0.9)
		my_min = np.min(self.GluN2Bs[:-1]) +np.log(0.9)
		mx_max = np.max(self.STGs) * 1.1# None # 2.6
		my_max = np.max(self.GluN2Bs) *1.1 # None # 10.1
		'''
		mx_min = 0
		my_min = 0
		mx_max = np.max(self.STGs[:-1]) * 1.05# None # 2.6
		my_max = np.max(self.GluN2Bs) *1.05 # None # 10.1
		
		colormap = c.cmap_phase_diagram7
		levels   = np.array([0.0,2.0,4.0])
		ax = self.prepare_plot()
		cs, cb = utils_panel.plot_a_panel(ax, self.phase_diagram, self.STGs, self.GluN2Bs, colormap, levels, \
			draw_border = True, \
			mx_min=mx_min, my_min=my_min, \
			mx_max=mx_max, my_max=my_max, \
			margin = 0.0)
		levels2		= [-1.5, 0.5, 1.5, 2.5]
		colormap2	= c.cmap_phase_diagram6
		utils_panel.plot_a_panel_overlay(ax, self.phase_diagram_two_condensates, \
			self.STGs, self.GluN2Bs,
			colormap2, levels2,\
			mx_min=mx_min, my_min=my_min, \
			mx_max=mx_max, my_max=my_max, \
			margin = 0.0)
		levels2		= [-0.5, 0.5, 1.5, 2.5]
		colormap2	= c.cmap_phase_diagram5
		utils_panel.plot_a_panel_overlay(ax, self.phase_diagram_partial_engulfment, \
			self.STGs, self.GluN2Bs,
			colormap2, levels2,\
			mx_min=mx_min, my_min=my_min, \
			mx_max=mx_max, my_max=my_max, \
			margin = 0.0)
		
		
class HandleConnectivityPhaseDiagramConcDependence(PlotConnectivity, PhaseDiagramConcDependence, SpecDatasets):
	def __init__( self, species, type_analysis ):
		PlotConnectivity.__init__(self, species, type_analysis )
		PhaseDiagramConcDependence.__init__(self)
		
class HandleCondVolumePhaseDiagramConcDependence(PlotCondVolume, PhaseDiagramConcDependence, SpecDatasets):
	def __init__( self, species ):
		PlotCondVolume.__init__(self, species )
		PhaseDiagramConcDependence.__init__(self)
		
if __name__ == '__main__':
	
	#'''
	pl = PlotPhaseDiagramConcDependence()
	pl.conc_dependence_merged()
	pl.plot()
	pl.save_plots()
	#'''
	
	'''
	species, type_analysis = 'CaMKII', 'average'
	#species, type_analysis = 'PSD95' , 'average'
	#species, type_analysis = 'PSD95' , 'ratio'
	# species, type_analysis = 'PSD95' , 'average_GluN2B'
	pl = HandleConnectivityPhaseDiagramConcDependence(species, type_analysis)
	pl.conc_dependence_merged()
	#pl.edit_data_save_them()
	pl.load_data()
	pl.plot_data()
	pl.save_plots()
	'''
	
	'''
	species = 'CaMKII' # 'CaMKII', 'STG'
	pl = HandleCondVolumePhaseDiagramConcDependence(species)
	pl.conc_dependence_merged()
	#pl.edit_data_save_them()
	pl.load_data()
	pl.plot_data()
	pl.save_plots()
	'''

