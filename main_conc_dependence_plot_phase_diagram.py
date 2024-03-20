
import os, glob, pickle, pprint, copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# import mpl_toolkits.axes_grid1
# from scipy.interpolate import griddata, RegularGridInterpolator

import utils
import colormap as c
import parameters as p

from skimage.measure import label

plt.rcParams.update( p.rc_param )


class PhaseDiagramConcDependence():
	def __init__( self ):
		
		# Parameters
		self.sigma = 2
		dir_target = 'conc_dependence'
		self.dir_edited_data = os.path.join('data3',dir_target)
		self.dir_imgs        = os.path.join('imgs3', dir_target,'phase_diagram')
		
		os.makedirs(self.dir_imgs, exist_ok=True)
		
		
		self.STGs    = np.array( p.STGs ) * 1000
		self.GluN2Bs = np.array( p.GluN2Bs ) * 100
		
		self.num_rows		= len( p.GluN2Bs )
		self.num_columns	= len( p.STGs )
		
	def prepare_plot( self ):
		self.fig  = plt.figure(figsize=(5, 5))
		self.fig.subplots_adjust(wspace=0.4,  hspace=0.6)
		ax = self.fig.add_subplot( 1, 1, 1 )
		ax.set_title( self.title )
		ax.set_xlabel('STG (beads / voxel) x 10-3')	
		ax.set_ylabel('GluN2B (beads / voxel) x 10-2')
		return ax
		
		
	def run( self ):
		#
		data = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		#
		for i, stg in enumerate(self.STGs):
			for j, glun in enumerate(self.GluN2Bs):
				id = i + j * len(self.STGs)
				prefix = str(id).zfill(3)
				print('Target file: ', prefix)
				d           = utils.load(self.dir_edited_data, prefix, self.suffix)
				data[i,j]   = self._modify_data(d)
		
		print('data ', data)
		ax = self.prepare_plot()
		cs, cb = utils.plot_a_panel(ax, data, self.STGs, self.GluN2Bs, self.colormap, self.levels)
		
		
	def save( self ):
		
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.svg' ) )
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.png' ) , dpi=150)
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
		
		
class PlotPhaseDiagram( PhaseDiagramConcDependence ):
	def __init__( self ):
		
		super().__init__()
		
		GluN2Bs = ([-570, 570, 1080, 4320, 6480, 8640, 10800, 12960, 17280])
		GluN2Bs.reverse()
		volume = np.prod(p.space_np)
		GluN2Bs = [ n / volume for n in GluN2Bs ]
		self.GluN2Bs  = np.array( GluN2Bs ) * 100
		self.num_rows = len( GluN2Bs )
		self.title    = 'Phase diagram'
		self.basename = 'phase_diagram_conc_dependence'
		# -1: Unclear
		# 1: Homogeneous LLPS (STG)
		# 2: PIPS
		# 3: Partial engulfment
		# 4: Homogeneous LLPS (CaMKII)
		
		phase_diagram = [\
			[ 1, 1, 1, 1,   1, 1], # 17280
			[ 1, 1, 1, 1,   2, 2], # 12960
			[ 1, 2, 2, 2,   2, 2], # 10800
			[ 2, 2, 2, 2,   2, 2], # 8640
			[ 2, 3, 2, 2,   2, 2], # 6480
			[ 3, 3, 3, 3,   3, 3], # 4320
			[ 4, 4, 4, 3,   3, 3], # 1080
			[ 4, 4, 4, 4,   4, 4], # 570
			[ 4, 4, 4, 4,   4, 4]] # -570
		self.phase_diagram = np.array(phase_diagram).T
		
	def plot( self ):
		
		colormap = c.cmap_phase_diagram1
		levels   = np.array([0.0,1.5,2.5,3.5,4.5])
		# cmap_gray_cr_pk_gray # c.cmap_white_green_universal, plt.colormaps['jet']# 'winter', etc
		
		ax = self.prepare_plot()
		cs, cb = utils.plot_a_panel(ax, self.phase_diagram, self.STGs, self.GluN2Bs, colormap, levels, draw_border = True)
		
class PlotConnectivityPhaseDiagramConcDependence(PlotConnectivity, PhaseDiagramConcDependence):
	def __init__( self, species, type_analysis ):
		PlotConnectivity.__init__(self, species, type_analysis )
		PhaseDiagramConcDependence.__init__(self)
		
class PlotCondVolumePhaseDiagramConcDependence(PlotCondVolume, PhaseDiagramConcDependence):
	def __init__( self, species ):
		PlotCondVolume.__init__(self, species )
		PhaseDiagramConcDependence.__init__(self)
		
if __name__ == '__main__':
	
	'''
	pl = PlotPhaseDiagram()
	pl.plot()
	pl.save()
	'''
	
	#species, type_analysis = 'CaMKII', 'average'
	species, type_analysis = 'PSD95' , 'average'
	#species, type_analysis = 'PSD95' , 'ratio'
	'''
	pl = PlotConnectivityPhaseDiagramConcDependence(species, type_analysis)
	pl.run()
	pl.save()
	'''
	
	species = 'STG' # 'CaMKII', 'STG'
	pl = PlotCondVolumePhaseDiagramConcDependence(species)
	pl.run()
	pl.save()
