
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


def prepare_plot():
	fig  = plt.figure(figsize=(5, 4))
	fig.subplots_adjust(left = 0.20, bottom=0.2)
	ax = fig.add_subplot( 1, 1, 1 )
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	return fig, ax

def plot_and_save( fig, dir_imgs, basename ):
	
	fig.savefig( os.path.join(dir_imgs, basename + '.svg' ) )
	fig.savefig( os.path.join(dir_imgs, basename + '.png' ) , dpi=150)
	plt.show()
	plt.clf()
	plt.close(fig=fig)


class ConcDependence():
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
		return data
		
		



class GetConnectivity():
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



class GetCondVolume():
	def __init__( self, species ):
		
		if species == 'CaMKII':
			self.title    = 'Volume of {} condensate'.format(species)
			self.basename = 'volume_largest_cond_CaMKII'
			self.colormap =  c.cmap_white_green_universal
			self.levels   = np.linspace(0,7e4,8)
		elif species == 'STG':
			self.title    = 'Volume of {} condensate'.format(species)
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
		
		
		
class GetConnectivityConcDependence(GetConnectivity, ConcDependence):
	def __init__( self, species, type_analysis ):
		GetConnectivity.__init__(self, species, type_analysis )
		ConcDependence.__init__(self)
		
class GetCondVolumeConcDependence(GetCondVolume, ConcDependence):
	def __init__( self, species ):
		GetCondVolume.__init__(self, species )
		ConcDependence.__init__(self)
		
if __name__ == '__main__':
	
	
	#species, type_analysis = 'PSD95' , 'ratio'
	
	species, type_analysis = 'PSD95' , 'average'
	species_vol = 'STG' # 'CaMKII', 'STG'
	
	#species, type_analysis = 'CaMKII', 'average'
	#species_vol = 'CaMKII'
	
	dir_target= 'conc_dependence'
	
	# Shared parameters
	prefix = 'conn_volume'
	suffix = species + '_' + species_vol
	dir_edited_data  = os.path.join('data3',dir_target, prefix)
	os.makedirs(dir_edited_data, exist_ok=True)
	dir_imgs = os.path.join('imgs3', dir_target, prefix)
	os.makedirs(dir_imgs, exist_ok=True)
	basename = suffix
	
	
	conn = GetConnectivityConcDependence(species, type_analysis)
	vol  = GetCondVolumeConcDependence(species_vol)
	
	'''
	numbers_connection = conn.run()
	volumes = vol.run()
	d = {'numbers_connection':numbers_connection, 'volumes':volumes}
	utils.save(dir_edited_data, prefix, suffix, d)
	'''
	
	
	d   = utils.load(dir_edited_data, prefix, suffix)
	numbers_connection = d['numbers_connection']
	volumes = d['volumes']
	
	numbers_connection = np.ravel(numbers_connection)
	volumes = np.ravel(volumes) / p.volume
	
	fig, ax = prepare_plot()
	ax.plot(numbers_connection, volumes * 100,'o', \
		markersize = 4, \
		color = c.cmap_universal_ratio[species_vol], \
		markerfacecolor = c.cmap_universal_ratio[species_vol] )
	#ax.set_title( 'Connectivity-Volume' )
	ax.set_xlabel(conn.title)	
	ax.set_ylabel('(% of total volume)')
	
	plot_and_save( fig, dir_imgs, basename )
	
