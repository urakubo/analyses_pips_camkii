
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
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import r2_score

from scipy.optimize import curve_fit
import warnings


plt.rcParams.update( p.rc_param )


def efron_rsquare(y, y_pred):
	n = float(len(y))
	t1 = np.sum((y - y_pred) * (y - y_pred))
	t2 = np.sum((y - np.average(y)) * (y - np.average(y)))
	return 1.0 - (t1 / t2)

def get_count_most_freq_outcome(y):
	num_0 = 0
	num_1 = 0
	for p in y:
		if p == 1.0:
			num_1 += 1
		else:
			num_0 += 1
	return float(max(num_0, num_1))

def count_adjusted_rsquare(y, y_pred, t=0.5):
    correct = get_num_correct(y, y_pred, t)
    total = float(len(y))
    n = get_count_most_freq_outcome(y)
    return (correct - n) / (total - n)


def get_num_correct(y, y_pred, t=0.5):
	y_correct = np.array([0.0 if p < t else 1.0 for p in y_pred])
	return sum([1.0 for p, p_pred in zip(y, y_correct) if p == p_pred])

def count_rsquare(y, y_pred, t=0.5):
	n = float(len(y))
	num_correct = get_num_correct(y, y_pred, t)
	return num_correct / n


def hill(x, a, b, c): # Hill sigmoidal equation from zunzun.com
	return  a * np.power(x, b) / (np.power(c, b) + np.power(x, b))


def prepare_plot():
	fig  = plt.figure(figsize=(5, 4))
	fig.subplots_adjust(left = 0.20, bottom=0.2)
	ax = fig.add_subplot( 1, 1, 1 )
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	return fig, ax

def save_plots( fig, dir_imgs, basename ):
	
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
		self.dir_edited_data = os.path.join('data4', dir_target)
		self.dir_imgs        = os.path.join('imgs4', dir_target, 'phase_diagram')
		
		os.makedirs(self.dir_imgs, exist_ok=True)
		
		
		self.STGs    = np.array( p.STGs ) * 1000
		self.GluN2Bs = np.array( p.GluN2Bs ) * 1000
		
		self.num_rows		= len( p.GluN2Bs )
		self.num_columns	= len( p.STGs )
		
		
	def run( self ):
		#
		data = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		for i, stg in enumerate(self.STGs):
			for j, glun in enumerate(self.GluN2Bs):
				id = i + j * len(self.STGs)
				prefix = str(id).zfill(3)
				print('Target file: ', prefix)
				
				d           = utils.load(self.dir_edited_data, prefix, self.suffix)
				data[i,j]   = self._modify_data(d, id)
		
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
	
	def _modify_data(self, d, id):
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
		elif species == 'PIPS':
			self.title    = 'Volume of {} condensate'.format(species)
			self.basename = 'PIPS'
			self.colormap = plt.colormaps['Greys']
			self.levels   = np.linspace(0,1.0,11)
		else:
			raise ValueError("Not implemented, species: ", species)
		
		self.species  = species
		self.basename = 'max_vol_{}'.format( self.species )
		self.suffix   = 'sigma_2'
		
	def _modify_data(self, d, id):
		#
		if self.species == 'PIPS':
			pips = [31,32,33,34,35,24,25,26,27,28,29,18,20,21,22,23,12]
			partial_engulf = [19,13,14,15,16,17,9,10,11]
			if id in pips:
				return 1.0
			elif id in partial_engulf:
				return 0.0
			else:
				return None
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
	#species_vol = 'PIPS'
	
	species, type_analysis = 'PSD95' , 'average'
	species_vol = 'STG' # 'CaMKII', 'STG'
	
	#species, type_analysis = 'CaMKII', 'average'
	#species_vol = 'CaMKII'
	
	dir_target= 'conc_dependence_merged'
	
	# Shared parameters
	prefix = 'conn_volume'
	suffix = species + '_' + species_vol
	dir_edited_data  = os.path.join('data4',dir_target, prefix)
	os.makedirs(dir_edited_data, exist_ok=True)
	dir_imgs = os.path.join('imgs4', dir_target, prefix)
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
	
	fig, ax = prepare_plot()
	ax.set_xlabel(conn.title)
	
	
	
	
	if species_vol == 'PIPS':
		
		nan = np.inf
		volumes = np.array([[nan, nan,  1.,  1.,  1., nan, nan, nan],\
		 	[nan, nan,  0.,  0.,  1.,  1., nan, nan],\
		 	[nan, nan,  0.,  1.,  1.,  1., nan, nan],\
		 	[nan,  0.,  0.,  1.,  1.,  1., nan, nan],\
		 	[nan,  0.,  0.,  1.,  1.,  1., 1., nan],\
		 	[nan,  0.,  0.,  1.,  1.,  1., 1., nan]])

		numbers_connection = np.ravel(numbers_connection)
		volumes = np.ravel(volumes)
		ids = np.nonzero(volumes != np.inf)
		x = numbers_connection[ids].reshape(-1, 1)
		y = volumes[ids].astype('int').reshape(-1, 1)
		
		model = LogisticRegression(penalty='none', solver = 'lbfgs')
		model.fit(x, y)
		# Model parameters
		print('w0: {:.4f}'.format( model.intercept_[0]) )
		print('w1: {:.4f}'.format( model.coef_[0][0]) )
		print('1/(1+exp[-(w0 + w1.x)])')
		
		y_pred = model.predict_proba(x)[:,1]
		y = np.ravel(y)
		y_pred = np.ravel(y_pred)
		w = np.array(model.coef_).transpose()
		
		print("Efron's  R2: {:.4f}".format(  efron_rsquare(y, y_pred) ) )
		print("Count R2   : {:.4f}".format(  count_rsquare(y, y_pred) ) )
		print("Adjust count R2: {:.4f}".format(  count_adjusted_rsquare(y, y_pred) ) )
		
		# https://datascience.oneoffcoder.com/psuedo-r-squared-logistic-regression.html
		# https://bookdown.org/egarpor/SSS2-UC3M/logreg-deviance.html
		xmax = np.max(x)*1.1
		xx = np.linspace(0.0,xmax,40).reshape(-1, 1)
		ax.plot([0, xmax], [0,0], ':', color = (0.5,0.5,0.5))
		ax.plot([0, xmax], [1,1], ':', color = (0.5,0.5,0.5))
		ax.plot(x, y,'o', \
			markersize = 4, \
			color = c.cmap_universal_ratio[species_vol], \
			markerfacecolor = c.cmap_universal_ratio[species_vol] )
		ax.plot(xx, model.predict_proba(xx)[:,1], '-', color = (0.5,0.5,0.5))
		ax.set_ylim([-0.1,1.1])
		ax.set_xlim([   0,xmax])
		ax.set_yticks([0,1])
		ax.set_yticklabels(['Partial \n engulfment','PIPS'])
		# ax.set_ylabel('(PIPS / No PIPS)')	
	else:
		
		x = np.ravel( numbers_connection )
		y = np.ravel( volumes ) / p.volume
		
		ids_in = np.nonzero(y < 0.033)
		xi = x[ids_in]
		yi = y[ids_in]
		ids_out = np.nonzero(y >= 0.033)
		xo = x[ids_out]
		yo = y[ids_out]
		
		xmax = np.max(conn.levels)
		print('xmax ', xmax)
		initialParameters = np.array([0.01, 3, xmax/10])
		bounds = ((0, 1, 1), (0.08, 4, 6))
		bounds = ((0, 0.1, 0.1), (0.08, 6, xmax))
		#warnings.filterwarnings("ignore")
		fittedParameters, pcov = curve_fit(hill, xi, yi, initialParameters, bounds = bounds )
		#fittedParameters = initialParameters 
		
		xx = np.linspace(0.0,xmax,40)
		yy = hill(xx, *fittedParameters)
		
		yi_pred =  hill(xi, *fittedParameters)
		r2 = r2_score(yi, yi_pred)
		
		print('fitted parameters')
		print(fittedParameters)
		print()
		print('R2 statistics')
		print(r2)
		
		ax.set_xlim([0, xmax])
		ax.plot(xx, yy * 100, '-', color = (0.5,0.5,0.5))
		ax.plot(xi, yi * 100 ,'o', \
			markersize = 4, \
			color = c.cmap_universal_ratio[species_vol], \
			markerfacecolor = c.cmap_universal_ratio[species_vol] )
		ax.plot(xo, yo * 100 ,'o', \
			markersize = 4, \
			color = c.cmap_universal_ratio[species_vol], \
			markerfacecolor = c.cmap_universal_ratio_light[species_vol] )
		ax.set_ylabel('(% of total volume)')
		ax.set_title('R2: {:.4f}, param: a = {:.4f}, b = {:.4f}, c = {:.4f}'.format(r2, *fittedParameters))
	save_plots( fig, dir_imgs, basename )
	
