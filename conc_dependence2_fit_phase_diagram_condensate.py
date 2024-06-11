
import os, glob, pickle, pprint, copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import lib.utils as utils
import lib.utils_fitting as utils_fitting

import lib.parameters as p
import lib.colormap as c


from skimage.measure import label
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import r2_score

from scipy.optimize import curve_fit
import warnings

from conc_dependence2_calc_connectivity_plot_phase_diagram import \
	HandleConnectivityPhaseDiagramConcDependence, \
	HandleCondVolumePhaseDiagramConcDependence, \
	PlotPhaseDiagram


plt.rcParams.update( p.rc_param )


	
def prepare_plot():
	# fig  = plt.figure(figsize=(5, 4))
	fig  = plt.figure(figsize=(4, 4))
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
	
	
	
def plot_fit_PIPS_partial_engulfment( fig, ax, title, ratios_connection ):
	
	pl = PlotPhaseDiagram()
	two_condensates    = pl.phase_diagram_two_condensates
	partial_engulfment = pl.phase_diagram_partial_engulfment
	two_condensates    = np.fliplr(two_condensates)
	partial_engulfment = np.fliplr(partial_engulfment)
	
	
	two_condensates    = two_condensates[:-1,:]
	partial_engulfment = partial_engulfment[:-1,:]
	ratios_connection  = ratios_connection[:-1,:]
	
	'''
	print('two_condensates')
	print(two_condensates)
	print('ratios_connection')
	print(ratios_connection)
	'''
	
	pips = two_condensates * (partial_engulfment == 0)
	
	
	x = ratios_connection.reshape(-1, 1)
	y = pips.astype('int').reshape(-1, 1)
	
	model = utils_fitting.logistic_regression(x,y)
	
	species_vol = 'PIPS'
	ax.set_xlabel(title)
	
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
	ax.set_yticklabels(['Others','PIPS'])
	
	
def hill(x, a, b, c): # Hill sigmoidal equation from zunzun.com
	return  a * np.power(x, b) / (np.power(c, b) + np.power(x, b))
	
	
def plot_fit_hill_volumes( fig, ax, title, numbers_connection, volumes, species_vol, xmax ):
	
	volumes = volumes[:-1,:]
	numbers_connection  = numbers_connection[:-1,:]
	x = np.ravel( numbers_connection )
	y = np.ravel( volumes ) / p.volume
	
	initialParameters = np.array([0.01, 3, xmax/10])
	bounds = ((0, 1, 1), (0.08, 4, 6))
	bounds = ((0, 0.1, 0.1), (0.08, 6, xmax))
	#warnings.filterwarnings("ignore")
	fitting_params, pcov = curve_fit(hill, x, y, initialParameters, bounds = bounds )
	#fittedParameters = initialParameters 
	
	xx = np.linspace(0.0,xmax,40)
	yy = hill(xx, *fitting_params)
	
	y_pred =  hill(x, *fitting_params)
	r2 = r2_score(y, y_pred)
	
	print('fitting parameters')
	print(fitting_params)
	print()
	print('R2 statistics')
	print(r2)
	
	ax.set_xlim([0, xmax])
	ax.plot(xx, yy * 100, '-', color = (0.5,0.5,0.5))
	ax.plot(x, y * 100 ,'o', \
		markersize = 4, \
		color = c.cmap_universal_ratio[species_vol], \
		markerfacecolor = c.cmap_universal_ratio[species_vol] )
	#ax.set_yticks([0,1,2,3,4])
	ax.set_yticks([0, 0.5, 1.0, 1.5])
	ax.set_ylabel('(% of total volume)')
	ax.set_title('R2: {:.4f}, param: a = {:.4f}, b = {:.4f}, c = {:.4f}'.format(r2, *fitting_params))
	
	
	
if __name__ == '__main__':
	
	#species, type_analysis, species_vol, xmax = 'CaMKII', 'average','CaMKII', 12
	#species, type_analysis, species_vol, xmax = 'PSD95' , 'average', 'STG', 3
	species, type_analysis = 'PSD95' , 'ratio'
	pl = HandleConnectivityPhaseDiagramConcDependence(species, type_analysis)
	pl.load_data()
	numbers_ratios_connection = pl.data	
	title    = pl.title
	dir_imgs = pl.dir_imgs
	basename = 'fit_' + species + '_' + type_analysis
	
	
	fig, ax = prepare_plot()
	
	if type_analysis == 'ratio':
		plot_fit_PIPS_partial_engulfment(fig, ax, title, numbers_ratios_connection )
	else:
		pl = HandleCondVolumePhaseDiagramConcDependence(species_vol)
		pl.load_data()
		volumes = pl.data
		plot_fit_hill_volumes(fig, ax, title, numbers_ratios_connection, volumes, species_vol, xmax )
	
	save_plots( fig, dir_imgs, basename )
