
import os, glob, pickle, pprint, copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from skimage.measure import label
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit
import warnings


import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
from specification_datasets import SpecDatasets


plt.rcParams.update( p.rc_param )


def exponential(x, tau, a): 
	return  a * np.exp(x/tau)


def soft_plus(x, tau, a, x0):
	return  a * np.log( 1+np.exp((x-x0)/tau) )


def linear(x, a, b):
	return  a * x + b


def prepare_plot():
	fig  = plt.figure(figsize=(3.0, 3.0)) # (3.5, 3.5)
	fig.subplots_adjust(left = 0.2, bottom = 0.2, wspace=0.4,  hspace=0.6)
	ax = fig.add_subplot( 1, 1, 1 )
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	
	return fig, ax

def plot_a_graph(ax, x, y, title, fitting_function, fitting_params, label):
	xlim  = [-3.2, 1]
	ylim  = [-3.2, 0]
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
	#ax.set_yticks([0,0.1,0.2])
	ax.set_xlabel('FRAP (log[10^9 MC steps])')
	ax.set_ylabel('Clustering coefficient')
	ax.set_title( title )
	
	xx = np.linspace(xlim[0],xlim[1],40)
	ax.plot( xx, fitting_function( xx,*fitting_params.tolist() ), '-', color = (0.5,0.5,0.5), label=label)
	ax.plot(x, y,'o', \
		markersize = 4, \
		color = 'k', \
		markerfacecolor = 'k' )
	ax.legend(frameon=False)

def save_plots( fig, dir_imgs, img_basename ):
	
	fig.savefig( os.path.join(dir_imgs, img_basename + '.svg' ) )
	fig.savefig( os.path.join(dir_imgs, img_basename + '.png' ) , dpi=150)
	plt.show()
	plt.clf()
	plt.close(fig=fig)



		
if __name__ == '__main__':
	
	t = SpecDatasets()
	t.CG_valency_length_only_local_move()
	
	dir_imgs = os.path.join(t.dir_imgs_root, 'fitting')
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	# Shared init
	filename_img = 'FRAP_clustering_coefficient_fitting'
	
	# Load data: FRAP
	prefix = 'FRAP_taus_CaMKII_merged'
	suffix = 'matrix'
	d_FRAP = utils.load(t.dir_edited_data, prefix, suffix)
	
	
	# Load data: Clustering coefficient.
	prefix = 'average_clustering_coefficient'
	d_cluster = utils.load(t.dir_edited_data, prefix, suffix)
	
	
	
	x = d_FRAP.flatten()
	x = np.log10(x)
	y = d_cluster.flatten()
	y = np.log10(y)
	
	
	# Fitting
	
	title = 'linear'
	fitting_function = linear
	inits = np.array([1, 0])
	
	fitting_params , pcov = curve_fit(fitting_function, x, y, p0= inits)
	
	y_pred = fitting_function(x, *fitting_params.tolist())
	r2_exp = r2_score(y, y_pred)
	
	label = 'R2: {:.3f}'.format(r2_exp)
	print('fitting parameters: ', fitting_params)
	print(label)
	print()

	# Plot a figure and save it.
	fig, ax = prepare_plot()
	plot_a_graph(ax, x, y, title, fitting_function, fitting_params, label)
	save_plots( fig, dir_imgs, filename_img + '_' + title)

