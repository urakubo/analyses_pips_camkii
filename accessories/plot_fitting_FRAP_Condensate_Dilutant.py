
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


def x_n(x, a, center, n):
	return  a * (x - center)**n


def prepare_plot():
	fig  = plt.figure(figsize=(3.8, 3.8)) # (3.5, 3.5)
	fig.subplots_adjust(left = 0.2, bottom = 0.2, wspace=0.4,  hspace=0.6)
	ax = fig.add_subplot( 1, 1, 1 )
	
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	return fig, ax

def plot_a_graph(ax, x, y, title):
	xlim  = [-4, 0.7]
	#ax.set_yscale('log')
	ax.set_xlim(xlim)
	ax.set_ylim(xlim)
	#ax.set_yticks([0,0.1,0.2])
	ax.set_xlabel('Recovery time on condensate \n (log[10^9 MC steps])')
	ax.set_ylabel('Recovery time on dilutant \n  (log[10^9 MC steps])')
	ax.set_title( title )
	
	ax.plot( xlim, xlim, '--', color = (0.5,0.5,0.5))
	ax.plot(x, y,'o', \
		markersize = 4, \
		color = 'k', \
		markerfacecolor = 'k' )
	ax.set_aspect('equal', 'box')
	ax.legend(frameon=False)

def save_plots( fig, dir_imgs, img_basename ):
	
	fig.savefig( os.path.join(dir_imgs, img_basename + '.svg' ) )
	fig.savefig( os.path.join(dir_imgs, img_basename + '.png' ) , dpi=150)
	plt.show()
	plt.clf()
	plt.close(fig=fig)



if __name__ == '__main__':
	
	# Shared init
	filename_img = 'FRAP_Condansate_dilutent'
	
	# Load data: FRAP condensate
	t = SpecDatasets()
	t.CG_valency_length_only_local_move()
	prefix = 'FRAP_taus_CaMKII_merged'
	suffix = 'matrix'
	d_FRAP = utils.load(t.dir_edited_data, prefix, suffix)
	
	
	# Load data: FRAP dilutant
	tt = SpecDatasets()
	tt.C_valency_length_FRAP_Control_fine_sampling(frap = True)
	prefix = 'FRAP_taus_CaMKII'
	suffix = 'matrix'
	d_cluster = utils.load(tt.dir_edited_data, prefix, suffix)
	dir_imgs = os.path.join(tt.dir_imgs_root, 'fitting')
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	x = d_FRAP.flatten()
	x = np.log10(x)
	y = d_cluster.flatten()
	y = np.log10(y)
	
	# Plot a figure and save it.
	fig, ax = prepare_plot()
	plot_a_graph(ax, x, y, filename_img)
	save_plots( fig, dir_imgs, filename_img)

