
import os, glob, pickle, pprint, copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# import mpl_toolkits.axes_grid1
# from scipy.interpolate import griddata, RegularGridInterpolator

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c

from skimage.measure import label
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit
import warnings


plt.rcParams.update( p.rc_param )


def exponential(x, tau, a): 
	return  a * np.exp(x/tau)


def soft_plus(x, tau, a, x0):
	return  a * np.log( 1+np.exp((x-x0)/tau) )


def x_n(x, a, center, n):
	return  a * (x - center)**n


def prepare_plot():
	fig  = plt.figure(figsize=(3.0, 3.5)) # (3.5, 3.5)
	fig.subplots_adjust(left = 0.2, bottom = 0.2, wspace=0.4,  hspace=0.6)
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



		
if __name__ == '__main__':
	
	# Shared init
	dir_target = 'small_colony2'
	dir_edited_data = os.path.join('data4', dir_target)
	suffix = 'matrix'
	
	dir_imgs = os.path.join('imgs4', dir_target, 'fitting')
	os.makedirs(dir_imgs, exist_ok=True)
	filename_img = 'FRAP_clustering_coefficient_fitting'
	
	# FRAP
	prefix = 'FRAP_merged'
	d_FRAP = utils.load(dir_edited_data, prefix, suffix)
	
	# Clustering coefficient.
	prefix = 'average_clustering_coefficient'
	d_cluster = utils.load(dir_edited_data, prefix, suffix)
	
	x = [d_FRAP[k] for k in d_FRAP.keys() if k not in ['04_000', '04_001', '04_006']]
	y = [d_cluster[k] for k in d_FRAP.keys() if k not in ['04_000', '04_001', '04_006']]
	x = np.array(x)
	y = np.array(y)
	
	
	fitting_params_exp , pcov = curve_fit(exponential, x, y, inits = np.array([0.5, 0.02]))
	fitting_params_soft, pcov = curve_fit(soft_plus  , x, y)
	
	y_pred =  exponential(x, *fitting_params_exp.tolist())
	r2_exp = r2_score(y, y_pred)
	
	y_pred  =  soft_plus(x, *fitting_params_soft.tolist())
	r2_soft = r2_score(y, y_pred)

	
	print('fitting parameters (exponential): ', fitting_params_exp)
	print('R2 statistics: '     , r2_exp)
	print()
	print('fitting parameters (soft plus)  : ', fitting_params_soft)
	print('R2 statistics: '     , r2_soft)
	print()

	
	title = 'Cluster coeeficient'
	xlim  = [-2.5, 2.5]
	fig, ax = prepare_plot()
	#ax.set_yscale('log')
	ax.set_xlim(xlim)
	ax.set_ylim([-0.03,0.3])
	ax.set_yticks([0,0.1,0.2,0.3])
	ax.set_xlabel('FRAP ()')
	ax.set_ylabel('Clustering coefficient')
	
	xx = np.linspace(xlim[0],xlim[1],40)
	#ax.plot( xx, exponential( xx,*fitting_params_exp.tolist() ), '-', color = 'b', label='R2 (exp): {:.3f}'.format(r2_exp) )
	ax.plot( xx, soft_plus( xx,*fitting_params_soft.tolist() ), '-', color = (0.5,0.5,0.5), label='R2 (soft): {:.3f}'.format(r2_soft))
	ax.plot(x, y,'o', \
		markersize = 4, \
		color = 'k', \
		markerfacecolor = 'k' )
	ax.legend(frameon=False)
	plot_and_save( fig, dir_imgs, filename_img )

