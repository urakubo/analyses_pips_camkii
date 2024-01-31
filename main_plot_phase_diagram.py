
import os, glob, pickle, pprint
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import mpl_toolkits.axes_grid1

from scipy.interpolate import griddata, RegularGridInterpolator
import utils
import colormap as c

plt.rcParams.update({
                    'pdf.fonttype' : 'truetype',
                    'svg.fonttype' : 'none',
                    'font.family' : 'sans-serif',
                    'font.sans-serif' : 'Arial',
                    'font.style' : 'normal'})


# Extrapolation works, but only for grid data.
def make_grid_using_RegularGridInterpolator(STG, GluN2B, oZ, mX, mY):
	m_points = np.array([mX.ravel(), mY.ravel()]).T
	f = RegularGridInterpolator((STG, GluN2B), oZ, bounds_error=False, fill_value=None)
	mZ = f(m_points).reshape(mX.shape)
	return mZ

# Interpolation only, but scattered data are applicable.
def make_grid_using_griddata(ox, oy, oz, mX, mY):
	mZ = griddata((ox, oy), oz, (mX, mY), method='cubic') #  ‘nearest’'cubic'
	return mZ


def plot_a_panel(ax, oZ, STG, GluN2B, colormap, levels, method = 'RegularGridInterpolator'):
	
	# Observed data arrangement
	oX, oY  = np.meshgrid(STG, GluN2B)
	ox      = np.ravel(oX)
	oy      = np.ravel(oY)
	oz      = np.ravel(oZ.T)
	oz_max  = np.max(oz)
	
	
	# Mesh grids for interpolation
	mx_max = np.max(STG)
	my_max = np.max(GluN2B)
	mx = np.linspace(mx_max*0.02, mx_max*1.1, 55)
	my =  np.linspace(my_max*0.02, my_max*1.1, 55)
	mX, mY = np.meshgrid(mx,my)
	
	if   method == 'RegularGridInterpolator':
		mZ = make_grid_using_RegularGridInterpolator(STG, GluN2B, oZ, mX, mY)
	elif method == 'griddata':
		mZ = make_grid_using_griddata(ox, oy, oz, mX, mY)
	else:
		raise ValueError('Invalid method')
	
	mZ[mZ > np.max(levels)] = np.max(levels)
	mZ[mZ < 0]      = 0
	
	# Plot
	cs = ax.contourf(mX, mY, mZ, levels=levels, alpha=0.3, \
				cmap= colormap ) # vmin=0, vmax=np.max(levels), 

	ax.scatter(ox, oy, c=oz, cmap=colormap, marker='o', edgecolors='k', s=16, vmin=0, vmax=np.max(levels))
	
	ax.set_xlim(0, mx_max*1.1)
	ax.set_ylim(0, my_max*1.1)
	
	ax.set_xlabel('STG (number)')
	ax.set_ylabel('GluN2B (number)')
	ax.set_box_aspect(1)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	
	divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
	cax = divider.append_axes('right', '5%', pad='3%')
	cb = plt.colorbar(cs, cax=cax)
	cb.ax.set_yticklabels(["{:.2f}".format(i) for i in cb.get_ticks()])
	
	return cs, cb
	
	
def plot_panels_partition_index(STG, GluN2B, partition_i_CaMKII, partition_i_STG, threshold ):
	
	num_levels = 10   #
	max_pi = np.max(partition_i_CaMKII)
	
	fig  = plt.figure(figsize=(8, 5))
	fig.subplots_adjust(wspace=0.4,  hspace=0.6)


	colormap   = plt.colormaps['jet'] # c.cmap_white_green_universal, plt.colormaps['jet']# 'winter', etc
	levels  = np.linspace(0, max_pi, 10)
	ax = fig.add_subplot( 2, 2, 1 )
	cs, cb = plot_a_panel(ax, partition_i_CaMKII, STG, GluN2B, colormap, levels)
	cb.set_label('(Ratio)')
	ax.set_title('Partition index of CaMKII')


	partition_i_CaMKII = (partition_i_CaMKII > threshold)
	colormap   = c.cmap_white_green_binary_universal # 'winter', etc
	levels  = np.array([0,0.5,1.0])
	ax = fig.add_subplot( 2, 2, 2 )
	cs, cb = plot_a_panel(ax, partition_i_CaMKII, STG, GluN2B, colormap, levels)
	cb.set_label('(Existence)')
	ax.set_title('Inhomogenous CaMKII')


	colormap = plt.colormaps['jet'] # c.cmap_white_red_universal, plt.colormaps['jet', 'winter]
	levels   = np.linspace(0, max_pi, 10)
	ax = fig.add_subplot( 2, 2, 3 )
	cs, cb = plot_a_panel(ax, partition_i_STG, STG, GluN2B, colormap, levels)
	cb.set_label('(Ratio)')
	ax.set_title('Partition index of STG')


	partition_i_STG = (partition_i_STG > threshold)
	colormap   = c.cmap_white_red_binary_universal # 'winter', etc
	levels  = np.array([0,0.5,1.0])
	ax = fig.add_subplot( 2, 2, 4 )
	cs, cb = plot_a_panel(ax, partition_i_STG, STG, GluN2B, colormap, levels)
	cb.set_label('(Existence)')
	ax.set_title('Inhomogenous STG')

	return fig


	
def plot_panels_watershed(STG, GluN2B, volume_watershed_CaMKII, volume_watershed_STG):
	
	num_levels = 10   #
	max_conc = np.max(volume_watershed_CaMKII)
	
	fig  = plt.figure(figsize=(8, 5))
	fig.subplots_adjust(wspace=0.4,  hspace=0.6)


	colormap   = plt.colormaps['jet'] # c.cmap_white_green_universal, plt.colormaps['jet']# 'winter', etc
	levels  = np.linspace(0, max_conc, 10)
	ax = fig.add_subplot( 2, 2, 1 )
	cs, cb = plot_a_panel(ax, volume_watershed_CaMKII, STG, GluN2B, colormap, levels)
	cb.set_label('(% Total volume)')
	ax.set_title('Basin volume of CaMKII')


	colormap   = c.cmap_white_green_binary_universal # 'winter', etc
	volume_watershed_CaMKII = (volume_watershed_CaMKII > 0.1)
	levels  = np.array([0,0.5,1.0])
	ax = fig.add_subplot( 2, 2, 2 )
	cs, cb = plot_a_panel(ax, volume_watershed_CaMKII, STG, GluN2B, colormap, levels)
	cb.set_label('(Existence)')
	ax.set_title('CaMKII basin')


	colormap   = plt.colormaps['jet'] # c.cmap_white_red_universal, plt.colormaps['jet']# 'winter', etc
	levels  = np.linspace(0, max_conc, 10)
	ax = fig.add_subplot( 2, 2, 3 )
	cs, cb = plot_a_panel(ax, volume_watershed_STG, STG, GluN2B, colormap, levels)
	cb.set_label('(% Total volume)')
	ax.set_title('Basin volume of STG')


	colormap   = c.cmap_white_red_binary_universal # 'winter', etc
	volume_watershed_STG = (volume_watershed_STG > 0.1)
	levels  = np.array([0,0.5,1.0])
	ax = fig.add_subplot( 2, 2, 4 )
	cs, cb = plot_a_panel(ax, volume_watershed_STG, STG, GluN2B, colormap, levels)
	cb.set_label('(Existence)')
	ax.set_title('STG basin')

	return fig

def load(dir_data, sigma, STG, GluN2B, target = 'volume_watershed'):
	
	CaMKII = np.zeros([len(STG), len(GluN2B)], dtype=float)
	STG    = np.zeros([len(STG), len(GluN2B)], dtype=float)
	tot_volume = utils.space[0] * utils.space[1] * utils.space[2]
	
	for i, stg in enumerate(STG):
		for j, glun in enumerate(GluN2B):
			id = i + j * len(STG)
			prefix = str(id).zfill(3)
			suffix = 'sigma_{}'.format(sigma)
			d = utils.load(dir_data, prefix, suffix)
			if target == 'volume_watershed':
				CaMKII[i,j] = np.sum( d['labels_watershed_in_grid_mesh']['CaMKII'] ) / tot_volume * 100
				STG[i,j]    = np.sum( d['labels_watershed_in_grid_mesh']['STG']    ) / tot_volume * 100
			if target == 'partition_index':
				cc = d['conc_condensate']
				target_conc, target_condensate, ref_condensate  = 'CaMKII', 'CaMKII', 'STG'
				CaMKII[i,j] = cc[target_condensate][target_conc] / cc[ref_condensate][target_conc]
				target_conc, target_condensate, ref_condensate  = 'STG',  'STG', 'CaMKII'
				STG[i,j]    = cc[target_condensate][target_conc] / cc[ref_condensate][target_conc]
	return CaMKII, STG



if __name__ == '__main__':
	
	# Shared parameters
	dir_data = 'data'
	dir_imgs = 'imgs/phase_diagram'
	os.makedirs(dir_imgs, exist_ok=True)
	sigma  = 4 # 4,3,2
	STG    = [500, 1000, 2000, 3000, 4000]
	GluN2B = [500, 2000, 4000, 6000, 8000, 12000]

	# Target 1: Watershed-based phase diagram
	'''
	target = 'volume_watershed'
	vCaMKII, vSTG = load(dir_data, sigma, STG, GluN2B, target = target)
	fig = plot_panels_watershed(STG, GluN2B, vCaMKII, vSTG )
	fig.savefig( os.path.join(dir_imgs, '{}_sigma_{}.svg'.format( target, sigma ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}_sigma_{}.png'.format( target, sigma ) ) , dpi=150)
	'''
	
	# Target 2: Partition-index based phase diagram
	#'''
	target = 'partition_index'
	vCaMKII, vSTG = load(dir_data, sigma, STG, GluN2B, target = target)
	
	threshold = 2
	fig = plot_panels_partition_index(STG, GluN2B, vCaMKII, vSTG, threshold )
	
	fig.savefig( os.path.join(dir_imgs, '{}_sigma_{}_th_{}.svg'.format( target, sigma, threshold ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}_sigma_{}_th_{}.png'.format( target, sigma, threshold ) ) , dpi=150)
	#plt.show()
	plt.clf()
	plt.close(fig=fig)
	#'''

