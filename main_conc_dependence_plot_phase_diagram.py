
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
	mx_min = np.min(STG)
	my_min = np.min(GluN2B)
	
	mx = np.linspace(0.0, mx_max*1.1, 55)
	my = np.linspace(0.0, my_max*1.1, 55)


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
	
	ax.set_xlim( np.min(mx), np.max(mx) )
	ax.set_ylim( np.min(my), np.max(my) )
	
	
	ax.set_xlabel('STG (beads / volume) x 10-3')
	ax.set_ylabel('GluN2B (beads / volume) x 10-2')
	ax.set_box_aspect(1)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	
	divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
	cax = divider.append_axes('right', '5%', pad='3%')
	cb = plt.colorbar(cs, cax=cax)
	cb.ax.set_yticklabels(["{:.2f}".format(i) for i in cb.get_ticks()])
	
	return cs, cb
	
	
if __name__ == '__main__':
	
	
	# Output files
	dir_imgs = os.path.join('imgs', 'conc_dependence','phase_diagram')
	filename_output = 'phase_diagram_conc_dependence'
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	# GluN2B = [1, 500, 1000, 2000, 4000, 6000, 8000, 12000, 16000, 20000]
	STG    = [1, 500, 1000, 1500, 2000, 2500, 3000] 
	GluN2B = [20000, 16000, 12000, 8000, 6000, 4000, 2000, 1000, 500, 1]
	
	volume = np.prod(utils.space_np)
	STG    = [ s / volume for s in STG    ]
	GluN2B = [ n / volume for n in GluN2B ]

	STG    = np.array( STG ) * 1000
	GluN2B = np.array( GluN2B ) * 100

	print('STG')
	print(STG)
	print('GluN2B')
	print(GluN2B)
	
	# -1: Unclear
	# 0: Homogeneous LLPS (STG)
	# 1: Homogeneous LLPS (CaMKII)
	# 2: Partial engulfment
	# 3: PIPS

	phase_diagram = [\
	[ 1, 1, 1, 1,   1, 1, 1], # 20000
	[ 1, 1, 1, 1,   2, 2, 2], # 16000
	[ 1, 1, 2, 2,   2, 2, 2], # 12000
	[ 1, 2, 2, 2,   2, 2, 3], # 8000
	[ 1, 2, 2, 2,   3, 3, 3], # 6000
	[ 1, 2, 2, 2,   3, 3, 3], # 4000
	[ 1, 2, 2, 2,   3, 3, 3], # 2000
	[ 1, 2, 2, 3,   4, 4, 4], # 1000
	[ 1, 3, 4, 4,   4, 4, 4], # 500
	[ 1, 4, 4, 4,   4, 4, 4]] # 1
	
	
	phase_diagram = np.array(phase_diagram).T
	
	
	num_levels = 4   #
	max_conc = np.max(phase_diagram)
	
	fig  = plt.figure(figsize=(5, 5))
	fig.subplots_adjust(wspace=0.4,  hspace=0.6)
	
	colormap   = plt.colormaps['jet'] # c.cmap_white_green_universal, plt.colormaps['jet']# 'winter', etc
	levels  = np.linspace(0, max_conc, num_levels)
	levels  = np.array([-0.5,0.5,1.5,2.5,3.5,4.5])
	ax = fig.add_subplot( 1, 1, 1 )
	cs, cb = plot_a_panel(ax, phase_diagram, STG, GluN2B, colormap, levels)
	# cb.set_label('(% Total volume)')
	ax.set_title('Phase diagram')
	
	
	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename_output ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename_output ) ) , dpi=150)
	
	plt.show()
	plt.clf()
	plt.close(fig=fig)
	
	
