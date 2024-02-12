
import os, glob, pickle, pprint, copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import mpl_toolkits.axes_grid1

from scipy.interpolate import griddata, RegularGridInterpolator
import utils
import colormap as c
import parameters as p

plt.rcParams.update( p.rc_param )


# Extrapolation works, but only for grid data.
def make_grid_using_RegularGridInterpolator(STG, GluN2B, oZ, mX, mY):
	m_points = np.array([mX.ravel(), mY.ravel()]).T
	f = RegularGridInterpolator((STG, GluN2B), oZ, bounds_error=False, fill_value=None)
	mZ = f(m_points).reshape(mX.shape)
	return mZ
	
	
def plot_a_panel(ax, oZ, STG, GluN2B, colormap, levels, method = 'RegularGridInterpolator'):
	
	# Observed data arrangement
	oX, oY  = np.meshgrid(STG, GluN2B)
	ox      = np.ravel(oX)
	oy      = np.ravel(oY)
	oz  = np.ravel(oZ.T)
	
	# Mesh grids for interpolation
	mx_max = np.max(STG)
	my_max = np.max(GluN2B)
	
	mx = np.linspace(0.0, mx_max*1.1, 55*4)
	my = np.linspace(0.0, my_max*1.1, 55*4)
	
	mX, mY = np.meshgrid(mx,my)
	
	#print('oz_graph.shape ', oz_graph.shape )
	
	oZ_panel = copy.deepcopy( oZ )
	oZ_panel[oZ_panel < 0] = 1
	mZ = make_grid_using_RegularGridInterpolator(STG, GluN2B, oZ_panel, mX, mY)
	
	# Plot
	colormap.set_bad(color='magenta')
	cs = ax.contourf(mX, mY, mZ, levels=levels, alpha=0.5, \
				cmap= colormap ) # vmin=0, vmax=np.max(levels), 
	ax.contour(cs, colors='k')
	ax.scatter(ox, oy, c=oz, cmap=colormap, marker='o', edgecolors='k', s=16, vmin=0, vmax=np.max(levels))
	
	
	# Overlay exception (not clean).
	oz_except = (oZ < 1)
	mZ_except = make_grid_using_RegularGridInterpolator(STG, GluN2B, oz_except, mX, mY)
	mZ_except[mZ_except > 0.5] = 1.0
	mZ_except[mZ_except <= 0.5] = np.nan
	print('np.unique(mZ_except) ' , np.unique(mZ_except) )
	cs = ax.contourf(mX, mY, mZ_except, vmin=0.3, vmax=0.4, cmap='binary' ) 
	#ax.set_facecolor("black")
	
	
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

	'''
	print('STG')
	print(STG)
	print('GluN2B')
	print(GluN2B)
	'''
	
	# -1: Unclear
	# 1: Homogeneous LLPS (STG)
	# 2: PIPS
	# 3: Partial engulfment
	# 4: Homogeneous LLPS (CaMKII)

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
	[-1, 4, 4, 4,   4, 4, 4]] # 1
	
	phase_diagram = np.array(phase_diagram).T
	
	
	num_levels = 40   #
	max_conc = np.max(phase_diagram)
	
	fig  = plt.figure(figsize=(5, 5))
	fig.subplots_adjust(wspace=0.4,  hspace=0.6)
	
	colormap   = c.cmap_phase_diagram1 # cmap_gray_cr_pk_gray # c.cmap_white_green_universal, plt.colormaps['jet']# 'winter', etc
	#levels  = np.linspace(0, max_conc, num_levels)
	levels  = np.array([0.0,1.5,2.5,3.5,4.5])
	ax = fig.add_subplot( 1, 1, 1 )
	cs, cb = plot_a_panel(ax, phase_diagram, STG, GluN2B, colormap, levels)
	# cb.set_label('(% Total volume)')
	ax.set_title('Phase diagram')
	
	
	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename_output ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename_output ) ) , dpi=150)
	
	plt.show()
	plt.clf()
	plt.close(fig=fig)
	
	
