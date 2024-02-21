
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
	
	
def plot_a_panel(ax, oZ, x, y, colormap1, levels, oZ_except, colormap2, method = 'RegularGridInterpolator'):
	
	# STG - x , GluN2B - y
	# linker_length - x, valency - y
	
	# Observed data arrangement
	oX, oY  = np.meshgrid(x, y)
	ox      = np.ravel(oX)
	oy      = np.ravel(oY)
	oz  = np.ravel(oZ.T)
	
	# Mesh grids for interpolation
	mx_max = np.max(x)
	my_max = np.max(y)
	
	mx = np.linspace(0.0, mx_max*1.1, 55*4)
	my = np.linspace(3.0, my_max*1.1, 55*4)
	
	mX, mY = np.meshgrid(mx,my)
	
	#print('oz_graph.shape ', oz_graph.shape )
	
	oZ_panel = copy.deepcopy( oZ )
	oZ_panel[oZ_panel < 0] = 1
	mZ = make_grid_using_RegularGridInterpolator(x, y, oZ_panel, mX, mY)
	mZ[mZ < np.min(oZ)] = np.min(oZ)
	
	# Plot
	cs = ax.contourf(mX, mY, mZ, levels=levels, alpha=0.5, \
				cmap= colormap1 ) # vmin=0, vmax=np.max(levels), 
	ax.contour(cs, colors='k')
	ax.scatter(ox, oy, c=oz, cmap=colormap1, marker='o', edgecolors='k', s=16, \
		vmin=np.min(levels), vmax=np.max(levels), \
		zorder = 1)
	
	
	# Overlay exception (not clean).
	#'''
	levels = [-0.5, 0.5, 1.5]
	mZ_except = make_grid_using_RegularGridInterpolator(x, y, oZ_except, mX, mY)
	ax.contour( mX, mY, mZ_except, levels=levels, colors='k', zorder = 3)
	mZ_except[mZ_except < 0.5] = np.nan
	ax.contourf(mX, mY, mZ_except, levels=levels, cmap=colormap2, zorder = 2) 
	
	oz_except  = np.ravel(oZ_except.T)
	ids_target = oz_except > 0
	ax.scatter(ox[ids_target], oy[ids_target], c=oz[ids_target], \
		cmap=colormap2, marker='o', edgecolors='k', s=16, \
		vmin=np.min(levels), vmax=np.max(levels), \
		zorder = 4)
	#'''
	
	# mZ_except[mZ_except >= 0.5] = 1.0
	# oz_except = (oZ_except < 1)
	
	ax.set_xlim( np.min(mx), np.max(mx) )
	ax.set_ylim( np.min(my), np.max(my) )
	
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
	target_dir = 'valency_linker_length'
	dir_imgs = os.path.join('imgs', target_dir,'phase_diagram')
	filename_output = 'phase_diagram'
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	#STG    = [1, 500, 1000, 1500, 2000, 2500, 3000] 
	#GluN2B = [20000, 16000, 12000, 8000, 6000, 4000, 2000, 1000, 500, 1]
	
	linker_length  = [1, 3, 5, 7, 9]
	valency = [12, 10, 8, 6, 4] # [4, 6, 8, 10, 12] 
	
	#'''
	print('Valency')
	print(valency)
	print('Linker_length')
	print(linker_length)
	#'''
	
	# 1: PIPS
	# 2: Partial engulfment
	# 3: iPIPS
	# 4: Homogeneous LLPS
	
	phase_diagram = [\
		[ 1, 1, 2, 2, 2], # 12
		[ 1, 2, 2, 2, 2], # 10
		[ 1, 1, 2, 2, 2], # 8
		[ 1, 1, 2, 2, 3], # 6
		[ 1, 1, 2, 2, 3]  # 4
		]
	
	homogeneous_LLPS = [\
		[ 0, 0, 0, 0, 0], # 12
		[ 0, 0, 0, 0, 0], # 10
		[ 0, 0, 0, 0, 0], # 8
		[ 0, 0, 0, 0, 0], # 6
		[ 0, 0, 1, 1, 1]  # 4
		]
	
	phase_diagram = np.array(phase_diagram).T
	levels        = np.array([0.5, 1.5, 2.5, 3.5, 4.5])
	homogeneous_LLPS = np.array(homogeneous_LLPS).T
	
	
	
	max_conc      = np.max(phase_diagram)
	
	#fig  = plt.figure(figsize=(5, 5))
	
	fig  = plt.figure(figsize=(4, 4))
	fig.subplots_adjust(wspace=0.4,  hspace=0.6)
	
	colormap1  = c.cmap_phase_diagram2
	colormap2  = c.cmap_phase_diagram3
	# cmap_gray_cr_pk_gray # c.cmap_white_green_universal, plt.colormaps['jet']# 'winter', etc
	# levels  = np.linspace(0, max_conc, num_levels)
	
	ax = fig.add_subplot( 1, 1, 1 )
	cs, cb = plot_a_panel(ax, phase_diagram, linker_length, valency, colormap1, levels, homogeneous_LLPS, colormap2)
	# STG: linker length, GluN2B: valency
	# cb.set_label('(% Total volume)')
	ax.set_title('Phase diagram')
	ax.set_xlabel('Linker length (l.u.)')
	ax.set_ylabel('Valency')
	
	
	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename_output ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename_output ) ) , dpi=150)
	
	plt.show()
	plt.clf()
	plt.close(fig=fig)
	
	
