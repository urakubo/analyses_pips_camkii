
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
	
	
if __name__ == '__main__':
	
		
	# Overlay exception (not clean).
	'''
	oz_except = (oZ < 1)
	mZ_except = make_grid_using_RegularGridInterpolator(STG, GluN2B, oz_except, mX, mY)
	mZ_except[mZ_except > 0.5] = 1.0
	mZ_except[mZ_except <= 0.5] = np.nan
	print('np.unique(mZ_except) ' , np.unique(mZ_except) )
	cs = ax.contourf(mX, mY, mZ_except, vmin=0.3, vmax=0.4, cmap='binary' )
	'''
	
	# Output files
	dir_imgs = os.path.join('imgs3', 'conc_dependence','phase_diagram')
	filename_output = 'phase_diagram_conc_dependence'
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	
	STG    = [540, 1620, 2160, 2700, 3240, 4520] 
	GluN2B = ([-570, 570, 1080, 4320, 6480, 8640, 10800, 12960, 17280])
	GluN2B.reverse()
	
	volume = np.prod(p.space_np)
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
	[ 1, 1, 1, 1,   1, 1], # 17280
	[ 1, 1, 1, 1,   2, 2], # 12960
	[ 1, 2, 2, 2,   2, 2], # 10800
	[ 2, 2, 2, 2,   2, 2], # 8640
	[ 2, 3, 2, 2,   2, 2], # 6480
	[ 3, 3, 3, 3,   3, 3], # 4320
	[ 4, 4, 4, 3,   3, 3], # 1080
	[ 4, 4, 4, 4,   4, 4], # 570
	[ 4, 4, 4, 4,   4, 4]] # -570
	
	
	
	phase_diagram = np.array(phase_diagram).T
	
	
	fig  = plt.figure(figsize=(5, 5))
	fig.subplots_adjust(wspace=0.4,  hspace=0.6)
	
	colormap   = c.cmap_phase_diagram1
	# cmap_gray_cr_pk_gray # c.cmap_white_green_universal, plt.colormaps['jet']# 'winter', etc
	
	
	levels  = np.array([0.0,1.5,2.5,3.5,4.5])
	ax = fig.add_subplot( 1, 1, 1 )
	cs, cb = utils.plot_a_panel(ax, phase_diagram, STG, GluN2B, colormap, levels)
	ax.set_title('Phase diagram')
	ax.set_xlabel('STG (beads / voxel) x 10-3')
	ax.set_ylabel('GluN2B (beads / voxel) x 10-2')
	
	
	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename_output ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename_output ) ) , dpi=150)
	
	plt.show()
	plt.clf()
	plt.close(fig=fig)
	
	
