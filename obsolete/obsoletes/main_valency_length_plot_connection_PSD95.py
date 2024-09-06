
import os, sys, glob, pickle, pprint, copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import mpl_toolkits.axes_grid1

from scipy.interpolate import griddata, RegularGridInterpolator
import utils
import colormap as c
import parameters as p


plt.rcParams.update(p.rc_param)

# Extrapolation works, but only for grid data.
def make_grid_using_RegularGridInterpolator(STG, GluN2B, oZ, mX, mY):
	m_points = np.array([mX.ravel(), mY.ravel()]).T
	f = RegularGridInterpolator((STG, GluN2B), oZ, bounds_error=False, fill_value=None)
	mZ = f(m_points).reshape(mX.shape)
	return mZ
	
	
def plot_a_panel(ax, oZ, STG, GluN2B, colormap, levels):
	
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
	
	oZ_panel = copy.deepcopy( oZ )
	oZ_panel[oZ_panel < 0] = 1
	mZ = make_grid_using_RegularGridInterpolator(STG, GluN2B, oZ_panel, mX, mY)
	
	# Plot
	# colormap.set_bad(color='magenta')
	cs = ax.contourf(mX, mY, mZ, levels=levels, alpha=0.5, \
				cmap= colormap, extend='both' ) # vmin=0, vmax=np.max(levels), 
	# ax.contour(cs, colors='k')
	vmin = 0
	vmax = np.max(levels)
	
	ax.scatter(ox, oy, c=oz, cmap=colormap, marker='o', edgecolors='k', s=16, vmin=vmin, vmax=vmax)
	
	
	# Overlay exception (not clean).
	'''
	oz_except = (oZ < 1)
	mZ_except = make_grid_using_RegularGridInterpolator(STG, GluN2B, oz_except, mX, mY)
	mZ_except[mZ_except > 0.5] = 1.0
	mZ_except[mZ_except <= 0.5] = np.nan
	print('np.unique(mZ_except) ' , np.unique(mZ_except) )
	cs = ax.contourf(mX, mY, mZ_except, vmin=0.3, vmax=0.4, cmap='binary' )
	'''
	#ax.set_facecolor("black")
	
	
	ax.set_xlim( np.min(mx), np.max(mx) )
	ax.set_ylim( np.min(my), np.max(my) )
	
	ax.set_box_aspect(1)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	
	divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
	cax = divider.append_axes('right', '5%', pad='3%')
	cb = plt.colorbar(cs, cax=cax, ticks=np.linspace(vmin, vmax, 5))
	#cb.ax.set_yticklabels(["{:.2f}".format(i) for i in cb.get_ticks()])
	
	return cs, cb
	
	
	

if __name__ == '__main__':

	
	# Dataset 1
	species = 'PSD95_connection' # 'STG','GluN2B', 'PSD95','CaMKII', 'PSD95_connection'
	
		
	# Output files
	dir_target =  'valency_length'
	dir_imgs = os.path.join('imgs3', dir_target)
	filename_output = 'phase_diagram'
	os.makedirs(dir_imgs, exist_ok=True)
	suffix = 'connection'
	dir_edited_data  = os.path.join('data3',dir_target)
	
	
	valency = list(range(2,14,2)) 
	length  = [1, 2, 3, 4, 5, 6, 9]
	
	fnames_valency       = { v: str(v).zfill(2) for v in valency }
	fnames_linker_length = {ll: str(i).zfill(3) for i,ll in enumerate(length) }
	
	num_rows		= len( valency )
	num_columns		= len( length )
	
	#
	ave_num_binding_GluN2B = np.zeros([num_columns, num_rows], dtype = 'float')
	ave_num_binding_GluN2B_CaMKII = np.zeros([num_columns, num_rows], dtype = 'float')
	ave_num_binding_STG    = np.zeros([num_columns, num_rows], dtype = 'float')
	ave_num_binding_PSD95_STG= np.zeros([num_columns, num_rows], dtype = 'float')
	ratio_binding_None     = np.zeros([num_columns, num_rows], dtype = 'float')
	ratio_binding_STG      = np.zeros([num_columns, num_rows], dtype = 'float')
	ratio_binding_GluN2B   = np.zeros([num_columns, num_rows], dtype = 'float')
	ratio_binding_Both     = np.zeros([num_columns, num_rows], dtype = 'float')
	#
	for i, v in enumerate(valency):
		for j, ll in enumerate(length):
			
			# Load data
			prefix = fnames_valency[v]+'_'+fnames_linker_length[ll]
			d      = utils.load(dir_edited_data, prefix, suffix)
			
			# sys.exit(0)
			title = prefix # prefix, None
			row    = num_rows-j-1
			column = i+1
			print('column: ', column, ', row: ', row)

			species = 'PSD95'
			ave_num_binding_GluN2B[j,i] = np.average(d[species]['nums']['GluN2B_PSD95'])
			ave_num_binding_STG[j,i] = np.average(d[species]['nums']['STG_PSD95'])

			species = 'CaMKII'
			ave_num_binding_GluN2B_CaMKII[j,i] = np.average(d[species]['nums']['GluN2B_CaMKII'])

			species = 'STG'
			ave_num_binding_PSD95_STG[j,i] = np.average(d[species]['nums']['STG_PSD95'])

			species = 'PSD95_connection'
			total_number = sum(d[species].values())
			ratio_binding_None[j,i]   = d[species]['None'] / total_number
			ratio_binding_STG[j,i]    = d[species]['STG only'] / total_number
			ratio_binding_GluN2B[j,i] = d[species]['GluN2B only'] / total_number
			ratio_binding_Both[j,i]   = d[species]['Both'] / total_number
	
	'''
	data  = ave_num_binding_GluN2B
	title = 'Number of GluN2B bound to one PSD95'
	filename_output = 'num_GluN2B_bound_to_one_PSD95'
	colormap   =  c.cmap_white_purple_universal
	levels  = np.linspace(0,3,10)
	
	data  = ave_num_binding_GluN2B_CaMKII
	title = 'Number of GluN2B bound to one CaMKII'
	filename_output = 'num_GluN2B_bound_to_one_CaMKII'
	colormap   =  c.cmap_white_green_universal
	levels  = np.linspace(0,12,7)
	
	
	data  = ave_num_binding_PSD95_STG
	title = 'Number of PSD95 bound to one STG'
	filename_output = 'num_PSD95_bound_to_one_STG'
	colormap   = c.cmap_white_red_universal
	levels  = np.linspace(0,4,6)
	'''
	
	data  = ave_num_binding_STG
	title = 'Number of STG bound to one PSD95'
	filename_output = 'num_STG_bound_to_one_PSD95'
	colormap   = c.cmap_white_red_universal
	levels  = np.linspace(0,3,6)
	
	#'''
	data  = ratio_binding_Both
	title = 'Ratio of PSD95 bound to both GluN2B and STG'
	filename_output = 'PSD95_bound_to_both_GluN2B_STG'
	colormap   =  plt.colormaps['Greys']
	# levels  = np.linspace(0,0.8,9)
	levels  = np.linspace(0,1.0,11)
	#'''
	
	
	#colormap   =  plt.colormaps['Oranges']
	
	
	fig  = plt.figure(figsize=(5, 5))
	fig.subplots_adjust(wspace=0.4,  hspace=0.6)
	
	
	# cmap_gray_cr_pk_gray # c.cmap_white_green_universal, plt.colormaps['jet']# 'winter', etc
	
	
	ax = fig.add_subplot( 1, 1, 1 )
	cs, cb = plot_a_panel(ax, data, length, valency, colormap, levels)
	ax.set_title( title )
	ax.set_xlabel('Linker length (l.u.)')
	ax.set_ylabel('Valency')
	
	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename_output ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename_output ) ) , dpi=150)
	
	plt.show()
	plt.clf()
	plt.close(fig=fig)
	
	
	