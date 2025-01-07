import os, glob, pickle, pprint, copy
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1
from scipy.interpolate import griddata, RegularGridInterpolator


# Extrapolation works, but only for grid data.
def make_grid_using_RegularGridInterpolator(x, y, oZ, mX, mY):
	m_points = np.array([mX.ravel(), mY.ravel()]).T
	f = RegularGridInterpolator((x, y), oZ, bounds_error=False, fill_value=None)
	mZ = f(m_points).reshape(mX.shape)
	return mZ
	
	
def plot_a_panel(ax, oZ, x, y, colormap, levels, draw_border = False, \
	mx_min = 0.0, my_min = 1.0, mx_max = None, my_max = None, ticks=None, margin = 0.05):
	
	# Observed data arrangement
	oX, oY  = np.meshgrid(x, y)
	ox      = np.ravel(oX)
	oy      = np.ravel(oY)
	oz  = np.ravel(oZ.T)
	
		
	# Mesh grids for interpolation
	if mx_max == None:
		mx_max = np.max(x)
	if mx_min == None:
		mx_min = np.min(x)
	if my_max == None:
		my_max = np.max(y)
	if my_min == None:
		my_min = np.min(y)
	
	mx_width = mx_max - mx_min
	my_width = my_max - my_min
	
	
	
	mx = np.linspace(mx_min-mx_width*margin, mx_max+mx_width*margin, 55*4)
	my = np.linspace(my_min-my_width*margin, my_max+my_width*margin, 55*4)
	
	'''
	# Mesh grids for interpolation
	if mx_max == None:
		mx_max = np.max(x)*1.1
	if my_max == None:
		my_max = np.max(y)*1.1
	mx = np.linspace(mx_min, mx_max, 55*4)
	my = np.linspace(my_min, my_max, 55*4)
	'''
	
	mX, mY = np.meshgrid(mx,my)
	
	oZ_panel = copy.deepcopy( oZ )
	oZ_panel[oZ_panel < 0] = 1
	mZ = make_grid_using_RegularGridInterpolator(x, y, oZ_panel, mX, mY)
	
	# Plot
	# colormap.set_bad(color='magenta')
	cs = ax.contourf(mX, mY, mZ, levels=levels, alpha=0.5, \
				cmap= colormap, extend='both' ) # vmin=0, vmax=np.max(levels)
	if draw_border == True:
		ax.contour(cs, colors='k')
	vmin = 0
	vmax = np.max(levels)
	
	
	
	## Exceptional auxiliary line
	'''
	levels = [0.852] # GluN2B_CaMKII
	#levels = [0.407] # _STG_PSD95
	ax.contour(mX, mY, mZ, levels = levels,linestyles = 'dashed', colors='gray',linewidths=2.0)
	'''
	##
	
	
	# Plot scatter 
	ax.scatter(ox, oy, c=oz, cmap=colormap, marker='o', edgecolors='k', s=16, vmin=vmin, vmax=vmax)
	
	
	#ax.set_facecolor("black")
	
	ax.set_xlim( np.min(mx), np.max(mx) )
	ax.set_ylim( np.min(my), np.max(my) )
	
	ax.set_box_aspect(1)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
		
	if ticks is None:
		ticks=np.linspace(vmin, vmax, 4)
	
	divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
	cax = divider.append_axes('right', '5%', pad='3%')
	cb = plt.colorbar(cs, cax=cax, ticks=ticks)
	#cb.ax.set_yticklabels(["{:.2f}".format(i) for i in cb.get_ticks()])
	
	return cs, cb
	
	
	
def plot_a_panel_log(ax, oZ, x, y, colormap, levels, draw_border = False, ticks_level=None, margin = 0.05 ):
	
	# Observed data arrangement
	oX, oY  = np.meshgrid(x, y)
	ox      = np.ravel(oX)
	oy      = np.ravel(oY)
	oz  = np.ravel(oZ.T)
	
	# Mesh grids for interpolation
	mx_max   = np.max(x)
	mx_min   = np.min(x)
	mx_width = mx_max - mx_min
	
	my_max = np.max(y)
	my_min = np.min(y)
	my_width = my_max - my_min
	
	mx = np.linspace(mx_min-mx_width*margin, mx_max+mx_width*margin, 55*4)
	my = np.linspace(my_min-my_width*margin, my_max+my_width*margin, 55*4)
	
	mX, mY = np.meshgrid(mx,my)
	
	oZ_panel = copy.deepcopy( oZ )
	mZ = make_grid_using_RegularGridInterpolator(x, y, oZ_panel, mX, mY)
	
	# Plot
	# colormap.set_bad(color='magenta')
	cs = ax.contourf(mX, mY, mZ, levels=levels, alpha=0.5, \
				cmap= colormap, extend='both' ) # vmin=0, vmax=np.max(levels)
	if draw_border == True:
		ax.contour(cs, colors='k')
	vmin = np.min(levels)
	vmax = np.max(levels)
	
	ax.scatter(ox, oy, c=oz, cmap=colormap, marker='o', edgecolors='k', s=16, vmin=vmin, vmax=vmax)
	
	
	# Overlay exception (not clean).
	'''
	oz_except = (oZ < 1)
	mZ_except = make_grid_using_RegularGridInterpolator(x, y, oz_except, mX, mY)
	mZ_except[mZ_except > 0.5] = 1.0
	mZ_except[mZ_except <= 0.5] = np.nan
	print('np.unique(mZ_except) ' , np.unique(mZ_except) )
	cs = ax.contourf(mX, mY, mZ_except, vmin=0.3, vmax=0.4, cmap='binary' )
	'''
	
	ax.set_xlim( np.min(mx), np.max(mx) )
	ax.set_ylim( np.min(my), np.max(my) )
	
	ax.set_box_aspect(1)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	
	if ticks_level is None:
		ticks_level=np.linspace(vmin, vmax, ticks_level)
	
	divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
	cax = divider.append_axes('right', '5%', pad='3%')
	cb = plt.colorbar(cs, cax=cax, ticks=ticks_level)
	#cb.ax.set_yticklabels(["{:.2f}".format(i) for i in cb.get_ticks()])
	
	return cs, cb
	
	
def plot_a_panel_overlay(ax, oZ, x, y, colormap, levels, draw_border = True, \
	mx_min = 0.0, my_min = 1.0, mx_max = None, my_max = None, margin = 0.05 ):
	
	
	# Observed data arrangement
	oX, oY  = np.meshgrid(x, y)
	ox      = np.ravel(oX)
	oy      = np.ravel(oY)
	oz  = np.ravel(oZ.T)
	
		
	# Mesh grids for interpolation
	if mx_max == None:
		mx_max = np.max(x)
	if mx_min == None:
		mx_min = np.min(x)
	if my_max == None:
		my_max = np.max(y)
	if my_min == None:
		my_min = np.min(y)
	
	mx_width = mx_max - mx_min
	my_width = my_max - my_min
	mx = np.linspace(mx_min-mx_width*margin, mx_max+mx_width*margin, 55*4)
	my = np.linspace(my_min-my_width*margin, my_max+my_width*margin, 55*4)
	
	
	mX, mY = np.meshgrid(mx,my)
	print('x ', x)
	print('y ', y)
	
	mZ = make_grid_using_RegularGridInterpolator(x, y, oZ, mX, mY)
	if draw_border == True:
		ax.contour( mX, mY, mZ, levels=levels, colors='k', zorder = 3)
	mZ[mZ < 0.5] = np.nan
	ax.contourf(mX, mY, mZ, levels=levels, cmap=colormap, zorder = 2) 
	
	
	ids_target = oz > 0
	ax.scatter(ox[ids_target], oy[ids_target], c=oz[ids_target], \
		cmap=colormap, marker='o', edgecolors='k', s=16, \
		vmin=np.min(levels), vmax=np.max(levels), \
		zorder = 4)



def plot_a_panel_intersection(ax, \
		oZ1, level1, \
		oZ2, level2, \
		x, y, draw_border = True, \
		mx_min = 0.0, my_min = 1.0, mx_max = None, my_max = None, ticks=None, margin = 0.05):
	
	# Observed data arrangement
	oX, oY  = np.meshgrid(x, y)
	ox      = np.ravel(oX)
	oy      = np.ravel(oY)
	oz1     = np.ravel(oZ1.T)
	oz2     = np.ravel(oZ2.T)
	
	# Mesh grids for interpolation
	if mx_max == None:
		mx_max = np.max(x)
	if mx_min == None:
		mx_min = np.min(x)
	if my_max == None:
		my_max = np.max(y)
	if my_min == None:
		my_min = np.min(y)
	mx_width = mx_max - mx_min
	my_width = my_max - my_min
	
	# Create meshgrid
	mx = np.linspace(mx_min-mx_width*margin, mx_max+mx_width*margin, 55*4)
	my = np.linspace(my_min-my_width*margin, my_max+my_width*margin, 55*4)
	mX, mY = np.meshgrid(mx,my)
	
	oZ1_panel = copy.deepcopy( oZ1 )
	mZ1 = make_grid_using_RegularGridInterpolator(x, y, oZ1_panel, mX, mY)
	oZ2_panel = copy.deepcopy( oZ2 )
	mZ2 = make_grid_using_RegularGridInterpolator(x, y, oZ2_panel, mX, mY)
	
	
	mZ = (mZ1 > level1)* (mZ2 > level2)
	
	# Plot
	# colormap.set_bad(color='magenta')
	cs = ax.contourf(mX, mY, mZ, \
				extend='both', cmap = 'binary', vmin=0, vmax=10, zorder=-1 ) # cmap= colormap, levels=levels, alpha=0.5, 
	if draw_border == True:
		ax.contour(cs, colors='k',linewidths=0.5)
	#vmin = 0
	#vmax = 10 # np.max(levels)
	'''
	colormap = np.array([[1.0,1.0,1.0], [0.0,0.0,0.0]])
	oZ = (oZ1 > level1) & (oZ2 > level2)
	oZ = oZ.T
	ax.scatter(ox, oy, c=oZ, marker='o', edgecolors='k', s=16, cmap='Greys') # cmap=colormap, , vmin=vmin, vmax=vmax
	'''
	
	#ax.set_facecolor("black")
	
	ax.set_xlim( np.min(mx), np.max(mx) )
	ax.set_ylim( np.min(my), np.max(my) )
	
	ax.set_box_aspect(1)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	
	
	if ticks is None:
		ticks=np.linspace(0, 10, 4)
	
	divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
	cax = divider.append_axes('right', '5%', pad='3%')
	cb = plt.colorbar(cs, cax=cax, ticks=ticks)
	
	
	return
	
	
