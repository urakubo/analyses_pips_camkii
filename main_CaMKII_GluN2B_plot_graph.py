
import os, sys, glob, pickle, pprint, copy

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1

import utils
import parameters as p
import colormap as c
from scipy import ndimage


plt.rcParams.update(p.rc_param)
	
	
def arrange_graph_bar(ax, panel_dx, y0, panel_size_x, panel_size_y):
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax_h, ax_w = ax.bbox.height, ax.bbox.width
	loc = ax.get_position()
	if y0 is None:
		y0 = loc.y0
	ax.set_position([loc.x0+panel_dx, y0, panel_size_x, panel_size_y])
	
	
def plot_bar(ax, data, color, ylim, ylegend, width=0.5):
	ax.bar(*zip(*data.items()), width=width, color=color)
	ax.set_ylabel(ylegend)
	ax.set_ylim(ylim)
	
	
def make_a_figure( d, data ):
	
	# Parameters
	vmax       = 0.7

	# Figure specification
	transps = [(0,1,2),(1,2,0),(2,0,1)]
	titles  = [ 'yz', 'xy', 'zx']
	num_columns = 10
	num_rows    = 3
	fig = plt.figure(figsize=(20, 8)) # , tight_layout=False
	
	left, right, bottom, top = 0.0, 0.95, 0.10, 0.99
	wspace, hspace = 0.2, 0.1
	plt.subplots_adjust(left, bottom, right, top, wspace, hspace)
	fig.suptitle( "Time step: {}, MC step: {:e}".format( d['time_frame'], d['mc_step'] ) )
	
	panel_size = 0.15
	
	# Plot profiles
	yloc = []
	for row,(title, transp) in enumerate(zip(titles, transps)):
		columns = {'CaMKII':1, 'GluN2B':2}
		axes = utils.plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp, pre_rotated=False )
		loc0  = axes[0].get_position()
		yloc.append(loc0.y0)
		for a in axes:
			loc = a.get_position()
			#print('loc ', loc)
			a.set_position([loc.x0, yloc[-1], panel_size, panel_size])
	
	panel_dx = 0.03
	
	# Plot concs in condensates
	r = 0.4
	color_gray = (r,r,r)
	color2 = c.cmap_universal_ratio['CaMKII'] # c.cmap_universal_ratio_light['CaMKII']
	
	
	column = 4
	ylim = [0,0.3]
	xlim = [-0.5,1.5]
	ylegend =  '(ratio)'
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*1 )
	plot_bar(ax, data['binding'], [color_gray, color2], ylim, ylegend)
	ax.set_xlim(xlim)
	yy = data['binding']['All']
	ax.plot(xlim, [yy,yy],':',color = color_gray)
	ax.set_title('Unconnected subunit \n of CaMKII')
	arrange_graph_bar(ax, panel_dx, yloc[1], panel_size/4, panel_size)
	#ax.tick_params(axis='x', rotation=45)
	
	title = 'Clustering \n coefficient'
	column = 5
	ylim = [0,0.12]
	ylegend =  '(1)'
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*1 )
	plot_bar(ax, {title :data['clustering']}, color_gray, ylim, ylegend, width=0.25)
	ax.set_xlim(-0.5,0.5)
	ax.set_title('Clustering \n coefficient')
	arrange_graph_bar(ax, panel_dx, yloc[1], panel_size/4, panel_size)
	
	return fig
	
	
if __name__ == '__main__':
	
	
	# Short GluN2B length 1
	dir_target  = 'small_colony'
	filenames_edited_data = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(3) for id_f in range(10)]
	
	ids = [8, 9, 18, 19]
	data = {filenames_edited_data[i]: [] for i in ids}
	
	
	data['00_008'] = {'binding': {'All': 0.185833, 'Surface':0.23968565815324164} , 'clustering': 0.111083 , 'average_shortest_path': 3.59527568}
	data['00_009'] = {'binding': {'All': 0.178   , 'Surface':0.2310405643738977}   , 'clustering': 0.1019189, 'average_shortest_path': 3.8130020}
	data['01_008'] = {'binding': {'All': 0.21708 , 'Surface':0.18005390835579516}     , 'clustering': 0.0215207, 'average_shortest_path': 2.9301503}
	data['01_009'] = {'binding': {'All': 0.208   , 'Surface':0.17779868297271872}      , 'clustering': 0.0147354, 'average_shortest_path': 3.0247695}
	
	
	
	# Shared init
	dir_edited_data	= os.path.join('data3', dir_target)
	dir_imgs = os.path.join('imgs3', dir_target,'profiles2')
	sigma = 2
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	for i in ids:
		filename = filenames_edited_data[i]
		
		# Load data
		prefix = filename
		suffix = 'sigma_{}'.format(sigma)
		print('Target: {}'.format(filename))
		d   = utils.load(dir_edited_data, prefix, suffix)
		
		# Make figure
		fig = make_a_figure(d, data[prefix])
		
		# Save figure
		fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename ) ) )
		fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename ) ) , dpi=150)
		#plt.show()
		plt.clf()
		plt.close(fig=fig)

