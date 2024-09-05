
import os, sys, glob, pickle, pprint, copy

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1

from scipy import ndimage

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
from specification_datasets import SpecDatasets



plt.rcParams.update(p.rc_param)
	
	
def arrange_graph_bar(ax, panel_dx, y0, panel_size_x, panel_size_y):
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax_h, ax_w = ax.bbox.height, ax.bbox.width
	loc = ax.get_position()
	if y0 is None:
		y0 = loc.y0
	ax.set_position([loc.x0+panel_dx, y0, panel_size_x, panel_size_y])
	
	
def plot_a_figure( d ):
	
	# Parameters
	#vmax       = 0.7

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
		axes = []
		column = 1
		ax    = utils.plot_regions_condenstate_from_a_direction_(fig, num_rows, num_columns, row, column, d, transp )
		loc0  = ax.get_position()
		axes.append(ax)
		columns = {'CaMKII':2, 'GluN2B':3, 'PSD95':4,'STG':5}
		axes1 = utils.plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp, pre_rotated=False )
		axes.extend(axes1)
		
		yloc.append(loc0.y0)
		for a in axes:
			loc = a.get_position()
			#print('loc ', loc)
			a.set_position([loc.x0, yloc[-1], panel_size, panel_size])
	
	panel_dx = 0.03
	
	
	# Plot concs in condensates
	column = 9
	targ = 'STG'
	ref  = 'CaMKII'
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*1 )
	utils.plot_conc_condensate_bar(ax, targ, ref, d)
	arrange_graph_bar(ax, panel_dx, yloc[1], panel_size/4, panel_size)
	
	targ = 'CaMKII'
	ref  = 'STG'
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*2 )
	utils.plot_conc_condensate_bar(ax, targ, ref, d)
	arrange_graph_bar(ax, panel_dx, yloc[2], panel_size/4, panel_size)
	
	targs        = ['STG','CaMKII']
	counterparts = ['CaMKII','STG']
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*0 )
	utils.plot_conc_ratio_condensate_bar(ax, targs, counterparts, d)
	arrange_graph_bar(ax, panel_dx, yloc[0], panel_size/4, panel_size )
	
	
	# Average concs in CaMKII/ condensate
	column = 10
	targets_condensate = ['CaMKII','STG']
	for i, t in enumerate(targets_condensate):
		ax = fig.add_subplot( num_rows, num_columns, column+num_columns*i )
		utils.plot_concs_condensate_bar(ax, t, d)
		arrange_graph_bar(ax, panel_dx, yloc[i], panel_size/2, panel_size)
	
	
	# Plot RDF
	ax = fig.add_subplot( num_rows, num_columns, column+num_columns*2 )
	
	
	errorbar= 'shaded' # errorbar='shaded', 'line', or 'final_frame_alone'
	utils.plot_a_rdf( ax, d, errorbar=errorbar, legend=True, ylim = (-0.05,0.5) ) # , ylim = (-0.006,0.46)
	#utils.plot_a_rdf( ax, d, errorbar=errorbar, legend=True, target_molecules = ['All', 'CaMKII', 'GluN2B'] , ylim = (-0.01,1.0) )
	
	
	arrange_graph_bar(ax, panel_dx, yloc[2], panel_size/2, panel_size)
	
	return fig
	
	
	
class PlotProfiles(SpecDatasets):
	def __init__( self ):
		
		pass
		
		
	def run( self ):
		
		# Shared init
		self.dir_imgs = os.path.join(self.dir_imgs_root, 'profiles')
		os.makedirs(self.dir_imgs, exist_ok=True)
		
		for filename in self.filenames_edited:
			fig = self.make_a_figure( filename )
			self.save_a_figure( fig, filename )
		
		
	def make_a_figure( self, filename ):
		
		# Load data
		prefix = filename
		suffix = '' # sigma_2
		print('Target: ', filename)
		d   = utils.load(self.dir_edited_data, prefix, suffix)
		
		# Make figure
		return plot_a_figure(d)
		
		
	def save_a_figure( self, fig, filename ):
		# Save figure
		fig.savefig( os.path.join(self.dir_imgs, '{}.svg'.format( filename ) ) )
		fig.savefig( os.path.join(self.dir_imgs, '{}.png'.format( filename ) ) , dpi=150)
		#plt.show()
		plt.clf()
		plt.close(fig=fig)
	
	
	
if __name__ == '__main__':
	
	obj = PlotProfiles()
	obj.inhibitor()  #  conc_dependence(), valency_length(), valency_length_CG(), boundary_conditions2()
	obj.run()
	
	
