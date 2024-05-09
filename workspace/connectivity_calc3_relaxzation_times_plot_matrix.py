
# https://stackoverflow.com/questions/59821151/plot-the-dendrogram-of-communities-found-by-networkx-girvan-newman-algorithm
# https://qiita.com/s-wakaba/items/a93f03f27137cff4a26c

import os, sys, glob, pprint,itertools
import numpy as np


# from itertools import chain, combinations
import matplotlib.pyplot as plt

from cycler import cycler
from scipy.optimize import curve_fit

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph

plt.rcParams.update(p.rc_param)
	
	
	
class PlotRelaxzationTimes():
	
	def __init__( self ):
		pass
	
	def plot_a_graph( self, row, column, d, title ):
		
		time_points, profs, taus, params, colors = utils_graph.get_profiles(d)
		
		# Plot figure
		ax = self.fig.add_subplot( self.num_rows, self.num_columns, row*self.num_columns+column )
		utils_graph.prep_fig(ax, time_points)
		for col, prof, param in zip( colors, profs, params ):
			ax.plot(time_points, prof, 'o', markersize=3, color= col )
			ax.plot(time_points, utils_graph.func_exponential(time_points,param[0],param[1],param[2]), '-', color= col )
		ax.set_title(title+', Tau: {:.4f}'.format(np.median(taus)))
		
		return taus
		
		
class MatrixValencyLength():
	
	def load_data(self, prefix_load_):
		data = {}
		for time_frame in range(self.num_time_frames):
			prefix_load = '{}_{}'.format(time_frame, prefix_load_)
			data[time_frame] = utils.load(self.dir_edited_data, prefix_load, self.suffix_load)
		return data

	def __init__( self ):
		
		self.valencies = list(range(4,14,2))
		self.lengths   = list(range(0,7))
		#self.valencies = list(range(8,14,2))
		#self.lengths   = list(range(0,3))
		self.num_time_frames = 80
		
		dir_target  = 'small_colony2'
		#suffix_load = 'greedy_modularity'
		self.suffix_load = 'louvain'
		self.basename    = 'louvain_matrix'
		
		# Shared init
		self.num_rows        = len( self.valencies )
		self.num_columns     = len( self.lengths )	
		self.dir_edited_data = os.path.join('data4', dir_target)
		self.dir_imgs = os.path.join('imgs4', dir_target, 'connectivity_relaxzation_time')
		os.makedirs(self.dir_imgs, exist_ok=True)
		
	def run( self ):
		
		vals = {}
		self.fig  = plt.figure(figsize=(10, 10), tight_layout=True)
		
		for i, v in enumerate(self.valencies):
			for j, l in enumerate(self.lengths):
				#
				row    = self.num_rows-i-1
				column = j+1
				# Load data
				prefix = str(v).zfill(2)+'_'+str(l).zfill(3)
				print('Target file: ', prefix, ', column: ', column, ', row: ', row)
				d = self.load_data(prefix)
				title = prefix
				vals[prefix] = self.plot_a_graph(row, column, d, title)
		return vals
		
		
	def save( self ):
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.svg' ) )
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.png' ), dpi=150 )
		plt.show()
		plt.clf()
		plt.close(fig=self.fig)
		
		
		
class PlotRelaxzationTimesMatrixValencyLength(PlotRelaxzationTimes, MatrixValencyLength):
	def __init__( self ):
		PlotRelaxzationTimes.__init__(self)
		MatrixValencyLength.__init__(self)
	
	
if __name__ == '__main__':
	
	relaxzation_times = PlotRelaxzationTimesMatrixValencyLength()
	values = relaxzation_times.run()
	relaxzation_times.save()
	
	dir_edited_data = relaxzation_times.dir_edited_data
	prefix = 'relaxzation_times'
	suffix = 'matrix'
	utils.save(dir_edited_data, prefix, suffix, values)

