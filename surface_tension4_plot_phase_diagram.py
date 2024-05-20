
import os, sys, glob, pickle, pprint
import numpy as np

import matplotlib.pyplot as plt

import networkx as nx

import lib.utils as utils
import lib.utils_fitting as utils_fitting
import lib.parameters as p
import lib.colormap as c


from surface_tension3_plot import calc_angle_and_distance_to_hub, calc_contraction_force


plt.rcParams.update(p.rc_param)
	
	
	
class MatrixValencyLength():
	
	def __init__( self ):
		# Parameters
		self.sigma = 2
		dir_target = 'CG_valency_length'
		self.dir_edited_data = os.path.join('data4', dir_target)
		self.dir_imgs        = os.path.join('imgs4', dir_target, 'phase_diagram')
		os.makedirs(self.dir_imgs, exist_ok=True)
		
		self.valencies = range(4,14,2)
		self.num_rows		= len( self.valencies )
		self.num_columns	= len( p.lengths )
		
		self.radiuses_condensate = utils.load(self.dir_edited_data, 'radiuses', 'CaMKII')
		
		self.ave_cos             = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		self.pull_force          = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		self.pull_force_per_area = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		
		
		self.radiuses_condensate_matrix = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		for i, valency in enumerate(self.valencies):
			for j, linker_length in enumerate(p.lengths):
				prefix      = p.fnames_valency[valency]+'_'+p.fnames_length[linker_length]
				self.radiuses_condensate_matrix[j,i] = self.radiuses_condensate[prefix]
		
		
	def run_calc( self ):
		
		for i, valency in enumerate(self.valencies):
			for j, max_linker_length in enumerate(p.lengths):
				
				# Load data
				
				prefix      = p.fnames_valency[valency]+'_'+p.fnames_length[max_linker_length]
				print('Target file: ', prefix)
				d           = utils.load(self.dir_edited_data, prefix, 'connectivity_graph')
				
				radius_condensate = self.radiuses_condensate[prefix]
				
				angles, distances_to_hub = calc_angle_and_distance_to_hub(d)
				
				ave_cos, contraction_force, contraction_force_per_area = \
						calc_contraction_force(angles, distances_to_hub, max_linker_length, radius_condensate)
				
				self.ave_cos[j,i]    = ave_cos
				self.pull_force[j,i] = contraction_force
				self.pull_force_per_area[j,i] = contraction_force_per_area
		
		
	def prepare_plot( self, title ):
		# self.fig  = plt.figure(figsize=(5, 5))
		self.fig  = plt.figure(figsize=(3.5, 3.5))
		self.fig.subplots_adjust(wspace=0.4,  hspace=0.6)
		ax = self.fig.add_subplot( 1, 1, 1 )
		ax.set_title(title)
		ax.set_xlabel('Linker length (l.u.)')
		ax.set_ylabel('Valency')
		return ax
		
	def save_data( self ):
		
		d = {}
		d['ave_cos']             = self.ave_cos
		d['pull_force']          = self.pull_force
		d['pull_force_per_area'] = self.pull_force_per_area
		
		prefix = 'Surface_tension'
		suffix = 'matrix'
		utils.save(self.dir_edited_data, prefix, suffix, d)
		
	def load_data( self ):
		
		prefix = 'Surface_tension'
		suffix = 'matrix'
		d = utils.load(self.dir_edited_data, prefix, suffix)
		
		self.ave_cos             = d['ave_cos']
		self.pull_force          = d['pull_force']
		self.pull_force_per_area = d['pull_force_per_area']
		
	def plot_phase_diagrams( self ):
		
		title    = 'Radiuses'
		filename = 'Radiuses'
		colormap =  plt.colormaps['Greys'] #  plt.get_cmap plt.colormaps
		levels   = np.linspace(0,36,35)
		ax = self.prepare_plot(title)
		cs, cb = utils.plot_a_panel(ax, self.radiuses_condensate_matrix, p.lengths, self.valencies, colormap, levels)
		self.save_a_fig( filename )
		
		title    = 'Cos similarity'
		filename = 'Cos_similarity'
		colormap =  c.cmap_white_green_universal
		levels   = np.linspace(0,1,10)
		ax = self.prepare_plot(title)
		cs, cb = utils.plot_a_panel(ax, self.ave_cos, p.lengths, self.valencies, colormap, levels)
		self.save_a_fig( filename )
		
		title    = 'Contraction force'
		filename = 'Contraction_force'
		colormap = c.cmap_white_red_universal
		levels   = np.linspace(0,3000,10)
		ax = self.prepare_plot(title)
		cs, cb = utils.plot_a_panel(ax, self.pull_force, p.lengths, self.valencies, colormap, levels)
		self.save_a_fig( filename )
		
		title    = 'Contraction force per area'
		filename = 'Contraction_force_per_area'
		colormap = plt.colormaps['Greys'] #  plt.get_cmap plt.colormaps
		levels   = np.linspace(0,0.6,10)
		ax = self.prepare_plot(title)
		cs, cb = utils.plot_a_panel(ax, self.pull_force_per_area, p.lengths, self.valencies, colormap, levels)
		self.save_a_fig( filename )
		
		
	def save_a_fig( self, filename ):
		
		self.fig.savefig( os.path.join(self.dir_imgs, filename + '.svg' ) )
		self.fig.savefig( os.path.join(self.dir_imgs, filename + '.png' ) , dpi=150)
		plt.show()
		#plt.clf()
		#plt.close(fig=self.fig)
		
		
	def plot_logistic_regression( self ):
		
		filename = 'logistic_regression'
		
		# 0: Partial engulfment
		# 1: PIPS
		# 2: Others
		
		# 1, 2, 3, 4, 5, 6, 9
		y = np.array( [\
			[ 1, 1, 1, 0,  0, 0, 0], # 12
			[ 1, 1, 1, 0,  0, 0, 0], # 10
			[ 1, 1, 1, 0,  0, 0, 0], # 8
			[ 1, 1, 1, 0,  0, 0, 0], # 6
			]).T
		x = self.pull_force_per_area[:,:4]
		
		print('y')
		print(y)
		print('x')
		print(x)
		
		
		x = x.reshape(-1, 1)
		y = y.reshape(-1, 1)

		model = utils_fitting.logistic_regression(x,y)


		self.fig  = plt.figure(figsize=(3.5, 3.5))
		self.fig.subplots_adjust(wspace=0.4,  hspace=0.4)
		ax = self.fig.add_subplot( 1, 1, 1 )
		ax.set_xlabel('Contraction force per area')
		xmax = np.max(x)*1.1
		xx = np.linspace(0.0,xmax,40).reshape(-1, 1)
		ax.plot([0, xmax], [0,0], ':', color = (0.5,0.5,0.5))
		ax.plot([0, xmax], [1,1], ':', color = (0.5,0.5,0.5))
		ax.plot(x, y,'o', \
			markersize = 4, \
			color = 'k', \
			markerfacecolor = 'k' )
		ax.plot(xx, model.predict_proba(xx)[:,1], '-', color = (0.5,0.5,0.5))
		ax.set_ylim([-0.1,1.1])
		ax.set_xlim([   0,xmax])
		ax.set_yticks([0,1])
		ax.set_yticklabels(['Partial \n engulfment','PIPS'])
		
		self.save_a_fig( filename )
		
		
if __name__ == '__main__':
	
	graph = MatrixValencyLength()
	#graph.run_calc()
	#graph.save_data()
	graph.load_data()
	#graph.plot_phase_diagrams()
	graph.plot_logistic_regression()
