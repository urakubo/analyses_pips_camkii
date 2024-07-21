
import os, sys, glob, pickle, pprint
import numpy as np

import matplotlib.pyplot as plt

import networkx as nx

import lib.utils as utils
import lib.utils_panel as utils_panel

import lib.utils_fitting as utils_fitting
import lib.utils_surface_tension as utils_surface_tension
import lib.parameters as p
import lib.colormap as c

from specification_datasets import SpecDatasets


plt.rcParams.update(p.rc_param)
	
	
	
class PlotSurfaceTensionPhaseDiagramValencyLength(SpecDatasets):
	
	def __init__( self ):
		
		pass
		
		
	def reflect_spec( self ):
		
		self.dir_imgs = os.path.join(self.dir_imgs_root, 'phase_diagram')
		os.makedirs(self.dir_imgs, exist_ok=True)
		
		
		self.num_rows		= len( self.surface_tension_valencies )
		self.num_columns	= len( self.surface_tension_lengths )
		
		self.ave_cos             = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		self.pull_force          = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		self.surface_tensions    = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		self.pull_force_per_area = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		
		
		self.radii_condensate = utils.load(self.dir_edited_data, 'radiuses', 'CaMKII')
		self.radii_condensate_matrix = np.zeros([self.num_columns, self.num_rows], dtype = 'float')
		for i, v in enumerate( self.surface_tension_valencies ):
			for j, l in enumerate( self.surface_tension_lengths ):
				prefix = self.filename_edited_matrix(v, l)
				self.radii_condensate_matrix[j,i] = self.radii_condensate[prefix]
		
		
	def run_calc( self ):
		
		for i, v in enumerate( self.surface_tension_valencies ):
			for j, l in enumerate(self.surface_tension_lengths ):
				
				
				filename_prefix    = self.filename_edited_matrix(v, l)
				d                  = utils.load(self.dir_edited_data, filename_prefix, 'connectivity_graph')
				real_linker_length = self.real_linker_length_from_filename[filename_prefix]
				print('Target file: ', filename_prefix)
				
				radiuse_condensate = self.radii_condensate[filename_prefix]
				
				angles_from_condensate_center, angles_from_hub, distances_to_hub = utils_surface_tension.calc_angle_and_distance_to_hub(d)
				
				ave_cos, contraction_force, surface_tension, contraction_force_per_area = \
						utils_surface_tension.calc_contraction_force(angles_from_condensate_center,\
						distances_to_hub, real_linker_length, radiuse_condensate)
				
				self.ave_cos[j,i]    = ave_cos
				self.pull_force[j,i] = contraction_force
				self.surface_tensions[j,i] = surface_tension
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
		d['surface_tensions']    = self.surface_tensions
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
		self.surface_tensions    = d['surface_tensions']
		self.pull_force_per_area = d['pull_force_per_area']
		
	def plot_phase_diagrams( self ):
		
		title    = 'Radii'
		filename = 'Radii'
		colormap =  plt.colormaps['Greys'] #  plt.get_cmap plt.colormaps
		levels   = np.linspace(0,36,35)
		ax = self.prepare_plot(title)
		cs, cb = utils_panel.plot_a_panel(ax, self.radii_condensate_matrix,\
			self.surface_tension_real_lengths, \
			self.surface_tension_valencies, \
			colormap, levels, \
			mx_min = None, my_min = None \
			)
		self.save_a_fig( filename )
		
		title    = 'Average_cos'
		filename = 'Average_cos'
		colormap =  c.cmap_white_green_universal
		levels   = np.linspace(0,1,10)
		ax = self.prepare_plot(title)
		cs, cb = utils_panel.plot_a_panel(ax, self.ave_cos, \
			self.surface_tension_real_lengths, \
			self.surface_tension_valencies, \
			colormap, levels, \
			mx_min = None, my_min = None \
			)
		self.save_a_fig( filename )
		
		title    = 'Contraction force'
		filename = 'Contraction_force'
		colormap = c.cmap_white_red_universal
		levels   = np.linspace(0,3000,10)
		ax = self.prepare_plot(title)
		cs, cb = utils_panel.plot_a_panel(ax, self.pull_force, \
			self.surface_tension_real_lengths, \
			self.surface_tension_valencies, \
			colormap, levels, \
			mx_min = None, my_min = None \
			)
		self.save_a_fig( filename )
		
		title    = 'Contraction force per area'
		filename = 'Contraction_force_per_area'
		colormap = plt.colormaps['Greys'] #  plt.get_cmap plt.colormaps
		levels   = np.linspace(0,0.6,10)
		ax = self.prepare_plot(title)
		cs, cb = utils_panel.plot_a_panel(ax, self.pull_force_per_area, \
			self.surface_tension_real_lengths, \
			self.surface_tension_valencies, \
			colormap, levels, \
			mx_min = None, my_min = None \
			)
		self.save_a_fig( filename )
		
		title    = 'Surface tension'
		filename = 'Surface_tension'
		colormap = plt.colormaps['Greys'] #  plt.get_cmap plt.colormaps
		levels   = np.linspace(0,25,8)
		ax = self.prepare_plot(title)
		ticks = [0, 5, 10, 15, 20, 25]
		cs, cb = utils_panel.plot_a_panel(ax, self.surface_tensions, \
			self.surface_tension_real_lengths, \
			self.surface_tension_valencies, \
			colormap, levels, ticks=ticks,\
			mx_min = None, my_min = None \
			)
		self.save_a_fig( filename )
		
		
		
	def save_a_fig( self, filename ):
		
		self.fig.savefig( os.path.join(self.dir_imgs, filename + '.svg' ) )
		self.fig.savefig( os.path.join(self.dir_imgs, filename + '.png' ) , dpi=150)
		plt.show()
		#plt.clf()
		#plt.close(fig=self.fig)
		
		
	def plot_logistic_regression( self ):
		
		
		filename = 'logistic_regression_surface_tension'
		x = self.surface_tensions
		
		#filename = 'regression_cos_similarity'
		#x = self.ave_cos
		
		#filename = 'regression_contraction_force'
		#x = self.pull_force
		
		x = np.fliplr(x)
		
		print('x')
		print(x)		
		
		x = x[:,:4]
		
		# 0: Partial engulfment
		# 1: PIPS
		# 2: Others
		
		# 1, 2, 3, 4, 5, 6, 9
		y = np.array( [\
			[ 1, 1, 1, 0,  0, 0, 0], # 12
			[ 1, 1, 1, 0,  0, 0, 0], # 10
			[ 1, 1, 1, 1,  0, 0, 0], # 8
			[ 1, 1, 1, 1,  0, 0, 0], # 6
			]).T
		
		# 2, 3, 4, 5, 6, 9
		y = np.array( [\
			[ 1, 1, 0,  0, 0, 0], # 12
			[ 1, 1, 0,  0, 0, 0], # 10
			[ 1, 1, 1,  0, 0, 0], # 8
			[ 1, 1, 1,  0, 0, 0], # 6
			]).T



		#'''
		print('y')
		print(y)
		print('x')
		print(x)
		#'''
		
		x = x.reshape(-1, 1)
		y = y.reshape(-1, 1)

		model = utils_fitting.logistic_regression(x,y)


		self.fig  = plt.figure(figsize=(2.6, 2.0))
		self.fig.subplots_adjust(wspace=0.4,  hspace=0.4)
		ax = self.fig.add_subplot( 1, 1, 1 )
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.set_xlabel('Surface tension')
		#xmax = np.max(x)*1.1
		xmax = 30
		ax.set_xticks([0,10,20,30])
		xx = np.linspace(0.0,xmax,40).reshape(-1, 1)
		ax.plot([0, xmax], [0,0], ':', color = (0.5,0.5,0.5))
		ax.plot([0, xmax], [1,1], ':', color = (0.5,0.5,0.5))
		ax.plot(x, y,'o', \
			markersize = 4, \
			color = 'k', \
			markerfacecolor = 'k' )
		ax.plot(xx, model.predict_proba(xx)[:,1], '-', color = (0.5,0.5,0.5))
		ax.set_ylim([-0.2,1.2])
		ax.set_xlim([   0,xmax])
		ax.set_yticks([0,1])
		ax.set_yticklabels(['Partial \n engulfment','PIPS'])
		
		self.save_a_fig( filename )
		
		
if __name__ == '__main__':
	
	graph = PlotSurfaceTensionPhaseDiagramValencyLength()
	graph.CG_valency_length()
	graph.reflect_spec()
	#graph.run_calc()
	#graph.save_data()
	graph.load_data()
	graph.plot_phase_diagrams()
	graph.plot_logistic_regression()
