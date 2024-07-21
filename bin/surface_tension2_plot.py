
import os, sys, glob, pickle, pprint
import numpy as np
import math

import matplotlib
import matplotlib.pyplot as plt

import networkx as nx

import lib.parameters as p
import lib.colormap as c
import lib.utils as utils
import lib.utils_surface_tension as utils_surface_tension


from specification_datasets import SpecDatasets



plt.rcParams.update(p.rc_param)

class PlotSurfaceTension(SpecDatasets):
		
	def __init__( self ):
		
		self.basename = 'radius_CaMKII'
		
		
	def apply_specification( self ):
		
		self.dir_imgs = os.path.join(self.dir_imgs_root, 'surface_tension')
		os.makedirs(self.dir_imgs, exist_ok=True)
		self.radii_condensate = utils.load(self.dir_edited_data, 'radiuses', 'CaMKII')
		
		
	def inspect_targets( self ):
		
		print('Target folder : {}\n'.format(self.dir_edited_data) )
		for filename, real_linker_length in self.real_linker_length_from_filename.items() :
			print('filename: {}, real linker length: {}'.format(filename, real_linker_length))
		print()
		
		
	def show_polar_graphs( self, targets = ['12_002','12_006'], mode = 'angle_from_condensate_center', mode_surrogate = True ):
		
		for target in targets:
			print(target)
			d = utils.load(self.dir_edited_data, target, 'connectivity_graph')
			real_linker_length = self.real_linker_length_from_filename[target]
			radius_condensate = self.radii_condensate[target]
			
			angles_from_condensate_center, angles_from_hub, distances_to_hub = \
				utils_surface_tension.calc_angle_and_distance_to_hub(d, surrogate = mode_surrogate)
			
			if mode == 'angle_from_condensate_center':
				filename_prefix = target + '_' + mode
				angles = angles_from_condensate_center
				utils_surface_tension.plot_polar_histogram(angles, distances_to_hub, real_linker_length, \
					self.dir_imgs, filename_prefix, step_color= 'Gray', center_direction = True ) # 'Greens'
				utils_surface_tension.plot_polar_scatter(angles, distances_to_hub, real_linker_length, \
					self.dir_imgs, filename_prefix, center_direction = True )

			elif mode == 'angle_from_hub':
				filename_prefix = target + '_' + mode
				angles = angles_from_hub
				utils_surface_tension.plot_polar_histogram(angles, distances_to_hub, real_linker_length, \
					self.dir_imgs, filename_prefix, cmap='Greens', step_color=c.light_green_universal_ratio,\
					center_direction = False) # c.cmap_white_green_universal 
				utils_surface_tension.plot_polar_scatter(angles, distances_to_hub, real_linker_length, \
					self.dir_imgs, filename_prefix, color=c.cmap_universal_ratio['CaMKII'], \
					center_direction = False)
		
		
	def single_run( self, filename_prefix = '12_002'):
		
		print(filename_prefix)
		
		# Load data CaMKII.
		d = utils.load(self.dir_edited_data, filename_prefix, 'connectivity_graph')
		real_linker_length = self.real_linker_length_from_filename[filename_prefix]
		radius_condensate = self.radii_condensate[filename_prefix]
		
		angles_from_condensate_center, angles_from_hub, distances_to_hub = utils_surface_tension.calc_angle_and_distance_to_hub(d)
		
		ave_cos, contraction_force, surface_tension, contraction_force_per_area = \
			utils_surface_tension.calc_contraction_force(angles_from_condensate_center, distances_to_hub, real_linker_length, radius_condensate)
		
		#print(' angles_from_hub ',  angles_from_hub )
		#print(' distances_to_hub ',  distances_to_hub )
		
		'''
		ave_cos, contraction_force, surface_tension, contraction_force_per_area = \
			utils_surface_tension.calc_contraction_force(angles_from_hub, distances_to_hub, real_linker_length, radius_condensate)
		'''
		
		
		return ave_cos, contraction_force, surface_tension, contraction_force_per_area
		
		
	def multiple_run_plot( self, filenames_prefix = [ '12_002', '12_006'] ):
		
		aves_cos = {}
		surface_tensions = {}
		pull_forces = {}
		pull_forces_per_area = {}
		for filename_prefix in filenames_prefix:
			
			ave_cos, contraction_force, surface_tension, contraction_force_per_area = \
				self.single_run( filename_prefix )
			
			l = str( self.real_linker_length_from_filename[filename_prefix] )
			surface_tensions[l]		= surface_tension
			aves_cos[l]				= ave_cos
			pull_forces[l]			= contraction_force
			pull_forces_per_area[l]	= contraction_force_per_area
		
		
		fig = utils_surface_tension.plot_graphs(aves_cos, pull_forces, surface_tensions, pull_forces_per_area)
		fig.savefig( os.path.join(self.dir_imgs, 'Surface_tension.svg' ) )
		fig.savefig( os.path.join(self.dir_imgs, 'Surface_tension.png' ) , dpi=150)
		plt.show()
		
		
		
if __name__ == '__main__':
	
	
	obj = PlotSurfaceTension()
	obj.CG_valency_length()
	obj.reflect_spec()
	obj.inspect_targets()
	obj.multiple_run_plot()
	obj.show_polar_graphs()
	
