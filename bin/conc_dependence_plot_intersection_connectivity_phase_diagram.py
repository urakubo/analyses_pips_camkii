
import os, glob, pickle, pprint, copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from skimage.measure import label

import lib.utils as utils
import lib.utils_panel as utils_panel
import lib.parameters as p
import lib.colormap as c
from specification_datasets import SpecDatasets

plt.rcParams.update( p.rc_param )

		
class PhaseDiagramConcDependence():
	def __init__( self ):
		
		GluN2Bs = np.array([270,540,1080,2160,4320,6480,8640,12960,17280])
		GluN2Bs = GluN2Bs / np.prod(p.space_np)
		GluN2Bs = GluN2Bs * 1000
		GluN2Bs = np.flip(GluN2Bs)
		#GluN2Bs = np.log(GluN2Bs)
		self.GluN2Bs = GluN2Bs
		
		
		STGs = np.array([108,216,432,576,864,1728,2592,3456,4320,5184])
		STGs = STGs / np.prod(p.space_np)
		STGs = STGs * 1000
		#STGs = np.log(STGs)
		self.STGs = STGs
		
		
		self.num_rows = len( GluN2Bs )
		self.title    = 'Phase diagram'
		self.basename = 'phase_diagram_conc_dependence'
		# -1: Unclear
		# 1: Homogeneous LLPS (CaMKII)
		# 2: PIPS
		# 3: Homogeneous LLPS (STG)
		
		phase_diagram = [\
			[ 1, 1, 1, 1,   1, 1, 1, 1,   1, 1],
			[ 1, 1, 1, 1,   1, 1, 1, 1,   1, 1],
			[ 1, 1, 1, 1,   1, 1, 1, 1,   1, 1],
			[ 1, 1, 1, 1,   1, 1, 1, 1,   1, 1],
			[ 1, 1, 1, 1,   1, 1, 1, 1,   1, 1],
			[ 3, 3, 3, 3,   3, 3, 3, 3,   3, 3],
			[ 3, 3, 3, 3,   3, 3, 3, 3,   3, 3],
			[ 3, 3, 3, 3,   3, 3, 3, 3,   3, 3],
			[ 3, 3, 3, 3,   3, 3, 3, 3,   3, 3]]
			
		# 1: Partial engulfment
		phase_diagram_partial_engulfment = [\
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0],
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0],
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0],
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0],
			[ 0, 0, 0, 0,   0, 0, 1, 1,   1, 1],
			[ 0, 0, 0, 0,   0, 0, 1, 1,   1, 1],
			[ 0, 0, 0, 0,   0, 0, 0, 1,   1, 1],
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0],
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0]]
		phase_diagram_two_condensates = [\
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0],
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0],
			[ 0, 0, 0, 0,   1, 1, 1, 1,   1, 1],
			[ 0, 0, 1, 1,   1, 1, 1, 1,   1, 1],
			[ 0, 0, 1, 1,   1, 1, 1, 1,   1, 1],
			[ 0, 0, 0, 0,   0, 1, 1, 1,   1, 1],
			[ 0, 0, 0, 0,   0, 0, 0, 1,   1, 1],
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0],
			[ 0, 0, 0, 0,   0, 0, 0, 0,   0, 0]]
			
			
		self.phase_diagram = np.array(phase_diagram).T
		self.phase_diagram_partial_engulfment = np.array(phase_diagram_partial_engulfment).T
		self.phase_diagram_two_condensates = np.array(phase_diagram_two_condensates).T
		


def set_colors( data, target_type ):
	
	## Set Color
	m = 255
	gray = [230/m,230/m,230/m]
	cols = []
	for two_phase, partial_engulfment, phase_diagram in \
		zip( data['two_phase_condensate'], data['partial_engulfment'], data['phase_diagram'] ):
		if target_type == 'phase_diagram':
			if (two_phase == 1) and (partial_engulfment == 0):
				cols.append([255/m,191/m/2,192/m/2]) # 'r'
			elif (two_phase == 1) and (partial_engulfment == 1):
				cols.append('w')
			elif (two_phase == 0) and (phase_diagram == 3):
				cols.append([77/m,196/m, 255/m]) # 'b'
			else:
				cols.append(gray)
		elif target_type == 'two_phase_condensate':
			cols.append('k' if two_phase == 1 else 'w')
		elif target_type == 'pips_partial':
			cols.append('k' if partial_engulfment == 0 else 'w')
		else:
			sys.exit('Unknown target type: {}'.format(target_type) )
			
	return cols
	

class PhaseDiagramConcDependenceIntersection(SpecDatasets, PhaseDiagramConcDependence):
	
	def __init__( self ):
		
		PhaseDiagramConcDependence.__init__(self)
		
		
	def prepare_plot( self ):
		#self.fig  = plt.figure(figsize=(5, 5))
		fig  = plt.figure(figsize=(4, 4))
		fig.subplots_adjust(wspace=0.4,  hspace=0.6)
		ax = fig.add_subplot( 1, 1, 1 )
		ax.set_xlabel('STG (beads / voxel) x 10-3')	
		ax.set_ylabel('GluN2B (beads / voxel) x 10-3')
		return fig, ax
		
		
		
	def load_data( self ):
		
		species1 = 'CaMKII'
		type_analysis1 = 'average'
		self.level1  = 0.852
		prefix1 = '{}_{}'.format( species1, type_analysis1 )
		
		species2 = 'PSD95'
		type_analysis2 = 'average_STG'
		self.level2    = 0.407
		prefix2   = '{}_{}'.format( species2, type_analysis2 )
		
		suffix = 'data'
		self.data1 = utils.load( self.dir_edited_data, prefix1, suffix  )
		self.data2 = utils.load( self.dir_edited_data, prefix2, suffix  )
		
		
	def plot_data( self ):
		mx_min = 0
		my_min = 0
		mx_max = np.max(self.concs_stg[:-1]) * 1.05# None # 2.6
		my_max = np.max(self.concs_glun2b)   * 1.05 # None # 10.1
		
		# Plot panel
		self.fig, self.ax = self.prepare_plot()
		utils_panel.plot_a_panel_intersection(self.ax, \
			self.data1, self.level1, \
			self.data2, self.level2, \
			self.concs_stg, self.concs_glun2b, \
			mx_min=mx_min, my_min=my_min, \
			mx_max=mx_max, my_max=my_max, \
			margin = 0
			)
		
		# Plot points
		
		
		
	def plot_overlay_points( self ):
		variable_names = [	'two_phase_condensate',
							'partial_engulfment',
							'phase_diagram' ]
							
		objects = [	self.phase_diagram_two_condensates,
					self.phase_diagram_partial_engulfment,
					self.phase_diagram]
		self.data = {vname: np.fliplr(obj).T[:,:-1].flatten() for obj, vname in zip( objects, variable_names) }
		
		
		
		target_type = 'phase_diagram' # 'phase_diagram', 'two_phase_condensate', 'pips_partial'
		cols = set_colors( self.data, target_type )
		
		xx, yy = np.meshgrid(self.concs_stg[:-1], self.concs_glun2b)
		xx = xx.flatten()
		yy = yy.flatten()
		for x, y, c in zip(xx, yy, cols):
			self.ax.plot(x, y,'o', \
				markersize = 4, \
				color = 'k', \
				markerfacecolor = c )
		
	def save_plots( self ):
		
		self.basename = 'CaMKII_STG_intersection'
		dir_imgs = os.path.join(self.dir_imgs_root, 'phase_diagram')
		os.makedirs(dir_imgs, exist_ok=True)
		self.fig.savefig( os.path.join(dir_imgs, self.basename + '.svg' ) )
		self.fig.savefig( os.path.join(dir_imgs, self.basename + '.png' ) , dpi=150)
		plt.show()
		plt.clf()
		plt.close(fig=self.fig)
		
		
	def _modify_data(self, d):
		#print('d ')
		#print(d[self.species][self.type_analysis])
		if self.species == 'CaMKII' and self.type_analysis == 'average':
			data      = d[self.species][self.type_analysis]['GluN2B']
		elif self.species == 'PSD95' and self.type_analysis == 'average':
			data      = d[self.species][self.type_analysis]['STG_PSD95']
		elif self.species == 'PSD95' and self.type_analysis == 'average_GluN2B':
			data      = d[self.species]['average']['GluN2B_PSD95']
		elif self.species == 'PSD95' and 'ratio' in self.type_analysis:
			both   = d[self.species][self.type_analysis]['Both']
			STG    = d[self.species][self.type_analysis]['STG only']
			GluN2B = d[self.species][self.type_analysis]['GluN2B only']
			No_bound= d[self.species][self.type_analysis]['None']
			num_total = both + STG + GluN2B + No_bound
			data      = both / num_total
		return data
		
		
		
if __name__ == '__main__':
	pass
