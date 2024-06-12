
import os, sys, glob, pickle, pprint, copy, pickle
import numpy as np

from scipy.optimize import curve_fit

import matplotlib
import matplotlib.pyplot as plt

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph


from specification_datasets import SpecDatasets
from lib.utils_matrix import MatrixValencyLength

plt.rcParams.update(p.rc_param)


class SingleValencyLength():
	def __init__( self ):
		
		self.num_rows = 1
		self.num_columns = 1
		
		
	def plot( self , valency = 12, length = 2 ):
		
		
		if self.target_molecule_fit == 'CaMKII' and length == 2:
			#self.dir_target  = 'small_colony2'
			xlim        = [-1,8]
			set_xticks  = [0,2,4,6,8]
		elif self.target_molecule_fit == 'GluN2B' and length == 2:
			#self.dir_target  = 'small_colony3'
			xlim        = [-0.01,0.08]
			set_xticks  = [0,0.04,0.08]
		elif self.target_molecule_fit in ['CaMKII','GluN2B'] and length == 6:
			#self.dir_target  = 'small_colony3'
			xlim        = [-0.01,0.08]
			set_xticks  = [0,0.04,0.08]
		else:
			print('Error!')
			sys.exit()
		
		
		# Load data
		prefix =  self.filename_edited_matrix(valency, length)
		d = utils.load(self.dir_edited_data, prefix, self.suffix)
		
		# Plot
		self.fig  = plt.figure(figsize=(1.3, 1.0))
		legend = True
		row = 0
		column = 1
		value, ax = self.plot_a_graph(row, column, d, self.basename, legend)
		ax.set_xlim( xlim )
		ax.set_xticks(set_xticks)
		ax.set_title('{},{:.4f}'.format( self.basename, value))
		#
		
		
	def save( self ):
	
		dir_imgs = os.path.join(self.dir_imgs_root, 'FRAP')
		os.makedirs(dir_imgs, exist_ok=True)
		self.fig.savefig( os.path.join(dir_imgs, self.basename + '.svg' ) )
		self.fig.savefig( os.path.join(dir_imgs, self.basename + '.png' ), dpi=150 )
		plt.show()
		plt.clf()
		plt.close(fig=self.fig)
	
	
	
def func_exponential_recovery(x, tau, a, b):
	return a*(-np.exp(-x/tau)/b + 1)
	
class PlotFRAP():
	
	def exponential_fitting( self, title, density, time_frame ):
		density_for_fit    = density[time_frame > 0]
		time_frame_for_fit = time_frame[time_frame > 0]
		
		min_a  , max_a   = 60, 100
		min_b  , max_b   = 1, 10
		if title in ['12_000','12_001','10_000','08_000']:
			max_tau, min_tau = 400, 0.01
			p0 = [100, 70, 2]
		elif ('006' in title) or ('005' in title)or ('004' in title):
			max_tau,min_tau = 0.5,0.01
			min_b  , max_b   = 1, 5
			p0 = [0.1, 70, 2]
		else:
			max_tau,min_tau = 10,0.01
			p0 = [5, 70, 2]

		pp, cov = curve_fit(func_exponential_recovery, time_frame_for_fit, density_for_fit,\
			p0=p0,\
			bounds = ((min_tau, min_a, min_b), (max_tau, max_a, max_b)),\
			maxfev=5000)
		return pp

	def __init__( self, target_molecule_fit ):
		
		self.suffix   = 'FRAP'
		self.basename = 'FRAP_' + target_molecule_fit
		self.target_molecule_fit = target_molecule_fit
		
		
	def plot_a_graph(self, row, column, d, title, legend = False):
		
		time_frame = d['time_steps'] /1e9
		molecular_concentration_in_target_area = d['molecular_concentration_in_target_area']
		
		ax = self.fig.add_subplot( self.num_rows, self.num_columns, row*self.num_columns+column )
		ax.set_title(title)
		
		ax.set_xlabel('Time (10^9 MC steps)')
		ax.set_ylabel('CaMKII (% baseline)')
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		
		for i, mname in enumerate(['GluN2B', 'CaMKII']):
			tmp = [ molecular_concentration_in_target_area[:,i] for i in p.molecules_with_all[mname]['id'] ]
			density = sum(tmp)
			
			#print(mname)
			#print('density[time_frame < 0] ', density[time_frame < 0])
			
			density = density / np.mean(density[time_frame < 0]) * 100
			
			if self.target_molecule_fit == 'GluN2B':
				density_  = density[time_frame < 0.08]
				time_frame_ = time_frame[time_frame < 0.08]
			else:
				density_  = density
				time_frame_ = time_frame
			
			ax.set_xlim([np.min(time_frame_), np.max(time_frame_)])
			ax.plot(\
				time_frame_, density_,\
				color = p.molecules_with_all[mname]['c'], label=mname)
			
			# Fitting
			if mname == self.target_molecule_fit:
				param = self.exponential_fitting( title, density_, time_frame_ )
				time_frame_for_fit = time_frame_[time_frame_ > 0]
				ax.plot(time_frame_for_fit, func_exponential_recovery( time_frame_for_fit, *param.tolist() ) ,\
					color = 'k')
		
		if legend == True:
			ax.legend(frameon=False, loc = 'lower right')
		
		return param[0], ax
		
		
class PlotFRAPMatrixValencyLength(PlotFRAP, MatrixValencyLength, SpecDatasets):
	def __init__( self, target_molecule = 'CaMKII' ):
		
		PlotFRAP.__init__( self, target_molecule )
		MatrixValencyLength.__init__( self )
		SpecDatasets.__init__( self )
		
		
		
class PlotFRAPSingleValencyLength(PlotFRAP, SingleValencyLength, SpecDatasets):
	def __init__( self, target_molecule = 'CaMKII' ):
		PlotFRAP.__init__(self, target_molecule )
		SingleValencyLength.__init__(self)


	
if __name__ == '__main__':
	
	
	# Calculate taus and save them.
	
	'''
	valnecy_length = PlotFRAPMatrixValencyLength( target_molecule = 'GluN2B')
	values = valnecy_length.run()
	valnecy_length.save()
	
	dir_edited_data = valnecy_length.dir_edited_data
	prefix = 'FRAP_GluN2B'
	suffix = 'matrix'
	utils.save(dir_edited_data, prefix, suffix, values)
	'''
	
	'''
	dir_target  = 'small_colony2' # 'small_colony2', 'small_colony3'
	valnecy_length = PlotFRAPMatrixValencyLength(dir_target  = dir_target)
	fitting_tau    = valnecy_length.run()
	valnecy_length.save()
	
	fitting_tau = {k: np.log10(v) for k, v in fitting_tau.items()}
	
	dir_edited_data = os.path.join('data4', dir_target)
	prefix          = 'FRAP'
	suffix          = 'matrix'
	utils.save(dir_edited_data, prefix, suffix, fitting_tau)
	'''
	
	
	# Single profiles of FRAP
	target_molecule = 'CaMKII' # 'GluN2B', 'CaMKII'
	length          = 6 # 2, 6
	
	valnecy_length = PlotFRAPSingleValencyLength(target_molecule = 'CaMKII',length = 6)
	values = valnecy_length.plot()
	valnecy_length.save()
	
	
	# Two FRAP files are merged and saved.
	'''
	prefix = 'FRAP'
	suffix = 'matrix'
	dir_edited_data = os.path.join('data4', 'small_colony2')
	d_colony2 = utils.load(dir_edited_data, prefix, suffix)
	dir_edited_data = os.path.join('data4', 'small_colony3')
	d_colony3 = utils.load(dir_edited_data, prefix, suffix)
	
	for valency_length in ['12_000','10_000','08_000','12_001','10_001','12_002']:
		d_colony3[valency_length] = d_colony2[valency_length]
	
	dir_edited_data = os.path.join('data4', 'small_colony2')
	prefix = 'FRAP_merged'
	utils.save(dir_edited_data, prefix, suffix, d_colony3)
	'''
	
