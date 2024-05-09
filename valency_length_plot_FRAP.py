
import os, sys, glob, pickle, pprint, copy, pickle
import numpy as np

from scipy.optimize import curve_fit

import matplotlib
import matplotlib.pyplot as plt

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c
import lib.utils_graph as utils_graph

plt.rcParams.update(p.rc_param)


class MatrixCGValencyLength():
	def __init__( self ):
		plt.rcParams.update( {'font.size': 6} )
		# Input files
		self.valencies = range(4,14,2)
		self.lengths   = p.lengths2
		
		self.filenames_edited    = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in self.valencies for id_f in self.lengths]
		dir_target  = 'small_colony2'
		self.suffix = 'FRAP'
		
		self.dir_edited_data = os.path.join('data4', dir_target)
		self.dir_imgs        = os.path.join('imgs4', dir_target, 'matrix')
		os.makedirs(self.dir_imgs, exist_ok=True)
		
		self.num_rows		= len( self.valencies )
		self.num_columns	= len( self.lengths   )
		
	def run( self ):
		
		vals = {}
		self.fig  = plt.figure(figsize=(8, 6), tight_layout=True)
		#fig.subplots_adjust(wspace=0.4,  hspace=0.6)
		
		for i, v in enumerate( self.valencies ):
			for j, l in enumerate( self.lengths ):
				# Load data
				prefix = p.fnames_valency2[v]+'_'+p.fnames_length2[l]
				row    = self.num_rows-i-1
				column = j+1
				print('Target file: ', prefix, ', column: ', column, ', row: ', row)
				d      = utils.load(self.dir_edited_data, prefix, self.suffix)
				
				legend = (v == 12) and (l == 1)
				
				title = prefix # prefix, None
				vv, _ = self.plot_a_graph(row, column, d, title, legend)
				vals[prefix] = np.log10(vv)
		return vals
		
		
	def save( self ):
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.svg' ) )
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.png' ), dpi=150 )
		plt.show()
		plt.clf()
		plt.close(fig=self.fig)
	
	
	
class SingleValencyLength():
	def __init__( self ):
		
		# Input files
		self.valency = 12
		self.length   = 2
		
		self.filename_edited = str(self.valency).zfill(2)+'_'+str(self.length).zfill(3)
		self.basename = self.filename_edited
		dir_target  = 'small_colony2'
		self.suffix = 'FRAP'
		
		self.dir_edited_data = os.path.join('data4', dir_target)
		self.dir_imgs        = os.path.join('imgs4', dir_target, 'FRAP')
		os.makedirs(self.dir_imgs, exist_ok=True)
		self.num_rows = 1
		self.num_columns = 1
		
	def plot( self ):
		d = utils.load(self.dir_edited_data, self.filename_edited, self.suffix)
		
		self.fig  = plt.figure(figsize=(1.3, 1.3))
		
		legend = True
		row = 0
		column = 1
		value, ax = self.plot_a_graph(row, column, d, self.basename, legend)
		ax.set_xlim([-1,8])
		ax.set_xticks([0,2,4,6,8])
		ax.set_title('{},{:.3f}'.format( self.basename, value))
		
	def save( self ):
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.svg' ) )
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.png' ), dpi=150 )
		plt.show()
		plt.clf()
		plt.close(fig=self.fig)
	
	
class PlotFRAP():
	def __init__( self ):
		
		self.suffix   = 'FRAP'
		self.basename = 'FRAP'
		
	def plot_a_graph(self, row, column, d, title, legend = False):
		
		time_frame = d['time_steps'] /1e9
		molecular_concentration_in_target_area = d['molecular_concentration_in_target_area']
		
		ax = self.fig.add_subplot( self.num_rows, self.num_columns, row*self.num_columns+column )
		ax.set_title(title)
		
		ax.set_xlabel('Time (109 MC steps)')
		ax.set_ylabel('CaMKII (% baseline)')
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.set_xlim([np.min(time_frame), np.max(time_frame)])
		
		for i, mname in enumerate(['GluN2B', 'CaMKII']):
			#ax.set_title(mname)
			tmp = [ molecular_concentration_in_target_area[:,i] for i in p.molecules_with_all[mname]['id'] ]
			density = sum(tmp)
			density = density / np.mean(density[time_frame < 0]) * 100
			ax.plot(\
				time_frame, density,\
				color = p.molecules_with_all[mname]['c'], label=mname)
			
			# Fitting!
			if mname == 'CaMKII':
				density_for_fit    = density[time_frame > 0]
				time_frame_for_fit = time_frame[time_frame > 0]
				
				def func_exponential_recovery(x, tau, a, b):
					return a*(-np.exp(-x/tau)/b + 1)
					
				min_a  , max_a   = 60, 100
				min_b  , max_b   = 1, 10
				if title in ['12_000','12_001','10_000','08_000']:
					max_tau,min_tau = 400,0.01
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
				ax.plot(time_frame_for_fit, func_exponential_recovery(time_frame_for_fit, pp[0],pp[1],pp[2]) ,\
					color = 'k')
				
			#def func_exponential(x, tau, a, b):
			#    return a*(np.exp(-x/tau) + b)
			
		if legend == True:
			ax.legend(frameon=False, loc = 'lower right')
		
		return pp[0], ax
	
class PlotFRAPMatrixCGValencyLength(PlotFRAP, MatrixCGValencyLength):
	def __init__( self ):
		PlotFRAP.__init__(self )
		MatrixCGValencyLength.__init__(self)


class PlotFRAPSingleValencyLength(PlotFRAP, SingleValencyLength):
	def __init__( self ):
		PlotFRAP.__init__(self )
		SingleValencyLength.__init__(self)


	
if __name__ == '__main__':
	
	##
	## Define Input
	##
	
	'''
	CG_valnecy_length = PlotFRAPMatrixCGValencyLength()
	values = CG_valnecy_length.run()
	CG_valnecy_length.save()
	
	dir_edited_data = CG_valnecy_length.dir_edited_data
	prefix = 'FRAP'
	suffix = 'matrix'
	utils.save(dir_edited_data, prefix, suffix, values)
	'''
	
	valnecy_length = PlotFRAPSingleValencyLength()
	values = valnecy_length.plot()
	valnecy_length.save()
	

