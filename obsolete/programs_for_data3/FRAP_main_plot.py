
import os, sys, glob, pickle, pprint, copy, pickle
import numpy as np


import matplotlib
import matplotlib.pyplot as plt

import utils
import parameters as p
import colormap as c


plt.rcParams.update(p.rc_param)
plt.rcParams.update( {'font.size': 6} )

class MatrixCGValencyLength():
	def __init__( self ):
		
		# Input files
		dir_target      = 'CG_valency_length'
		self.dir_edited_data = os.path.join('data3', dir_target)
		self.dir_imgs        = os.path.join('imgs3', dir_target, 'matrix')
		os.makedirs(self.dir_imgs, exist_ok=True)
		
		self.num_rows		= len( p.valencies )
		self.num_columns	= len( p.lengths   )
		
	def run( self ):
		
		vals = np.zeros([self.num_rows, self.num_columns])
		self.fig  = plt.figure(figsize=(8, 8), tight_layout=True)
		#fig.subplots_adjust(wspace=0.4,  hspace=0.6)
		
		for i, v in enumerate( p.valencies ):
			for j, ll in enumerate( p.lengths ):
				# Load data
				prefix = p.fnames_valency[v]+'_'+ p.fnames_length[ll]
				row    = self.num_rows-i-1
				column = j+1
				print('Target file: ', prefix, ', column: ', column, ', row: ', row)
				d      = utils.load(self.dir_edited_data, prefix, self.suffix)
				
				legend = (v == 12) and (ll == 1)
				
				title = prefix # prefix, None
				vals[row, column-1] = self.plot_a_graph(row, column, d, title, legend)
		
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
		
		time_frame = d['time_steps']
		molecular_concentration_in_target_area = d['molecular_concentration_in_target_area']
		types_     = d['types']
		positions_ = d['positions']
		
		ax = self.fig.add_subplot( self.num_rows, self.num_columns, row*self.num_columns+column )
		ax.set_title(title)
		
		ax.set_xlabel('Time (MC moves)')
		ax.set_ylabel('Num / voxel (% baseline)')
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.set_xlim([np.min(time_frame), np.max(time_frame)])
		
		for i, mname in enumerate(['GluN2B', 'CaMKII']):
			#ax.set_title(mname)
			tmp = [ molecular_concentration_in_target_area[:,i] for i in p.molecules_with_all[mname]['id'] ]
			density = sum(tmp)
			baseline = np.mean(density[time_frame < 0])
			ax.plot(\
				time_frame, density / baseline * 100,\
				color = p.molecules_with_all[mname]['c'], label=mname)
		if legend == True:
			ax.legend(frameon=False, loc = 'lower right')
		
		return True
	
class PlotFRAPMatrixCGValencyLength(PlotFRAP, MatrixCGValencyLength):
	def __init__( self ):
		PlotFRAP.__init__(self )
		MatrixCGValencyLength.__init__(self)
	
if __name__ == '__main__':
	
	##
	## Define Input
	##
	
	target = 'conc_GluN2B'
	plot_CG_valnecy_length = PlotFRAPMatrixCGValencyLength() # PlotConcMatrixConcDependence
	values = plot_CG_valnecy_length.run()
	plot_CG_valnecy_length.save()
	
	
	# CG_ Valency length
	'''
	filenames_input = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
	filenames_input = ['12_000']
	dir_target    = 'CG_valency_length'
	pprint.pprint(filenames_input)
	
	#print('Input data')
	#pprint.pprint(filenames_input)
	
	# Shared init
	dir_edited_data	= os.path.join('data3', dir_target)
	dir_imgs = os.path.join('imgs3', dir_target,'FRAP')
	os.makedirs(dir_imgs, exist_ok=True)
	
	##
	##
	##
	##
	filename_input  = filenames_input[-1] # -34
	
	# Load data
	prefix = filename_input
	suffix = 'FRAP'
	print('Prefix: {}, suffix: {}'.format(prefix, suffix))
	d   = utils.load(dir_edited_data, prefix, suffix)
	
	time_frame	= d['modified_time_frame']
	molecular_concentration_in_target_area = d['molecular_concentration_in_target_area']
	types_ = d['types']
	positions_ = d['positions']
	
	
	print("Plot the FRAP")
	fig = plt.figure(figsize=(8,8))
	#fig.suptitle(filename)
	#fig.subplots_adjust(hspace=0.6, wspace=0.4)
	#gs = fig.add_gridspec(3,1)
	
	
	def init_panel(ax):
		ax.set_xlabel('Time (frame)')
		ax.set_ylabel('Num per grid volume (% baseline)')
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.set_xlim([np.min(time_frame), np.max(time_frame)])
		#ax.set_xlim([np.min(time_frame), 60])

	ax = fig.add_subplot(2,2,1)
	for i, mname in enumerate(['CaMKII','GluN2B']):
		init_panel(ax)
		#ax.set_title(mname)
		density  = molecular_concentration_in_target_area[:,p.molecules_with_all[mname]['id']]
		baseline = np.mean(density[time_frame < 0])
		ax.plot(\
			time_frame, density / baseline * 100,\
			color = p.molecules_with_all[mname]['c'], label=mname)

	ax.legend(frameon=False, loc = 'lower right')
	
	fig.savefig(os.path.join(dir_imgs,prefix+'_FRAP.pdf'))
	fig.savefig(os.path.join(dir_imgs,prefix+'_FRAP.png'),dpi=150)
	plt.show()
	plt.close()
	'''

