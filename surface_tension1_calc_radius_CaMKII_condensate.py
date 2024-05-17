
##
## Plot profile - CaMKII and GluN2B only
## Need d['rdf']['CaMKII_bead']
## 


import os, sys, glob, pickle, pprint, copy

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1

import lib.utils as utils
import lib.parameters as p
import lib.colormap as c

from scipy import ndimage

plt.rcParams.update(p.rc_param)
plt.rcParams.update( {'font.size': 6} )
	
	
def peak(x, c):
    return np.exp(-np.power(x - c, 2) / 16.0)
	
def lin_interp(x, y, i, half):
    return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))
	
def half_max_x(x, y):
    half = max(y)/2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    return lin_interp(x, y, zero_crossings_i[-1], half)
	
	
def plot_a_rdf( ax, r, y, legend=False,  ylim = (-0.006,0.66) ):

	color = c.cmap_universal_ratio['CaMKII']
	ax.step( r, y, color=color, label='CaMKII bead')
	if legend==True:
		ax.legend(frameon=False)
	
	ax.set_xlabel('Distance from \n center-of-mass (l.u.)')
	ax.set_ylabel('(beads / voxel)')
	ax.set_xlim(0,40)
	ax.set_ylim(*ylim)
	ax.set_xticks(range(0,50,10))
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	
	return
	
	
def arrange_graph_bar(ax, panel_size_x, panel_size_y):
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax_h, ax_w = ax.bbox.height, ax.bbox.width
	loc = ax.get_position()
	y0 = loc.y0
	ax.set_position([loc.x0, y0, panel_size_x, panel_size_y])
	
	
def plot_bar(ax, data, color, ylim, ylegend, width=0.5):
	ax.bar(*zip(*data.items()), width=width, color=color)
	ax.set_ylabel(ylegend)
	ax.set_ylim(ylim)
	

class HandleRDFCaMKII():
		
	def __init__( self ):
		dir_target = 'CG_valency_length'
		subdir     = 'radius_CaMKII'
		self.basename = 'radius_CaMKII'
		
		self.dir_lammpstrj    = os.path.join('..', 'lammpstrj4', dir_target)
		self.dir_edited_data  = os.path.join('data4', dir_target)
		self.dir_imgs = os.path.join('imgs4', dir_target, subdir)
		
		os.makedirs(self.dir_edited_data, exist_ok=True)
		os.makedirs(self.dir_imgs, exist_ok=True)
		
		self.valencies = range(2,14,2)
		self.lengths   = range(7)
		
		self.num_rows	 = len( self.valencies )
		self.num_columns = len( self.lengths   ) 
		
		
	def run( self ):
		self.fig, self.axes  =  plt.subplots(self.num_rows, self.num_columns, figsize=(10, 10), tight_layout=True)
		self.hmws = {}
		for i, v in enumerate(self.valencies):
			for j, l in enumerate(self.lengths):
				# Filename
				filename_input = os.path.join('val{}'.format(v), 'R2_{}.lammpstrj'.format(str(l).zfill(3)) )
				prefix = '{}_{}'.format(str(v).zfill(2), str(l).zfill(3))
				suffix = 'sigma_2'
				# d      = utils.load(dir_edited_data, prefix, suffix)
				
				# Load lammpstrj
				print("\n"+filename_input)
				sampling_frame = utils.get_num_frames(self.dir_lammpstrj, filename_input)
				rdf, rdf_bins  = \
					utils.get_rdf_CaMKII(self.dir_lammpstrj, filename_input, sampling_frame )
				
				# Make figure
				row    = self.num_rows-i-1
				column = j
				self.hmws[prefix] = self.plot_a_graph(row, column, rdf, rdf_bins, prefix)
		
		
	def plot_a_graph(self, row, column, rdf, rdf_bins, title):
		#print('rdf')
		#print(rdf)
		# Plot RDF
		rdf_bins = rdf_bins[4:-1]
		rdf      = rdf['CaMKII_bead'][4:]
		r = rdf_bins
		y = rdf
		
		hmw = half_max_x(r, y)
		print('hmw ', hmw)
		plot_a_rdf( self.axes[row][column], r, y ) 
		self.axes[row][column].set_title( '{}, r = {:.3f}'.format(title, hmw) )
		
		return hmw
		
	def save_figs( self ):
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.svg' ) )
		self.fig.savefig( os.path.join(self.dir_imgs, self.basename + '.png' ), dpi=150 )
		plt.show()
		plt.clf()
		plt.close(fig=self.fig)
		
	def save_hmws( self ):
		prefix = 'radiuses'
		suffix = 'CaMKII'
		print('{}_{}'.format(prefix, suffix))
		utils.save(self.dir_edited_data, prefix, suffix, self.hmws)		
		
		
if __name__ == '__main__':
	
	obj = HandleRDFCaMKII()
	obj.run()
	obj.save_figs()
	obj.save_hmws()
	
	
