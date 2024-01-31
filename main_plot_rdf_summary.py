import os, glob

import numpy as np
import matplotlib.pyplot as plt

from scipy import ndimage
import pickle
import utils
import colormap as c

plt.rcParams.update({
                    'pdf.fonttype' : 'truetype',
                    'svg.fonttype' : 'none',
                    'font.family' : 'sans-serif',
                    'font.sans-serif' : 'Arial',
                    'font.style' : 'normal'})
	

def plot_hist(ax, rdf, rdf_bins):
	for k, v in rdf.items():
		ax.hist(rdf_bins[:-1], rdf_bins, weights=v[:,-1], color=c.cmap_universal_ratio[k], histtype="step", label=k)
	
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_xlabel('Distance from center-of-mass (l.u.)')
	ax.set_ylabel('(beads per volume)')
	
	
if __name__ == '__main__':
	
	# Parameters
	filenames_prefix = [str(i).zfill(3) for i in range(30)]
	sigma = 2 # Same data are stored in 2,3,4
	filename_suffix = 'sigma_{}'.format(sigma)
	dir_data	    = 'data'
	dir_imgs 		= 'imgs/rdf_summary'
	
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	print('Plot')
	fig = plt.figure(figsize=(15, 15), tight_layout=True)
	for i, filename_prefix in enumerate( filenames_prefix ):
		
		# Load file
		d   = utils.load(dir_data, filename_prefix, filename_suffix)
		rdf_bins 			= d['rdf_bins']
		rdf_target_frames 	= d['rdf_target_frames']
		rdfs 				= d['rdf']
		
		# print('i: {}, dists {}'.format(i, dists))
		ax = fig.add_subplot( 6, 5, i+1 )
		plot_hist(ax, rdfs, rdf_bins)
		ax.set_title( filename_prefix )
		if i == 0:
			ax.legend(frameon=False)
	
	fig.savefig( os.path.join(dir_imgs, 'rdf_summary.svg') )
	fig.savefig( os.path.join(dir_imgs, 'rdf_summary.png') , dpi=150)
	plt.show()
	
	
