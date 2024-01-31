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
		ax.hist(rdf_bins[:-1], rdf_bins, weights=v, color=c.cmap_universal_ratio[k], histtype="step", label=k)
	
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_xlabel('Distance from center-of-mass (l.u.)')
	ax.set_ylabel('(beads per volume)')
	
	
if __name__ == '__main__':
	
	# Parameters
	filenames_input = ['R3_{}.lammpstrj'.format( str(i).zfill(3) ) for i in range(30)]
	dir_data        = './../231018Graphs/data'
	dir_imgs 		= 'imgs/rdf_summary'
	reference_molecule_for_centering = 'All'
	os.makedirs(dir_imgs, exist_ok=True)
	
	# For rdf
	rdf_bins        = list(range(0,30))
	rdf_grid_points = utils.get_lattice_grids()
	
	
	rdfs = {}
	for i, filename_input in enumerate( filenames_input ):
		# Set the last frame as a target frame
		target_frame = utils.get_num_frames(dir_data, filename_input)
		types, positions,_ = utils.load_data( dir_data, filename_input, target_frame )
		# Get a rdf
		rdfs[i] = utils.get_rdf(types, positions, \
			rdf_grid_points, rdf_bins, \
			reference_molecule_for_centering)
	
	print('Plot')
	fig = plt.figure(figsize=(15, 15), tight_layout=True)
	for i, filename_input in enumerate( filenames_input ):
		# print('i: {}, dists {}'.format(i, dists))
		ax = fig.add_subplot( 6, 5, i+1 )
		plot_hist(ax, rdfs[i], rdf_bins)
		ax.set_title( os.path.splitext(filename_input)[0] )
		if i == 0:
			ax.legend(frameon=False)
	
	fig.savefig( os.path.join(dir_imgs, 'rdf_summary.svg') )
	fig.savefig( os.path.join(dir_imgs, 'rdf_summary.png') , dpi=150)
	plt.show()
	
	
