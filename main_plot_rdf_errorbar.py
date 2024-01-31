import os, glob

import numpy as np
import matplotlib.pyplot as plt
import utils
import colormap as c

plt.rcParams.update({
                    'pdf.fonttype' : 'truetype',
                    'svg.fonttype' : 'none',
                    'font.family' : 'sans-serif',
                    'font.sans-serif' : 'Arial',
                    'font.style' : 'normal'})
	


if __name__ == '__main__':
	
	# Parameters
	i = 18
	filename_input = 'R3_{}.lammpstrj'.format( str(i).zfill(3) )
	dir_data        = './../231018Graphs/data'
	dir_imgs 		= 'imgs/rdf_errorbar'
	reference_molecule_for_centering = 'All'
	
	os.makedirs(dir_imgs, exist_ok=True)
	prefix = os.path.splitext(filename_input)[0]
	
	
	# For rdf
	rdf_bins        = np.arange(0, 35) # You can set "np.arange(0, 35, 2)"
	rdf_grid_points = utils.get_lattice_grids()
	
	
	# Set target sampling frames
	num_target_frames = 10
	sampling_interval = 2
	num_frames        = utils.get_num_frames(dir_data, filename_input)
	rdf_target_frames = list(range( num_frames - num_target_frames*sampling_interval, \
							num_frames,\
							sampling_interval))
	print("Target frames ", rdf_target_frames)
	
	
	rdfs = utils.get_rdf_from_multiple_frame(dir_data, filename_input, \
			rdf_target_frames, rdf_bins, rdf_grid_points, reference_molecule_for_centering)
	
	
	# Plot the histogram
	fig = plt.figure( figsize=(5, 4), tight_layout=True )
	ax = fig.add_subplot( 1, 1, 1 )
	
	for k, rdf in rdfs.items():
		rdf_mean = np.mean(rdf, axis=1)
		rdf_std  = np.std(rdf, axis=1)
		color    = c.cmap_universal_ratio[k]
		
		ax.hist(rdf_bins[:-1], rdf_bins, weights=rdf_mean, \
			color=color, \
			histtype="step", label=k)
		ax.errorbar(rdf_bins[:-1]+0.5,rdf_mean, yerr=rdf_std,\
			ecolor=color,\
			linestyle='',
			alpha = 0.4 )
	
	ax.set_title( prefix )
	ax.legend(frameon=False)
	
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_xlabel('Distance from center-of-mass (l.u.)')
	ax.set_ylabel('(beads / volume)')
	
	plt.savefig( os.path.join(dir_imgs, 'rdf_{}.svg'.format(prefix) ) )
	plt.savefig( os.path.join(dir_imgs, 'rdf_{}.png'.format(prefix) ) , dpi=150)
	plt.show()
	
	
