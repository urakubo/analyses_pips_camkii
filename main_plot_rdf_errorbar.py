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
	dir_data	= 'data'
	dir_imgs 	= 'imgs/rdf_errorbar'
	i = 18
	filename_prefix = str(i).zfill(3)
	sigma = 2 # Same data are stored in 2,3,4
	filename_suffix = 'sigma_{}'.format(sigma)
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	# Load file
	d   = utils.load(dir_data, filename_prefix, filename_suffix)
	rdf_bins 			= d['rdf_bins']
	rdf_target_frames 	= d['rdf_target_frames']
	rdfs 				= d['rdf']
	
	
	print("Target time frames ", rdf_target_frames)
	
	
	
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
	
	ax.set_title( filename_prefix )
	ax.legend(frameon=False)
	
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_xlabel('Distance from center-of-mass (l.u.)')
	ax.set_ylabel('(beads / volume)')
	
	plt.savefig( os.path.join(dir_imgs, 'rdf_{}.svg'.format(filename_prefix) ) )
	plt.savefig( os.path.join(dir_imgs, 'rdf_{}.png'.format(filename_prefix) ) , dpi=150)
	plt.show()
	
	
