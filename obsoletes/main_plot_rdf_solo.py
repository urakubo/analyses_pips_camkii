import os, glob

import numpy as np
import matplotlib.pyplot as plt
import utils


plt.rcParams.update({
                    'pdf.fonttype' : 'truetype',
                    'svg.fonttype' : 'none',
                    'font.family' : 'sans-serif',
                    'font.sans-serif' : 'Arial',
                    'font.style' : 'normal'})
	


if __name__ == '__main__':
	
	# Parameters
	dir_data	= 'data'
	dir_imgs 	= 'imgs/rdf_solo'
	i = 18
	filename_prefix = str(i).zfill(3)
	sigma = 2 # Same data are stored in 2,3,4
	filename_suffix = 'sigma_{}'.format(sigma)
	os.makedirs(dir_imgs, exist_ok=True)
	
	errorbar= 'shaded' # errorbar='shaded', 'line', or 'final_frame_alone'
	
	# Load file
	d   = utils.load(dir_data, filename_prefix, filename_suffix)
	
	
	# Plot the histogram
	fig = plt.figure( figsize=(5, 4), tight_layout=True )
	ax = fig.add_subplot( 1, 1, 1 )
	utils.plot_a_rdf( ax, d, errorbar=errorbar, legend=True  )
	
	ax.set_title( filename_prefix )
	
	fig.savefig( os.path.join(dir_imgs, 'rdf_{}_{}.svg'.format(filename_prefix, errorbar) ) )
	fig.savefig( os.path.join(dir_imgs, 'rdf_{}_{}.png'.format(filename_prefix, errorbar) ) , dpi=150)
	plt.show()
	
	
