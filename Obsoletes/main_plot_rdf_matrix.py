import os, glob

import numpy as np
import matplotlib.pyplot as plt

from scipy import ndimage
import pickle
import utils


plt.rcParams.update({
                    'pdf.fonttype' : 'truetype',
                    'svg.fonttype' : 'none',
                    'font.family' : 'sans-serif',
                    'font.sans-serif' : 'Arial',
                    'font.style' : 'normal'})
	
	
if __name__ == '__main__':
	
	# Parameters
	filenames_prefix = [str(i).zfill(3) for i in range(30)]
	sigma = 2 # Same data are stored in 2,3,4
	filename_suffix = 'sigma_{}'.format(sigma)
	dir_data	    = 'data'
	dir_imgs 		= 'imgs/rdf_matrix'
	
	os.makedirs(dir_imgs, exist_ok=True)
	
	errorbar= 'final_frame_alone' # errorbar='shaded', 'line', or 'final_frame_alone'
	
	print('Plot')
	fig = plt.figure(figsize=(15, 15), tight_layout=True)
	for i, filename_prefix in enumerate( filenames_prefix ):
		
		# Load file
		d   = utils.load(dir_data, filename_prefix, filename_suffix)
		
		ax = fig.add_subplot( 6, 5, i+1 )
		legend = (i == 0)
		utils.plot_a_rdf( ax, d, errorbar=errorbar, legend=legend )
		ax.set_title( filename_prefix )

	
	fig.savefig( os.path.join(dir_imgs, 'rdf_matrix_{}.svg'.format(errorbar)) )
	fig.savefig( os.path.join(dir_imgs, 'rdf_matrix_{}.png'.format(errorbar)) , dpi=150)
	plt.show()
	
	
