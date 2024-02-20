
import os, glob, pickle, pprint
import numpy as np
import utils
import parameters as p
import matplotlib.pyplot as plt
plt.rcParams.update(p.rc_param)



if __name__ == '__main__':
	
	target = 'rdf' # 'region_condensates', 'conc_CaMKII', 'conc_STG', 'rdf'
	
	
	# Files
	dir_edited_data 		= os.path.join('data', 'conc_dependence')
	filenames_edited_data 	= [str(i).zfill(3) for i in range(3) ] # 30
	target_dir = 'valency_linker_length'
	
	
	dir_edited_data	= os.path.join('data', target_dir)
	dir_imgs = os.path.join('imgs', target_dir,'matrix')
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	valency = [4, 6, 8, 10, 12] 
	linker_length  = [3, 5, 7, 9]
	
	num_rows		= len( valency )
	num_columns		= len( linker_length )
	
	
	# Other params
	errorbar= 'shaded' # errorbar='shaded', 'line', or 'final_frame_alone'
	sigma  = 2
	
	fig  = plt.figure(figsize=(8, 8), tight_layout=True)
	#fig.subplots_adjust(wspace=0.4,  hspace=0.6)
	
	for i, v in enumerate(valency):
		for j, ll in enumerate(linker_length):
			# Load data
			id = i + j * len(valency)
			prefix = str(id).zfill(3)
			suffix = 'sigma_{}'.format(sigma)
			d      = utils.load(dir_edited_data, prefix, suffix)
			print('Target: {}, sigma: {}'.format(prefix, sigma))
			
			
			# Specify row and column
			row    = num_rows-i-1
			column = j+1
			print('column: ', column, ', row: ', row)
			scalebar = True if (row == 0) and (column == 1) else False
			
			
			if target == 'region_condensates':
				utils.plot_regions_condenstate_from_a_direction(fig, num_rows, num_columns, row, column, d, title=False, scalebar=scalebar )
			elif target == 'conc_CaMKII':
				columns = {'CaMKII':column}
				utils.plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, \
					title=None, colorbar=False, scalebar=scalebar, pre_rotated=True )
			elif target == 'conc_STG':
				columns = {'STG':column}
				utils.plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, \
					title=None, colorbar=False, scalebar=scalebar, pre_rotated=True )
			elif target == 'rdf':
				plt.rcParams.update( {'font.size': 6} )
				ax    = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
				utils.plot_a_rdf( ax, d, errorbar=errorbar, legend=scalebar )
			else:
				raise ValueError("Target function has not been implemented: {}".format(target))
			#  transps = [(0,1,2),(1,2,0),(2,0,1)]
	
	# Save figure
	fig.savefig( os.path.join(dir_imgs, 'matrix_{}_sigma_{}.svg'.format( target, sigma ) ) )
	fig.savefig( os.path.join(dir_imgs, 'matrix_{}_sigma_{}.png'.format( target, sigma ) ) , dpi=150)
	#plt.show()
	plt.clf()
	plt.close(fig=fig)


