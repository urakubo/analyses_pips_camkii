
import os, glob, pickle, pprint
import numpy as np
import utils
import parameters as p
import matplotlib.pyplot as plt
plt.rcParams.update(p.rc_param)



if __name__ == '__main__':
	
	target   =  'region_condensates' # 'region_condensates', 'conc_CaMKII', 'conc_STG', 'watershed_CaMKII', 'rdf'
	sigma    = 2 # 2 or 4
	
	# Input files
	dir_edited_data 		= os.path.join('data2', 'conc_dependence')
	filenames_edited_data 	= [str(i).zfill(3) for i in range(3) ] # 30
	# Output files
	dir_imgs = os.path.join('imgs2', 'conc_dependence','matrix')
	os.makedirs(dir_imgs, exist_ok=True)
	
	#STG    = [1, 500, 1000, 1500, 2000, 2500, 3000] 
	#GluN2B = [1, 500, 1000, 2000, 4000, 6000, 8000, 12000, 16000, 20000]
	
	STG    = [540 , 1620, 2160,  2700,  3240, 4520] 
	GluN2B = [1080, 4320, 8640, 10800, 12960]
	
	volume = np.prod(p.space_np)
	STG    = [ s / volume for s in STG    ]
	GluN2B = [ n / volume for n in GluN2B ]
	
	num_rows		= len( GluN2B )
	num_columns		= len( STG )
	
	#volume_watershed_CaMKII = np.zeros([len(STG), len(GluN2B)], dtype=float)
	#volume_watershed_STG    = np.zeros([len(STG), len(GluN2B)], dtype=float)
	
	
	fig  = plt.figure(figsize=(10, 10), tight_layout=True)
	#fig.subplots_adjust(wspace=0.4,  hspace=0.6)
	
	for i, stg in enumerate(STG):
		for j, glun in enumerate(GluN2B):
			# Load data
			id = i + j * len(STG)

			prefix = str(id).zfill(3)
			suffix = 'sigma_{}'.format(sigma)
			d      = utils.load(dir_edited_data, prefix, suffix)
			print('Target: {}, sigma: {}'.format(prefix, sigma))
			
			if id == (len(GluN2B)-1)*len(STG):
				scalebar=True
			else:
				scalebar=False
			
			title = prefix # prefix, None
			errorbar= 'shaded' # errorbar='shaded', 'line', or 'final_frame_alone'
			
			row    = num_rows-j-1
			column = i+1
			print('column: ', column, ', row: ', row)
			if target == 'region_condensates':
				utils.plot_regions_condenstate_from_a_direction(fig, num_rows, num_columns, row, column, d, title=title, scalebar=scalebar )
			elif target == 'conc_CaMKII':
				columns = {'CaMKII':column}
				utils.plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, \
					title=title, colorbar=False, scalebar=scalebar, pre_rotated=True )
			elif target == 'conc_STG':
				columns = {'STG':column}
				utils.plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, \
					title=title, colorbar=False, scalebar=scalebar, pre_rotated=True )
			elif target == 'watershed_CaMKII':
				columns = {'CaMKII':column}
				utils.plot_watershed_region_from_a_direction(fig, num_rows, num_columns, row, columns, d, title=False, scalebar=scalebar )
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


