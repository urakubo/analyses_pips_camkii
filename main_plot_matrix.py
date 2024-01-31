
import os, glob, pickle, pprint
import numpy as np
import utils
import matplotlib.pyplot as plt
plt.rcParams.update({
                    'pdf.fonttype' : 'truetype',
                    'svg.fonttype' : 'none',
                    'font.family' : 'sans-serif',
                    'font.sans-serif' : 'Arial',
                    'font.style' : 'normal'})



if __name__ == '__main__':

	dir_data = 'data'
	dir_imgs = 'imgs/matrix'
	target   = 'watershed_CaMKII' # 'region_condensates', 'profile_CaMKII', 'profile_STG', 'watershed_CaMKII'
	sigma    = 2 # 2, 3, or 4
	
	os.makedirs(dir_imgs, exist_ok=True)
	
	STG    = [500, 1000, 2000, 3000, 4000]
	GluN2B = [500, 2000, 4000, 6000, 8000, 12000]
	
	volume_watershed_CaMKII = np.zeros([len(STG), len(GluN2B)], dtype=float)
	volume_watershed_STG    = np.zeros([len(STG), len(GluN2B)], dtype=float)
	
	num_rows		= len( GluN2B )
	num_columns		= len( STG )
	
	
	fig  = plt.figure(figsize=(10, 10), tight_layout=True)
	#fig.subplots_adjust(wspace=0.4,  hspace=0.6)
	
	for i, stg in enumerate(STG):
		for j, glun in enumerate(GluN2B):
			# Load data
			id = i + j * len(STG)
			prefix = str(id).zfill(3)
			suffix = 'sigma_{}'.format(sigma)
			d      = utils.load(dir_data, prefix, suffix)
			print('Target: {}, sigma: {}'.format(prefix, sigma))
			
			row    = num_rows-j-1
			column = i+1
			print('column: ', column, ', row: ', row)
			if target == 'region_condensates':
				utils.plot_regions_condenstate_from_a_direction(fig, num_rows, num_columns, row, column, d, title=False )
			elif target == 'profile_CaMKII':
				columns = {'CaMKII':column}
				utils.plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, title=False, colorbar=False )
			elif target == 'profile_STG':
				columns = {'STG':column}
				utils.plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, title=False, colorbar=False )
			elif target == 'watershed_CaMKII':
				columns = {'CaMKII':column}
				utils.plot_watershed_region_from_a_direction(fig, num_rows, num_columns, row, columns, d, title=False )
			else:
				raise ValueError("Target function has not been implemented: {}".format(target))
			#  transps = [(0,1,2),(1,2,0),(2,0,1)]
	
	# Save figure
	fig.savefig( os.path.join(dir_imgs, 'matrix_{}_sigma_{}.svg'.format( target, sigma ) ) )
	fig.savefig( os.path.join(dir_imgs, 'matrix_{}_sigma_{}.png'.format( target, sigma ) ) , dpi=150)
	#plt.show()
	plt.clf()
	plt.close(fig=fig)


