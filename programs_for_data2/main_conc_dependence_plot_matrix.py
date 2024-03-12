
import os, glob, pickle, pprint
import numpy as np
import utils
import parameters as p
import matplotlib.pyplot as plt
plt.rcParams.update(p.rc_param)
plt.rcParams.update( {'font.size': 6} )


if __name__ == '__main__':
	
	target = 'concs_in_STG'
	# 'region_condensates', 'conc_CaMKII', 'conc_STG', 'rdf', 'concs_in_CaMKII', 'concs_in_STG'
	# 'energy_anisotropic_STG', 'energy_anisotropic_CaMKII'
	# 'energy_isotropic_STG', 'energy_isotropic_CaMKII', 
	# 'energy_anisotropic_dilute', 'energy_isotropic_dilute'
	

	sigma    = 2 # 2 or 4
	
	# Input files
	dir_edited_data 		= os.path.join('data2', 'conc_dependence')
	filenames_edited_data 	= [str(i).zfill(3) for i in range(3) ] # 30
	# Output files
	dir_imgs = os.path.join('imgs2', 'conc_dependence','matrix')
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	STG    = [540 , 1620, 2160,  2700,  3240, 4520] 
	GluN2B = [1080, 4320, 8640, 10800, 12960]
	
	volume = np.prod(p.space_np)
	STG    = [ s / volume for s in STG    ]
	GluN2B = [ n / volume for n in GluN2B ]
	
	num_rows		= len( GluN2B )
	num_columns		= len( STG )
	
	
	vals = np.zeros([num_rows, num_columns])
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
			
			title = prefix # prefix, None
			row    = num_rows-j-1
			column = i+1
			print('column: ', column, ', row: ', row)
			vals[row, column-1] = utils.select_plot(target, fig, num_rows, num_columns, row, column, d, title)
			
			
	# Save figure
	fig.savefig( os.path.join(dir_imgs, 'matrix_{}_sigma_{}.svg'.format( target, sigma ) ) )
	fig.savefig( os.path.join(dir_imgs, 'matrix_{}_sigma_{}.png'.format( target, sigma ) ) , dpi=150)
	#plt.show()
	plt.clf()
	plt.close(fig=fig)
	print(target)
	print(vals)

