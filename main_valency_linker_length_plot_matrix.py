
import os, glob, pickle, pprint
import numpy as np
import utils
import parameters as p
import matplotlib.pyplot as plt
plt.rcParams.update(p.rc_param)
plt.rcParams.update( {'font.size': 6} )


if __name__ == '__main__':
	
	target =  'energy_isotropic_STG'
	# 'region_condensates', 'conc_CaMKII', 'conc_STG', 'rdf', 'concs_in_CaMKII', 'concs_in_STG'
	# 'energy_anisotropic_STG', 'energy_anisotropic_CaMKII'
	# 'energy_isotropic_STG', 'energy_isotropic_CaMKII', 
	
	# Files
	dir_target  = 'valency_linker_length'
	dir_edited_data	= os.path.join('data2', dir_target)
	dir_imgs = os.path.join('imgs2', dir_target, 'matrix')
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	valency = list(range(2,14,2)) 
	linker_length  = [1, 3, 5, 7, 9]
	
	fnames_valency       = { v: str(v).zfill(2) for v in valency }
	fnames_linker_length = {ll: str(i).zfill(3) for i,ll in enumerate(linker_length) }
	
	num_rows		= len( valency )
	num_columns		= len( linker_length )
	
	
	# Other params
	sigma  = 2
	
	
	vals = np.zeros([num_rows, num_columns])
	fig  = plt.figure(figsize=(8, 8), tight_layout=True)
	#fig.subplots_adjust(wspace=0.4,  hspace=0.6)
	
	for i, v in enumerate(valency):
		for j, ll in enumerate(linker_length):
			# Load data
			prefix = fnames_valency[v]+'_'+fnames_linker_length[ll]
			suffix = 'sigma_{}'.format(sigma)
			d      = utils.load(dir_edited_data, prefix, suffix)
			print('Target: {}, sigma: {}'.format(prefix, sigma))
			
			# Specify row and column
			row    = num_rows-i-1
			column = j+1
			print('column: ', column, ', row: ', row)
			
			title = prefix # prefix, None
			vals[row, column-1] = utils.select_plot(target, fig, num_rows, num_columns, row, column, d, title)
			
			
	# Save figure
	fig.savefig( os.path.join(dir_imgs, 'matrix_{}_sigma_{}.svg'.format( target, sigma ) ) )
	fig.savefig( os.path.join(dir_imgs, 'matrix_{}_sigma_{}.png'.format( target, sigma ) ) , dpi=150)
	#plt.show()
	plt.clf()
	plt.close(fig=fig)
	print(target)
	print(vals)


