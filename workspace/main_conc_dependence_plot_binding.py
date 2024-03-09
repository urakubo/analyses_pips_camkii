
import os, sys, glob, pickle, pprint
import numpy as np
import math

import matplotlib
import matplotlib.pyplot as plt

sys.path.append('../')

import utils
import parameters as p
import colormap as c

plt.rcParams.update(p.rc_param)
plt.rcParams.update( {'font.size': 6} )

def equal_list(lst1, lst2):
    lst = lst1.copy()
    for element in lst2:
        try:
            lst.remove(element)
        except ValueError:
            break
    else:
        if not lst:
            return True
    return False

def plot_a_binding(ax, d, prefix, species = 'GluN2B'):
	
	if species == 'GluN2B':
		binding_types = [ [0, 0], \
			[0, 'GluN2B_CaMKII'], ['GluN2B_CaMKII', 'GluN2B_CaMKII'], \
			[0, 'GluN2B_PSD95'], ['GluN2B_PSD95', 'GluN2B_PSD95'], ['GluN2B_CaMKII', 'GluN2B_PSD95']]
		titles_binding_type = ['None', 'CaMKII, None', 'CaMKII, CaMKII', 'PSD95, None',  'PSD95, PSD95', 'CaMKII, PSD95' ]
		ymax = 6000
		
	elif species == 'STG':
		binding_types = [\
			[0, 0, 0, 0], \
			[0, 0, 0, 'STG_PSD95'], \
			[0, 0, 'STG_PSD95', 'STG_PSD95'], \
			[0, 'STG_PSD95', 'STG_PSD95', 'STG_PSD95'], \
			['STG_PSD95', 'STG_PSD95', 'STG_PSD95', 'STG_PSD95'], \
			]
		titles_binding_type = ['None', '1 PSD95', '2 PSD95', '3 PSD95', '4 PSD95' ]
		ymax = 6000
	elif species == 'PSD95':
		binding_types = [\
			[0, 0, 0], \
			[0, 0, 'STG_PSD95'],[0, 'STG_PSD95', 'STG_PSD95'],['STG_PSD95', 'STG_PSD95', 'STG_PSD95'], \
			[0, 0, 'GluN2B_PSD95'], [0, 'GluN2B_PSD95', 'GluN2B_PSD95'],['GluN2B_PSD95', 'GluN2B_PSD95', 'GluN2B_PSD95'],\
			[0, 'STG_PSD95', 'GluN2B_PSD95'],\
			['STG_PSD95', 'STG_PSD95', 'GluN2B_PSD95'],\
			['STG_PSD95', 'GluN2B_PSD95', 'GluN2B_PSD95']\
			]
		titles_binding_type = ['None', \
				'1 STG', '2 STG', '3 STG', \
				'1 GluN2B', '2 GluN2B', '3 GluN2B', \
				'1 STG, 1 GluN2B', '2 STG, 1 GluN2B', '1 STG, 2 GluN2B'
				 ]
		ymax = 3000

	numbers_binding_type = np.zeros(len(binding_types), dtype='int')
	for id, v in  d[species].items():
		binding = v['binding_types']
		for j, binding_ref in enumerate( binding_types ):
			if equal_list(binding, binding_ref):
				numbers_binding_type[j] += 1
	
	ax.set_title(prefix)
	ax.bar(titles_binding_type, numbers_binding_type, width=0.5)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_ylim(0,ymax)
	ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
	ax.set_ylabel('(Number)')
	# ax.set_aspect(1.0 / ax.get_data_ratio())


if __name__ == '__main__':

	
	# Dataset 1
	species = 'PSD95' # 'STG','GluN2B', 'PSD95'
	
	dir_target = 'conc_dependence'
	dir_edited_data  = 'conc_dependence_cluster'
	dir_edited_data  = os.path.join('..','data2',dir_edited_data)
	
	dir_imgs = os.path.join('imgs2', dir_target)
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
	#
	for i, stg in enumerate(STG):
		for j, glun in enumerate(GluN2B):
			# Load data
			id = i + j * len(STG)
			prefix = str(id).zfill(3)
			suffix = '_'
			d      = utils.load(dir_edited_data, prefix, suffix)
			print('Target: {}'.format(prefix))
			
			# sys.exit(0)
			title = prefix # prefix, None
			row    = num_rows-j-1
			column = i+1
			print('column: ', column, ', row: ', row)
			ax = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
			plot_a_binding(ax, d, prefix, species)
	
	# Save figure
	filename = 'binding_partner_{}'.format(species)
	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename ) ) , dpi=150)
	plt.show()
	plt.clf()
	plt.close(fig=fig)
	
	
	