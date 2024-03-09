
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

def plot_a_binding(ax, d, prefix, species, vv ):
	
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
		ymax = 3000
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
		ymax = 1200
	elif species == 'CaMKII':
		titles_binding_type = ['{} GluN2B'.format(i) for i in range(vv+1)]
		binding_types = [ [0]*(vv - i + 1) +['GluN2B_CaMKII'] * i for i in range(vv+1)]
		ymax = 6000 / (vv-1)
		if vv == 2:
			binding_types = [ [0, 0], [0, 'GluN2B_CaMKII'], ['GluN2B_CaMKII','GluN2B_CaMKII'] ]
			ymax = 5000


	numbers_binding_type = np.zeros(len(binding_types), dtype='int')
	for id, v in  d[species].items():
		binding = v['binding_types']
		for j, binding_ref in enumerate( binding_types ):
			#if vv == 2:
			#	print('Compare ', binding, binding_ref)
			if equal_list(binding, binding_ref):
				numbers_binding_type[j] += 1
	ax.grid()
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
	species = 'CaMKII' # 'STG','GluN2B', 'PSD95', 'CaMKII'
	

	dir_target  = 'valency_linker_length'
	dir_edited_data	= os.path.join('..','data2', dir_target)
	dir_imgs = os.path.join('imgs2', dir_target, 'binding')
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	
	valency = list(range(2,14,2)) 
	linker_length  = [1, 3, 5, 7, 9]
	
	fnames_valency       = { v: str(v).zfill(2) for v in valency }
	fnames_linker_length = {ll: str(i).zfill(3) for i,ll in enumerate(linker_length) }
	
	num_rows		= len( valency )
	num_columns		= len( linker_length )
	
	
	vals = np.zeros([num_rows, num_columns])
	fig  = plt.figure(figsize=(10, 10), tight_layout=True)
	#
	
	for i, v in enumerate(valency):
		for j, ll in enumerate(linker_length):
			# Load data
			prefix = fnames_valency[v]+'_'+fnames_linker_length[ll]
			suffix = '_'
			d      = utils.load(dir_edited_data, prefix, suffix)
			print('Target: {}'.format(prefix))
			# sys.exit(0)
			title = prefix # prefix, None
			row    = num_rows-i-1
			column = j+1
			print('column: ', column, ', row: ', row)
			ax = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
			print('valency  ', v)
			plot_a_binding(ax, d, prefix, species, v)
	
	# Save figure
	filename = 'binding_partner_{}'.format(species)
	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename ) ) , dpi=150)
	plt.show()
	plt.clf()
	plt.close(fig=fig)
	
	
	