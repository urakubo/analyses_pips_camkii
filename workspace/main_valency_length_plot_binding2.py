
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


def plot_a_binding(ax, d, prefix, species):
	

	if species == 'GluN2B':
		# np.unique(d['GluN2B']['nums']['GluN2B_CaMKII'], return_counts=True)
		ave_CaMKII = np.average(d[species]['nums']['GluN2B_CaMKII'])
		ave_PSD95 = np.average(d[species]['nums']['GluN2B_PSD95'])
		ave_numbers_binding = [
			ave_CaMKII,
			ave_PSD95,
			]
		titles_binding_type = ['CaMKII \n {:.2f}'.format(ave_CaMKII), 'PSD95 \n {:.2f}'.format(ave_PSD95)]
		ymax = 2
	elif species == 'PSD95':
		ave_GluN2B = np.average(d[species]['nums']['GluN2B_PSD95'])
		ave_STG = np.average(d[species]['nums']['STG_PSD95'])
		ave_numbers_binding = [
			ave_GluN2B,
			ave_STG
			]
		titles_binding_type = ['GluN2B \n {:.2f}'.format(ave_GluN2B), 'STG \n {:.2f}'.format(ave_STG)]
		ymax = 3
	elif species == 'CaMKII':
		ave_GluN2B = np.average(d[species]['nums']['GluN2B_CaMKII'])
		ave_numbers_binding = [ ave_GluN2B ]
		titles_binding_type = ['GluN2B \n {:.2f}'.format(ave_GluN2B)]
		ymax = 10
	elif species == 'STG':
		ave_PSD95 = np.average(d[species]['nums']['STG_PSD95'])
		ave_numbers_binding = [ ave_PSD95 ]
		titles_binding_type = ['PSD95 \n {:.2f}'.format(ave_PSD95)]
		ymax = 4
		
		
	ax.set_title(prefix)
	ax.bar(titles_binding_type, ave_numbers_binding, width=0.5)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_ylim(0,ymax)
	# ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
	ax.set_ylabel('(Number)')
	# ax.set_aspect(1.0 / ax.get_data_ratio())



if __name__ == '__main__':

	
	# Dataset 1
	species = 'STG' # 'STG','GluN2B', 'PSD95', 'CaMKII'
	

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
			plot_a_binding(ax, d, prefix, species)
	
	# Save figure
	filename = 'ave_num_binding_{}'.format(species)
	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename ) ) , dpi=150)
	plt.show()
	plt.clf()
	plt.close(fig=fig)
	
	
	