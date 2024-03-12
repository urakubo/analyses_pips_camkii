
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
		ymax = 12
	elif species == 'STG':
		ave_PSD95 = np.average(d[species]['nums']['STG_PSD95'])
		ave_numbers_binding = [ ave_PSD95 ]
		titles_binding_type = ['PSD95 \n {:.2f}'.format(ave_PSD95)]
		ymax = 4
	elif species == 'PSD95_connection':
		ave_numbers_binding = list(d[species].values())
		titles_binding_type = list(d[species].keys()) #['None', 'STG only', 'GluN2B only','Both']
		ymax = 3000
		
	ax.set_title(prefix)
	ax.bar(titles_binding_type, ave_numbers_binding, width=0.5)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_ylim(0,ymax)
	# ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
	ax.set_ylabel('(Number)')
	# ax.set_aspect(1.0 / ax.get_data_ratio())

	if species == 'PSD95_connection':
		ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

if __name__ == '__main__':

	
	# Dataset 1
	species = 'PSD95_connection' # 'STG','GluN2B', 'PSD95','CaMKII', 'PSD95_connection'
	
	dir_target = 'conc_dependence'
	dir_edited_data  = 'conc_dependence_cluster'
	dir_edited_data  = os.path.join('..','data2',dir_edited_data)
	
	dir_imgs = os.path.join('imgs2', dir_target)
	os.makedirs(dir_imgs, exist_ok=True)
	suffix = 'connection'
	
	
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
	filename = 'ave_num_binding_{}'.format(species)
	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename ) ) , dpi=150)
	plt.show()
	plt.clf()
	plt.close(fig=fig)
	
	
	