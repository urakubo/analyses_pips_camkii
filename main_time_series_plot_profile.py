
import os, glob, pickle, pprint, copy

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1

import utils
import colormap as c
from scipy import ndimage


plt.rcParams.update({
                    'pdf.fonttype' : 'truetype',
                    'svg.fonttype' : 'none',
                    'font.family' : 'sans-serif',
                    'font.sans-serif' : 'Arial',
                    'font.style' : 'normal'})
	
	
def make_a_figure( d, target_molecules ):
	
	# Parameters
	vmax       = 0.7

	# Figure specification
	transps = [(0,1,2),(1,2,0),(2,0,1)]
	titles  = [ 'yz', 'xy', 'zx']
	num_columns = 10
	num_rows    = 3
	fig = plt.figure(figsize=(20, 8)) # , tight_layout=False
	
	left, right, bottom, top = 0.0, 0.95, 0.10, 0.99
	wspace, hspace = 0.2, 0.1
	plt.subplots_adjust(left, bottom, right, top, wspace, hspace)
	fig.suptitle( "Time step: {}, MC step: {:e}".format( d['time_frame'], d['mc_step'] ) )
	
	panel_size = 0.15
	
	# Plot profiles
	yloc = []
	for row,(title, transp) in enumerate(zip(titles, transps)):
		axes = []
		columns = target_molecules
		axes1 = utils.plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp, vmax = vmax )
		axes.extend(axes1)
		loc0  = axes1[0].get_position()
		yloc.append(loc0.y0)
		for a in axes:
			loc = a.get_position()
			#print('loc ', loc)
			a.set_position([loc.x0, yloc[-1], panel_size, panel_size])
	
	
	panel_dx = 0.03
	
	# Plot RDF
	column = 1
	ax = fig.add_subplot( num_rows, num_columns, num_columns*2 )
	errorbar= 'shaded' # errorbar='shaded', 'line', or 'final_frame_alone'
	utils.plot_a_rdf( ax, d, errorbar=errorbar, legend=True, target_molecules = target_molecules.keys(), ylim = (-vmax/20,vmax) )
	utils.arrange_graph_bar(ax, panel_dx, yloc[2], panel_size/2, panel_size)
	
	return fig
	
	
if __name__ == '__main__':
	
	# Input files 1
	dir_edited_data 		= 'data'
	filenames_edited_data 	= [	'CG',\
							'SP',\
							'SPG1',\
							'SPG2']
	targets_molecules 		= [ {'CaMKII':2, 'GluN2B':3} ,\
							{'STG':4,'PSD95':5},\
							{'GluN2B':3, 'STG':4,'PSD95':5},\
							{'GluN2B':3, 'STG':4,'PSD95':5}]
	sampling_interval = 20
	sampling_time_frames_ = list(range(0,200+sampling_interval,sampling_interval))
	sampling_time_framess = 4*[sampling_time_frames_]


	# Input files 2
	dir_edited_data 		= 'data'
	filenames_edited_data = ['CaMKIIalone',\
							'GluN2Balone',\
							'STGalone']
	targets_molecules 		= [ {'CaMKII':2, 'GluN2B':3,'STG':4,'PSD95':5} ] * 3
	sampling_time_framess = [ [0, 100, 200], [0, 100, 187] , [0, 100, 195] ]


	# Output files
	dir_imgs = os.path.join('imgs', 'time_series','profiles')
	os.makedirs(dir_imgs, exist_ok=True)
	
	# Parameters
	
	
	for filename, target_molecules, sampling_time_frames in zip(filenames_edited_data, targets_molecules, sampling_time_framess):
		print('Recorded time frame: ', sampling_time_frames)
		for target_frame in sampling_time_frames:
			# Load data
			prefix = filename
			suffix = 'frame_{}'.format(target_frame)
			print('Target: {}, {}'.format(filename, suffix))
			d   = utils.load(dir_edited_data, prefix, suffix)
			
			# Process data
			# labels_watershed, volume_ratios_watershed = watershed_segmentation( d )
			fig = make_a_figure(d, target_molecules)
			
			# Save figure
			fig.savefig( os.path.join(dir_imgs, '{}_frame_{}.svg'.format( prefix, str(target_frame).zfill(4)) ) )
			fig.savefig( os.path.join(dir_imgs, '{}_frame_{}.png'.format( prefix, str(target_frame).zfill(4)) ) , dpi=150)
			#plt.show()
			plt.clf()
			plt.close(fig=fig)
