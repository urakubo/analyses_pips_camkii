
import os, sys, glob, pickle, pprint
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import utils
import colormap as c


	
	
if __name__ == '__main__':
	# Parameters
	
	# Dataset 1
	dir_input       = './../lammpstrj/example_profiles'
	filenames_input = [		'length3.lammpstrj',\
							'length9.lammpstrj',\
							'Partial_engulfment.lammpstrj',\
							'Homogeneous.lammpstrj']
	
	dir_output       = 'data'
	filenames_output = [	'PIPS',\
							'iPIPS',\
							'PartialE',\
							'Homo']
	
	i = 0
	reference_molecule_for_centering = 'All'
	filename_input, filename_output = filenames_input[i], filenames_output[i]
	
	
	# Load data
	print("\n"+filename_input)
	num_frames     = utils.get_num_frames(dir_input, filename_input)
	print("num_frames ", num_frames )
	types, positions,ids_molecule = utils.load_data( dir_input, filename_input, num_frames )
	print("The last timeframe was loaded." )
	# Centering
	center    = utils.get_center(types, positions, reference_molecule_for_centering)
	positions = utils.centering(positions, center)	
	#
	unique_ids_molecule = np.unique( ids_molecule )
	ids_CaMKII = [id for id in unique_ids_molecule if utils.subunits['CaMKII Hub']['id'] in types[ids_molecule == id] ]
	ids_GluN2B = [id for id in unique_ids_molecule if utils.subunits['GluN2Bc']['id']    in types[ids_molecule == id] ]
	ids_STG    = [id for id in unique_ids_molecule if utils.subunits['STGc1']['id']    in types[ids_molecule == id] ]
	ids_PSD    = [id for id in unique_ids_molecule if utils.subunits['PSD1']['id']    in types[ids_molecule == id] ]
	ids = {'CaMKII':ids_CaMKII, 'GluN2B':ids_GluN2B , 'STG':ids_STG, 'PSD95': ids_PSD }
	
	
	# Plot
	m = 25
	rng = np.random.default_rng(seed=13)
	
	fig = plt.figure(figsize=(10, 7), tight_layout=True)
	fig.subplots_adjust(wspace=0.6, hspace=0.6)
	for i, t in enumerate([ 'CaMKII','PSD95']): # 'CaMKII','PSD95','GluN2B','STG'
		ax = fig.add_subplot(1,2,i+1, projection='3d')
		ax.set_title(t)
		for id in ids[t]:
			pos = positions[ids_molecule == id]
			dist = np.max( np.linalg.norm(pos, axis=1) )
			if dist < m:
				if rng.random() > 0.8:
					ax.plot(pos[:,0], pos[:,1], pos[:,2],\
						'-' ,\
						color = c.cmap_universal_ratio[t], \
						linewidth = 0.5)
	

		ax.set_xlim(-m,m)
		ax.set_ylim(-m,m)
		ax.set_zlim(-m,m)
		ax.set_box_aspect([1,1,1])
		ax.set_xlabel('(l.u.)')
		ax.set_ylabel('(l.u.)')
		ax.set_zlabel('(l.u.)')	

	
	plt.savefig('my-plot.png' , dpi=150)
	plt.show()
