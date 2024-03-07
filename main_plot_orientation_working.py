
import os, sys, glob, pickle, pprint
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import utils
import parameters as p
import colormap as c

	
	
if __name__ == '__main__':
	
	i = 2
	
	# Dataset 1
	dir_input       = './../lammpstrj2/Feb_Figure2'
	fnames = [27,21,15,9]
	filenames_input  = [ 'R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in fnames ]
	dir_target       = 'conc_dependence'
	
	filenames_output = [ format(str(i).zfill(3)) for i in fnames ]
	dir_imgs = os.path.join('imgs2', dir_target,'prof3d')
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	reference_molecule_for_centering = 'All'
	filename_input, filename_output = filenames_input[i], filenames_output[i]
	
	
	# Load data
	print("\n"+filename_input)
	num_frames     = utils.get_num_frames(dir_input, filename_input)
	print("num_frames ", num_frames )
	types, positions_grid_coord,ids_molecule, time_stamp = utils.load_lammpstrj( dir_input, filename_input, num_frames )
	print("The last timeframe was loaded." )
	# Centering
	center    = utils.get_center_of_mass(types, positions_grid_coord)
	positions_real_coord = utils.centering(positions_grid_coord, center)
	#
	unique_ids_molecule = np.unique( ids_molecule )
	ids_CaMKII = [id for id in unique_ids_molecule if p.subunits['CaMKII hub']['id'] in types[ids_molecule == id] ]
	ids_GluN2B = [id for id in unique_ids_molecule if p.subunits['GluN2B binding site']['id']    in types[ids_molecule == id] ]
	ids_STG    = [id for id in unique_ids_molecule if p.subunits['STG binding site']['id']    in types[ids_molecule == id] ]
	ids_PSD    = [id for id in unique_ids_molecule if p.subunits['PSD']['id']    in types[ids_molecule == id] ]
	ids = {'CaMKII':ids_CaMKII, 'GluN2B':ids_GluN2B , 'STG':ids_STG, 'PSD95': ids_PSD }
	
	
	# Plot
	m = 25
	rng = np.random.default_rng(seed=13)
	
	fig = plt.figure(figsize=(10, 7), tight_layout=True)
	fig.subplots_adjust(wspace=0.6, hspace=0.6)
	for i, t in enumerate(['CaMKII','GluN2B','PSD95']): # 'CaMKII','PSD95','GluN2B','STG'
		ax = fig.add_subplot(2,2,i+1, projection='3d')
		ax.set_title(t)
		for id in ids[t]:
			pos = positions_real_coord[ids_molecule == id]
			dist = np.max( np.linalg.norm(pos, axis=1) )
			if dist < m:
				if rng.random() > 0.8:
					if t == 'CaMKII':
						for i in range(12):
							ax.plot([pos[0,0],pos[i+1,0]], [pos[0,1], pos[i+1,1]], [pos[0,2],pos[i+1,2]],\
								'-' ,\
								color = c.cmap_universal_ratio[t], \
								linewidth = 0.5)
					elif t in ['PSD95','GluN2B']:
						ax.plot(pos[:,0], pos[:,1], pos[:,2], '-' ,\
							color = c.cmap_universal_ratio[t], \
							linewidth = 0.5)
	

		ax.set_xlim(-m,m)
		ax.set_ylim(-m,m)
		ax.set_zlim(-m,m)
		ax.set_box_aspect([1,1,1])
		ax.set_xlabel('(l.u.)')
		ax.set_ylabel('(l.u.)')
		ax.set_zlabel('(l.u.)')	

	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename_output ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename_output ) ) , dpi=150)
	plt.show()
