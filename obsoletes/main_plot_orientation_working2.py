
import os, sys, glob, pickle, pprint
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


from scipy.spatial import distance

import utils
import parameters as p
import colormap as c
	
if __name__ == '__main__':
	linear =False
	
	#'''
	target = 'CGSP' #  ['CGSP', 'CG','CPG','PG','SP'] # ,'SPG'
	filename_input 	= '{}.lammpstrj'.format(target)
	dir_target      = 'mixtures'
	filename_edited = target
	filename_output = target
	#'''
	
	#'''
	subdirs    = ['CG_con', 'CG_len9', 'CG_lin']
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(10)]
	filenames_input  = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_edited = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(3) for id_f in range(10)]
	dir_target       = 'small_colony'
	#'''
	
	#i = 19
	#i = 9
	i = 19
	#'''
	i = 29
	linear = True
	#'''
	
	'''
	# Velocity
	subdirs    = ['val_{}'.format(i) for i in range(2,14,2)]
	filenames  = ['R2_{}.lammpstrj'.format(str(i).zfill(3)) for i in range(7)]
	filenames_input  = [ os.path.join(d, f) for d in subdirs for f in filenames]
	filenames_edited = [ str(id_d).zfill(2)+'_'+str(id_f).zfill(3) for id_d in range(2,14,2) for id_f in range(7) ]
	dir_target       = 'valency_length'
	i = 24 # val_8\R2_003
	i = 35 # val_12\R2_000
	i = 36 # val_12\R2_001
	i = 38 # val_12\R2_003
	i = 15 # val_6\R2_001
	i = 9  # val_4\R2_002
	i = 13 # val_4\R2_006
	
	i = 9  # val_4\R2_002
	i = 23 # val_8\R2_002
	i = 37 # val_12\R2_002
	
	i = 41 # val_12\R2_006
	i = 27 # val_8\R2_006
	i = 13 # val_4\R2_006
	'''
	
	
	filename_input  = filenames_input[i]
	filename_output = filename_input
	filename_edited = filenames_edited[i]
	
	# Shared part of initialization
	dir_lammpstrj    = os.path.join('..', 'lammpstrj3', dir_target)
	dir_edited_data  = os.path.join('data3',dir_target)
	suffix = 'connectivity_graph'
	
	dir_imgs  = os.path.join('imgs3', dir_target,'prof3d')
	os.makedirs(dir_imgs, exist_ok=True)
	
	
	# Load lammpstrj file
	print("\n"+filename_input)
	sampling_frame = utils.get_num_frames(dir_lammpstrj, filename_input)
	print("The last timeframe was loaded: ", sampling_frame )
	types, positions_grid_coord,ids_molecule, time_stamp = utils.load_lammpstrj( dir_lammpstrj, filename_input, sampling_frame )
	
	# Load graph
	d = utils.load(dir_edited_data, filename_edited, suffix)
	multi_graph = d['multi_graph']
	
	# Centering
	center    = utils.get_center_of_mass(types, positions_grid_coord)
	print('center ', center)
	positions_real_coord = utils.centering(positions_grid_coord, center)
	#
	unique_ids_molecule = np.unique( ids_molecule )
	ids_CaMKII = [id for id in unique_ids_molecule if p.subunits['CaMKII binding site']['id'] in types[ids_molecule == id] ]
	ids_GluN2B = [id for id in unique_ids_molecule if p.subunits['GluN2B binding site']['id']    in types[ids_molecule == id] ]
	ids = {	'CaMKII':ids_CaMKII, \
			'GluN2B':ids_GluN2B}
	
	
	targs = ['CaMKII'] # ['CaMKII','GluN2B','PSD95','Shared PSD95']
	
	
	# Plot
	m = 50
	rng = np.random.default_rng(seed=13)
	
	
	fig = plt.figure(figsize=(10, 7)) # , tight_layout=True
	fig.subplots_adjust(wspace=0.6, hspace=0.6)
	for i, t in enumerate(targs):
		ax = fig.add_subplot(1, len(targs), i+1, projection='3d')
		ax.set_title('{} from {}'.format(t, filename_edited))
		#cmap =  matplotlib.colormaps.get_cmap('hsv', 255)
		for id in ids[t]:
			pos = positions_real_coord[ids_molecule == id]
			dist = np.max( np.linalg.norm(pos, axis=1) )
			if dist < m:
				if rng.random() > 0.75:
					if t == 'CaMKII' and linear == False:

						ax.plot(pos[0,0], pos[0,1],pos[0,2],'o', \
							markersize= 2.0, \
							markerfacecolor = 'k',\
							color = 'k')

						for i in range(pos.shape[0]-1):
							#'''
							ax.plot([pos[0,0],pos[i+1,0]], [pos[0,1], pos[i+1,1]], [pos[0,2],pos[i+1,2]],\
								'-' ,\
								color = c.cmap_universal_ratio[t], \
								linewidth = 0.5)
							#'''
							ax.plot(pos[i+1,0], pos[i+1,1],pos[i+1,2],'o', \
								markersize= 2.0, \
								markerfacecolor = 'w',\
								markeredgecolor = c.cmap_universal_ratio[t],\
								color = c.cmap_universal_ratio[t] ) # cmap(id % 256)

					elif t == 'CaMKII' and linear == True:
						for i in range(pos.shape[0]):
							'''
							ax.plot(pos[i,0], pos[i,1],pos[i,2],'o', \
								markersize= 2.0, \
								markerfacecolor = 'w',\
								markeredgecolor = c.cmap_universal_ratio[t],\
								color = c.cmap_universal_ratio[t] ) # cmap(id % 256)
							'''
							if i < pos.shape[0] - 1:
								ax.plot([pos[i,0],pos[i+1,0]], [pos[i,1], pos[i+1,1]], [pos[i,2],pos[i+1,2]],\
									'-' ,\
									color = c.cmap_universal_ratio[t], \
									linewidth = 0.5)


					elif t in ['PSD95','GluN2B']:
						
						ax.plot(pos[:,0], pos[:,1], pos[:,2], 'o-' ,\
							color = c.cmap_universal_ratio[t], \
							linewidth = 0.5, markersize= 2.0)
					elif t in ['Unshared PSD95', 'Shared PSD95']:
						
						ax.plot(pos[:,0], pos[:,1], pos[:,2], 'o-' ,\
							color = c.cmap_universal_ratio['PSD95'], \
							linewidth = 0.5, markersize= 2.0)
						
						
	

		ax.set_xlim(-m,m)
		ax.set_ylim(-m,m)
		ax.set_zlim(-m,m)
		ax.set_box_aspect([1,1,1])
		ax.set_xlabel('(l.u.)')
		ax.set_ylabel('(l.u.)')
		ax.set_zlabel('(l.u.)')	

	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename_edited ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename_edited ) ) , dpi=150)
	plt.show()
	
	
	
	'''
	locs  = pos[1:13,:]
	dists = distance.cdist(locs, locs, 'euclidean')
	for l in range(12):
		dists[l,l] += 100
	for s in range(12):
		j,k = np.unravel_index(np.argsort(dists, axis=None)[s], dists.shape)
		ax.plot([locs[k,0],locs[j,0]], [locs[k,1], locs[j,1]], [locs[k,2],locs[j,2]],\
				'-' ,\
				color = c.cmap_universal_ratio[t], \
				linewidth = 0.5)
	'''
