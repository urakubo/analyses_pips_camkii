
import os, sys, glob, pickle, pprint
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import utils
import parameters as p
import colormap as c

	
	
def get_ids_PSD95_shared_by_STG_GluN2B(multi_graph):
	
	species = 'PSD95'
	
	ids = [ i for i, attr in multi_graph.nodes('species') if attr == species ]
	
	edges_from_molecules = [ multi_graph.edges(id, keys=True, data=True) for id in ids ]
	connections_from_one_molecules = [ [e[3]['type_connection'] for e in es] for es in edges_from_molecules ]
	
	ids_molecule = [id for id, c in zip(ids, connections_from_one_molecules) if ('STG_PSD95' in c) and ('GluN2B_PSD95' in c)]
	
	ids_bead = [multi_graph.nodes[id]['ids_bead'] for id in ids_molecule]
	ids_bead = np.ravel(ids_bead).tolist()
	
	return ids_molecule, ids_bead
	
	
if __name__ == '__main__':
	
	
	#'''
	target = 'CGSP' #  ['CGSP', 'CG','CPG','PG','SP'] # ,'SPG'
	filename_input 	= '{}.lammpstrj'.format(target)
	dir_target      = 'mixtures'
	filename_edited = target
	filename_output = target
	#'''
	
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
	# Find ids_bead for the PSD95 shared by STG and GluN2B
	ids_molecule_shared_PSD, ids_bead_shared_PSD = get_ids_PSD95_shared_by_STG_GluN2B(multi_graph)
	
	# Centering
	center    = utils.get_center_of_mass(types, positions_grid_coord)
	positions_real_coord = utils.centering(positions_grid_coord, center)
	#
	unique_ids_molecule = np.unique( ids_molecule )
	ids_CaMKII = [id for id in unique_ids_molecule if p.subunits['CaMKII hub']['id'] in types[ids_molecule == id] ]
	ids_GluN2B = [id for id in unique_ids_molecule if p.subunits['GluN2B binding site']['id']    in types[ids_molecule == id] ]
	ids_STG    = [id for id in unique_ids_molecule if p.subunits['STG binding site']['id']    in types[ids_molecule == id] ]
	ids_PSD    = [id for id in unique_ids_molecule if p.subunits['PSD']['id']    in types[ids_molecule == id] ]
	ids_unshared_PSD  = set( ids_PSD )
	ids_unshared_PSD.discard(set(ids_molecule_shared_PSD))
	ids = {	'CaMKII':ids_CaMKII, \
			'GluN2B':ids_GluN2B ,\
			'STG':ids_STG,    \
			'PSD95': ids_PSD ,\
			'Shared PSD95': ids_molecule_shared_PSD, \
			'Unshared PSD95': list(ids_unshared_PSD)}
	
	
	# Plot
	m = 30
	rng = np.random.default_rng(seed=13)
	
	fig = plt.figure(figsize=(10, 7)) # , tight_layout=True
	fig.subplots_adjust(wspace=0.6, hspace=0.6)
	for i, t in enumerate(['PSD95','Unshared PSD95','Shared PSD95']): # ['CaMKII','GluN2B','PSD95','Shared PSD95']
		ax = fig.add_subplot(1,3,i+1, projection='3d')
		ax.set_title(t)
		for id in ids[t]:
			pos = positions_real_coord[ids_molecule == id]
			dist = np.max( np.linalg.norm(pos, axis=1) )
			if dist < m:
				if rng.random() > 0.9:
					if t == 'CaMKII':
						for i in range(12):
							ax.plot([pos[0,0],pos[i+1,0]], [pos[0,1], pos[i+1,1]], [pos[0,2],pos[i+1,2]],\
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

	fig.savefig( os.path.join(dir_imgs, '{}.svg'.format( filename_output ) ) )
	fig.savefig( os.path.join(dir_imgs, '{}.png'.format( filename_output ) ) , dpi=150)
	plt.show()
