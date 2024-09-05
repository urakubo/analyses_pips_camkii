import os, glob, pickle, pprint, copy
import numpy as np
from skimage import morphology
from skimage import measure
from skimage.segmentation import watershed

from scipy import ndimage

import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1


# from sklearn import mixture
#https://matsci.org/t/compatibility-issue-between-python-ovito-library-and-matplotlib/50794
os.environ['OVITO_GUI_MODE'] = '1'
from ovito.io import import_file

import lib.parameters as p
import lib.colormap as c


plt.rcParams.update(p.rc_param)


############### Save, load, and decoding
	
def save(dir_data, prefix, suffix, data):
	filename = prefix + '_'+suffix+'.pickle'
	os.makedirs(dir_data, exist_ok=True)
	with open(os.path.join(dir_data, filename), mode='wb') as f:
		pickle.dump(data, f)
	return
	
	
def load(dir_data, prefix, suffix):
	filename = prefix + '_'+suffix+'.pickle'
	with open(os.path.join(dir_data, filename), mode='rb') as f:
		data = pickle.load(f)
	return data
	
	
def get_num_frames(dir_data, filename_data):
	data_all            = import_file(os.path.join(dir_data, filename_data), input_format= "lammps/dump" )
	num_frames          = data_all.source.num_frames
	return num_frames
	
	
def get_num_frames_mc_steps(dir_data, filename_data):
	data_all   = import_file(os.path.join(dir_data, filename_data), input_format= "lammps/dump" )
	num_frames = data_all.source.num_frames
	mc_steps   = data_all.compute(num_frames).attributes['Timestep']
	return num_frames, mc_steps
	
	
def load_data(dir_data, filename_data, id_frame):
	print('Load data.')
	data_all            = import_file(os.path.join(dir_data, filename_data), input_format= "lammps/dump" )
	data_target_frame   = data_all.compute(id_frame)
	types, positions, id_molecule = decode_data(data_target_frame)
	return types, positions, id_molecule
	
	
def load_lammpstrj(dir_data, filename_data, id_frame):
	print('Load types, positions, ids_molecule, time_stamp.')
	data_all            = import_file(os.path.join(dir_data, filename_data), input_format= "lammps/dump" )
	data_target_frame   = data_all.compute(id_frame)
	types, positions, ids_molecule = decode_data(data_target_frame)
	mc_steps 			= data_all.compute(id_frame).attributes['Timestep']
	# data_target_frame.particles['bp']
	
	return types, positions, ids_molecule, mc_steps
	
	
def load_lammpstrj_binding_partners(dir_data, filename_data, id_frame):
	print('Load binding partners.')
	data_all            = import_file(os.path.join(dir_data, filename_data), input_format= "lammps/dump" )
	data_target_frame   = data_all.compute(id_frame)
	bp 	= np.array( data_target_frame.particles['bp'] ).astype('int')
	return bp
	
	
def load_lammpstrj_binding_energy(dir_data, filename_data, id_frame):
	print('Load binding energies.')
	data_all            = import_file(os.path.join(dir_data, filename_data), input_format= "lammps/dump" )
	data_target_frame   = data_all.compute(id_frame)
	energy_isotropic        = np.array( data_target_frame.particles['energy_isotropic'] )
	energy_anisotropic      = np.array( data_target_frame.particles['energy_anisotropic'] )
	energy_anisotropic_self = np.array( data_target_frame.particles['energy_anisotropic_self'] )
	energy = {\
		'energy_isotropic'	: energy_isotropic,\
		'energy_anisotropic': energy_anisotropic,\
		'energy_anisotropic_self': energy_anisotropic_self}
	return energy
	
	
def decode_data(data_frame):
	type         = np.array( data_frame.particles['Particle Type'] )
	position     = np.array( data_frame.particles['Position'] )
	id_molecule  = np.array( data_frame.particles['Molecule Identifier'] )
	# data_target_frame.particles['bp']
	return type, position, id_molecule
	
	
def decode_species(types, positions):
	types_binary = {k: [True if t in v['id'] else False for t in types] for k, v in p.molecules_with_all.items() }
	types_positions = {k: positions[types_binary[k],:] for k in p.molecules_with_all.keys() }
	return types_positions
	
	
	
############### Edit data
	
	## It provides the same operation but offers shorter time: 
	## ids_loc = np.nonzero( get_hist(locs_subunits) )
def get_locs_in_grid_coord(locs):
	locs_int = np.floor(locs).astype(int)
	loc0 = np.array( [p.edge0.index(i) for i in locs_int[:,0]] )
	loc1 = np.array( [p.edge1.index(i) for i in locs_int[:,1]] )
	loc2 = np.array( [p.edge2.index(i) for i in locs_int[:,2]] )
	locs_in_grid_space = np.vstack([loc0, loc1, loc2]).T
	return locs_in_grid_space


def get_periphery_in_grid_mesh():
	radius = (np.min(p.space_np) - 1) / 2
	ball   = morphology.ball(radius)
	periphery = (ball == 0)
	return periphery
	

def get_local_mins(data, mask = None):
	
	axes = [ 0, 1, 2, (0,1),  (0,1), (1,2), (1,2) , (2,0), (2, 0), (0,1,2), (0,1,2) , (0,1,2) , (0,1,2)]
	vecs = [ 1, 1, 1, (1,1), (1,-1), (1,1), (1,-1), (1,1), (1,-1), (1,1,1), (1,1,-1), (1,-1,1), (-1,1,1)]
	diff = np.ones_like(data, dtype=bool)
	for a, v in zip(axes, vecs):
		vv= np.array(v)
		data1 = np.roll(data,   vv, axis=a)
		data2 = np.roll(data,  -vv, axis=a)
		diff  = diff&(data1-data >= 0)&(data2-data >= 0)
		#diff  = diff&(data1-data > 0)&(data2-data > 0)
	
	if mask is not None:
		diff  = diff*mask
	return diff
	
	
def get_high(conc, th = 0.5):
	return (conc > np.max(conc)*th)
	
	
def get_hist(positions):

	H, (xedges, yedges, zedges) = np.histogramdd( positions, bins=(p.edges0, p.edges1, p.edges2) )
	return H
	
	
def get_sum_energy(positions, energy):
	H, (xedges, yedges, zedges) = np.histogramdd( positions, bins=(p.edges0, p.edges1, p.edges2), weights=energy )
	return H
	
	
def get_min_local_mins(local_mins, conc_smooth, mask):
	num_local_mins   = np.sum( local_mins )
	if num_local_mins != 0:
		locs_local_min   = np.nonzero(local_mins)
		local_min_values = conc_smooth[locs_local_min]
		j = np.argmin( local_min_values )
		loc_local_min    = [ locs_local_min[0][j], locs_local_min[1][j], locs_local_min[2][j] ]
		local_min_value  = local_min_values[j]
	else:
		loc_local_min    = None
		local_min_value  = np.min( conc_smooth[mask] )
	
	return {'location': loc_local_min, 'value':local_min_value}
	
	
def get_rotation_matrix_in_CaMKII_PSD95_direction( locs_in_real_coord ): 
	# Rotation
	# https://stackoverflow.com/questions/14607640/rotating-a-vector-in-3d-space
	
	targs_molecule  = p.molecules_with_all.keys()
	ids_target = {t: np.linalg.norm( locs_in_real_coord[t], axis = 1 ) < np.min(p.center_np) for t in targs_molecule}
	targets_locs_in_real_coord = {t: locs_in_real_coord[t][ids_target[t]] for t in targs_molecule}
	
	directions = {t: np.mean(targets_locs_in_real_coord[t], axis = 0) for t in targs_molecule}
	#print(directions)
	direction = directions['PSD95'] - directions['CaMKII']
	direction = direction / np.linalg.norm(direction)
	
	x = direction[0]
	y = direction[1]
	z = direction[2]
	x2_y2= np.sqrt(x*x+y*y)
	theta_xy = np.arctan2(y, x)
	theta_xz = np.arctan2(x2_y2, z)
	r1 = np.array([[np.cos(theta_xy), np.sin(theta_xy), 0],[-np.sin(theta_xy), np.cos(theta_xy), 0],[0,0,1]])
	r2 = np.array([[np.cos(theta_xz), 0, -np.sin(theta_xz)],[0, 1, 0],[np.sin(theta_xz), 0, np.cos(theta_xz)]])
	rot_matrix = np.dot(r2, r1)
	return rot_matrix
	
	
	
def get_ids_PSD95_shared_by_STG_GluN2B(multi_graph, shared_or_unshared = 'shared'):
	
	species = 'PSD95'
	
	ids = [ i for i, attr in multi_graph.nodes('species') if attr == species ]
	
	edges_from_molecules = [ multi_graph.edges(id, keys=True, data=True) for id in ids ]
	connections_from_one_molecules = [ [e[3]['type_connection'] for e in es] for es in edges_from_molecules ]
	
	if shared_or_unshared == 'unshared':
		ids_molecule = [id for id, c in zip(ids, connections_from_one_molecules) if ('STG_PSD95' not in c) or ('GluN2B_PSD95' not in c)]
	else:
		ids_molecule = [id for id, c in zip(ids, connections_from_one_molecules) if ('STG_PSD95' in c) and ('GluN2B_PSD95' in c)]
	
	ids_bead = [multi_graph.nodes[id]['ids_bead'] for id in ids_molecule]
	ids_bead = np.ravel(ids_bead).tolist()
	
	return ids_molecule, ids_bead
	
	
def get_ids_PSD95_shared_only_by_GluN2B(multi_graph):
	
	species = 'PSD95'
	
	ids = [ i for i, attr in multi_graph.nodes('species') if attr == species ]
	
	edges_from_molecules = [ multi_graph.edges(id, keys=True, data=True) for id in ids ]
	connections_from_one_molecules = [ [e[3]['type_connection'] for e in es] for es in edges_from_molecules ]
	
	ids_molecule = [id for id, c in zip(ids, connections_from_one_molecules) if ('STG_PSD95' not in c) and ('GluN2B_PSD95' in c)]
	
	ids_bead = [multi_graph.nodes[id]['ids_bead'] for id in ids_molecule]
	ids_bead = np.ravel(ids_bead).tolist()
	
	return ids_molecule, ids_bead
	
	
	#
	# grid_coord: Molecular coordinate in the grid space (0, 1,..., 119) 
	# (numpy uint) [(x0,y0,z0), (x1,y1,z1), ..., (xn,yn,zn)]
	#
	# real_coord: Molecular coordinate in the real space [-60, 60) 
	# (numpy float) [(x0,y0,z0), (x1,y1,z1), ..., (xn,yn,zn)]
	#
	# grid_mesh : molecular existance in the 3D space
	# (numpy bool in 3D space) (p.space[0], p.space[1], p.space[2])
	#
	
def get_concs_and_condensates(types, positions, ids_molecule, multi_graph = None, sigma=2):
	
	# Parameters
	targs_molecule  = p.molecules_with_all.keys() # ['GluN2B', 'CaMKII', 'STG', 'PSD95', 'All']
	
	# Get the locations of target molecules separately.
	flag_type          = {t: [True if k in p.molecules_with_all[t]['id'] else False for k in types] for t in targs_molecule}
	locs_in_real_coord = {k: positions[flag_type[k],:] for k in targs_molecule } 
	locs_in_grid_mesh  = {k: get_hist(locs_in_real_coord[k]) for k in targs_molecule}
	
	# Get a peripheral region (a region outside of a sphere) 
	# and obtain the concentrations in this area (diluted region concentration for a baseline).
	region_periphery = get_periphery_in_grid_mesh()
	concs_periphery  = {t: np.sum(locs_in_grid_mesh[t] * region_periphery ) / np.sum( region_periphery ) for t in targs_molecule}
	
	# Get the condenate regions of targs_molecule.
	concs_in_grid_mesh = {t: ndimage.gaussian_filter(locs_in_grid_mesh[t], sigma = sigma) for t in targs_molecule}
	regions_condensate_in_grid_mesh = {t: get_high(concs_in_grid_mesh[t]-concs_periphery[t]) for t in targs_molecule}
	regions_condensate_in_grid_mesh['dilute'] = region_periphery
	
	def get_concs_condensate(ref_molecule):
		rs     = {t: np.sum( locs_in_grid_mesh[t] * regions_condensate_in_grid_mesh[ref_molecule] ) for t in targs_molecule }
		rs_sum = {t: np.sum( regions_condensate_in_grid_mesh[ref_molecule] ) for t in targs_molecule }
		concs_region = {t: rs[t] / rs_sum[t] if rs_sum[t] != 0 else rs[t] for t in targs_molecule }
		return concs_region
	concs_condensate = {t: get_concs_condensate(t) for t in targs_molecule}
	
	
	# Summary
	d = {
			'locs_in_grid_mesh':	locs_in_grid_mesh,
			'concs_in_grid_mesh':	concs_in_grid_mesh,
			'region_condensate_in_grid_mesh': regions_condensate_in_grid_mesh,
			#
			'conc_periphery'	:	concs_periphery,
			'conc_condensate'	:	concs_condensate,
		}
	
	
	# Optionals
	# Rotate particles and obtain their concentrations
	if any( flag_type['CaMKII'] ) and any( flag_type['PSD95'] ) :
		print('Rotated')
		rot_matrix = get_rotation_matrix_in_CaMKII_PSD95_direction( locs_in_real_coord )
		rotated_in_real_coord = {t: rot_matrix.dot(locs_in_real_coord[t].T).T  for t in targs_molecule}
	
		rotated_in_grid_mesh  = {k: get_hist(rotated_in_real_coord[k]) for k in targs_molecule}
		rotated_concs_in_grid_mesh = \
			{t: ndimage.gaussian_filter(rotated_in_grid_mesh[t], sigma = sigma) for t in targs_molecule}
		rotated_regions_condensate_in_grid_mesh = \
			{t: get_high(rotated_concs_in_grid_mesh[t]-concs_periphery[t]) for t in targs_molecule}
	
		d['rotated_concs_in_grid_mesh'] = rotated_concs_in_grid_mesh
		d['rotated_region_condensate_in_grid_mesh'] = rotated_regions_condensate_in_grid_mesh
	
	
	if multi_graph is not None:
		_, ids_bead_shared_PSD   = get_ids_PSD95_shared_by_STG_GluN2B(multi_graph, shared_or_unshared = 'shared')
		_, ids_bead_unshared_PSD = get_ids_PSD95_shared_by_STG_GluN2B(multi_graph, shared_or_unshared = 'unshared')
		_, ids_bead_PSD_shared_by_GluN2B = get_ids_PSD95_shared_only_by_GluN2B(multi_graph)
		
		for m, ids_bead in zip( ['Shared PSD95', 'Unshared PSD95', 'PSD95 shared only by GluN2B'], \
								[ids_bead_shared_PSD, ids_bead_unshared_PSD, ids_bead_PSD_shared_by_GluN2B] ):
			loc_in_grid_mesh  = get_hist(positions[ids_bead,:])
			conc_in_grid_mesh = ndimage.gaussian_filter(loc_in_grid_mesh, sigma = sigma)
			d['locs_in_grid_mesh'][m]  = loc_in_grid_mesh
			d['concs_in_grid_mesh'][m] = conc_in_grid_mesh
			
			rotated_in_real_coord = rot_matrix.dot(positions[ids_bead,:].T).T
			rotated_in_grid_mesh  = get_hist(rotated_in_real_coord)
			rotated_conc_in_grid_mesh = ndimage.gaussian_filter(rotated_in_grid_mesh, sigma = sigma)
			d['rotated_concs_in_grid_mesh'][m] = rotated_conc_in_grid_mesh
		
	return d
	
	
############### Watershed segmentation
	
def watershed_segmentation( d ):
	markers = np.zeros(p.space, dtype = int)
	loc = int(p.space[0]/2), int(p.space[1]/2), int(p.space[2]/2)
	markers[ loc ]    = 1
	markers[10,10,10] = 2
	
	tot_volume = p.space[0] * p.space[1] * p.space[2]
	targets_watershed = ['CaMKII', 'STG']
	
	labels_watershed_in_grid_mesh = {}
	ratio_volumes_watershed = {}
	for j, t_watershed in enumerate( targets_watershed ):
		concs_in_grid_mesh = d['concs_in_grid_mesh'][t_watershed]
		label_watershed = watershed(concs_in_grid_mesh, markers=markers) == 1
		labels_watershed_in_grid_mesh[t_watershed] = label_watershed
		ratio_volumes_watershed[t_watershed] = np.sum( label_watershed ) / tot_volume
	return labels_watershed_in_grid_mesh, ratio_volumes_watershed
	
	
	
############### Centering

def get_center_of_mass_simple(positions):
	
	# print('position_ref: ', position_ref)
	x = (positions[:,0] / p.space[0]) * 2 * np.pi
	y = (positions[:,1] / p.space[1]) * 2 * np.pi
	z = (positions[:,2] / p.space[2]) * 2 * np.pi
	x0, x1 = np.cos(x), np.sin(x)
	y0, y1 = np.cos(y), np.sin(y)
	z0, z1 = np.cos(z), np.sin(z)

	x0, x1 = np.mean(x0), np.mean(x1)
	y0, y1 = np.mean(y0), np.mean(y1)
	z0, z1 = np.mean(z0), np.mean(z1)

	theta_x = np.arctan2(x1, x0)
	theta_y = np.arctan2(y1, y0)
	theta_z = np.arctan2(z1, z0)

	center_of_mass = np.array( [theta_x * p.space[0], theta_y * p.space[1] , theta_z * p.space[2] ] ) + np.pi
	center_of_mass /= (2 * np.pi)

	for dim in [0,1,2]:
		center_of_mass[dim] += \
						- p.space[dim] * (center_of_mass[dim]  >=  p.space[dim]) \
						+ p.space[dim] * (center_of_mass[dim]  <  0)

	return center_of_mass
	


def get_center_of_mass(types_, positions_, reference_molecule_for_centering = p.reference_molecule_for_centering):
	
	types = [True if t in p.molecules_with_all[reference_molecule_for_centering]['id'] else False for t in types_ ]
	position_ref = positions_[types,:]
	
	# print('position_ref: ', position_ref)
	x = (position_ref[:,0] / p.space[0]) * 2 * np.pi
	y = (position_ref[:,1] / p.space[1]) * 2 * np.pi
	z = (position_ref[:,2] / p.space[2]) * 2 * np.pi
	x0, x1 = np.cos(x), np.sin(x)
	y0, y1 = np.cos(y), np.sin(y)
	z0, z1 = np.cos(z), np.sin(z)

	x0, x1 = np.mean(x0), np.mean(x1)
	y0, y1 = np.mean(y0), np.mean(y1)
	z0, z1 = np.mean(z0), np.mean(z1)

	theta_x = np.arctan2(x1, x0)
	theta_y = np.arctan2(y1, y0)
	theta_z = np.arctan2(z1, z0)

	center_of_mass = np.array( [theta_x * p.space[0], theta_y * p.space[1] , theta_z * p.space[2] ] ) + np.pi
	center_of_mass /= (2 * np.pi)

	return center_of_mass

	# Reference
	# https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions


def centering(positions, center): 
	p_centered = positions - center
	for dim in [0,1,2]:
		over  = (p_centered[:,dim]  >  p.space[dim]/2)
		under = (p_centered[:,dim] <= -p.space[dim]/2)
		p_centered[over ,dim] -= p.space[dim]
		p_centered[under,dim] += p.space[dim]
	return p_centered


def centering_lattice_space(positions, center): 
	p_centered = (positions - center) + p.center_np.astype('int')
	for dim in [0,1,2]:
		over  = (p_centered[:,dim]  >=  p.space[dim])
		under = (p_centered[:,dim]  <  0)
		p_centered[over ,dim] -= p.space[dim]
		p_centered[under,dim] += p.space[dim]
	return p_centered


def get_average_center_of_mass( dir_lammpstrj, filename_input, sampling_time_frames ):
	centers = np.zeros( (len(sampling_time_frames), 3) )
	
	for i, target_frame in enumerate( sampling_time_frames ):
		types, positions_grid_coord,ids_molecule = load_data( dir_lammpstrj, filename_input, target_frame )
		positions_real_coord = centering(positions_grid_coord, p.center_np)
		centers[i,:] = get_center_of_mass(types, positions_real_coord) - p.center_np
	
	center = np.mean(centers, axis = 0)
	center -= p.space_np * (center >  p.center_np)
	center += p.space_np * (center <= -p.center_np)
	
	return center


############### Radial distrbution function

def get_lattice_grids():
	x, y, z = np.meshgrid(p.edge0, p.edge1, p.edge2)
	x, y, z = np.ravel(x), np.ravel(y), np.ravel(z)
	grids   = np.vstack((x,y,z)).T
	return grids


def get_a_rdf(types, positions, rdf_grid_points, center = None, multi_graph=None):
	
	# Centering
	if center is None:
		center = get_center_of_mass(types, positions)
	
	positions           = centering(positions, center)
	positions_grid_centered = centering(rdf_grid_points, center)
	
	# Decode species
	types_positions     = decode_species(types, positions)
	
	if multi_graph is not None:
		_, ids_bead_shared_PSD   = get_ids_PSD95_shared_by_STG_GluN2B(multi_graph, shared_or_unshared = 'shared')
		_, ids_bead_unshared_PSD = get_ids_PSD95_shared_by_STG_GluN2B(multi_graph, shared_or_unshared = 'unshared')
		types_positions['Shared PSD95']   = positions[ids_bead_shared_PSD,:]
		types_positions['Unshared PSD95'] = positions[ids_bead_unshared_PSD,:]
	
	# Get distances from the center
	dists = {k: np.linalg.norm(v, axis=1) for k, v in types_positions.items()}
	dists_grid = np.linalg.norm(positions_grid_centered, axis=1)

	# Get radial distribution function (rdf)
	num_grid_around_center, _  = np.histogram(dists_grid , bins=p.rdf_bins)
	rdf = {}
	for k, v in dists.items():
		num_molecule_around_center, _  = np.histogram(v , bins=p.rdf_bins)
		rdf[k] = num_molecule_around_center / num_grid_around_center
	
	return rdf

def get_rdf_CaMKII( dir_lammpstrj, filename_input, sampling_frame ):

	# Load data
	types, positions, _ = load_data( dir_lammpstrj, filename_input, sampling_frame )
	
	# Intial setting
	rdf_grid_points         = get_lattice_grids()
	center = get_center_of_mass(types, positions, reference_molecule_for_centering = 'CaMKII')
	print('get_rdf_CaMKII center: ', center)
	
	positions               = centering(positions, center)
	positions_grid_centered = centering(rdf_grid_points, center)
	
	type_CaMKII = [True if t in p.molecules_with_all['CaMKII']['id'] else False for t in types]
	position_CaMKII = positions[type_CaMKII,:]
	
	
	# Get distances from the center
	dists_CaMKII = np.linalg.norm(position_CaMKII, axis=1)
	dists_grid   = np.linalg.norm(positions_grid_centered, axis=1)
	
	
	# Get radial distribution function (rdf)
	num_grid_around_center, _      = np.histogram(dists_grid   , bins=p.rdf_bins)
	num_molecule_around_center, _  = np.histogram(dists_CaMKII , bins=p.rdf_bins)
	rdf = {}
	rdf['CaMKII_bead'] = num_molecule_around_center / num_grid_around_center
	
	return rdf, p.rdf_bins
	

def get_rdf_CaMKII_hub_binding( dir_lammpstrj, filename_input, sampling_frame ):

	# Load data
	types, positions, _ = load_data( dir_lammpstrj, filename_input, sampling_frame )
	
	# Intial setting
	rdf_grid_points         = get_lattice_grids()
	center = get_center_of_mass(types, positions, reference_molecule_for_centering = 'CaMKII')
	print('get_rdf_CaMKII center: ', center)
	
	positions               = centering(positions, center)
	positions_grid_centered = centering(rdf_grid_points, center)
	
	type_CaMKII_hub     = [True if t == p.subunits['CaMKII hub']['id'] else False for t in types]
	type_CaMKII_binding = [True if t == p.subunits['CaMKII binding site']['id'] else False for t in types]
	position_CaMKII_hub = positions[type_CaMKII_hub,:]
	position_CaMKII_binding = positions[type_CaMKII_binding,:]
	
	
	# Get distances from the center
	dists_CaMKII_hub = np.linalg.norm(position_CaMKII_hub, axis=1)
	dists_CaMKII_binding = np.linalg.norm(position_CaMKII_binding, axis=1)
	dists_grid   = np.linalg.norm(positions_grid_centered, axis=1)
	
	
	# Get radial distribution function (rdf)
	num_grid_around_center, _      = np.histogram(dists_grid   , bins=p.rdf_bins)
	num_CaMKII_binding_around_center, _  = np.histogram(dists_CaMKII_binding , bins=p.rdf_bins)
	num_CaMKII_hub_around_center, _      = np.histogram(dists_CaMKII_hub , bins=p.rdf_bins)
	rdf = {}
	rdf['CaMKII_binding'] = num_CaMKII_binding_around_center / num_grid_around_center
	rdf['CaMKII_hub'] = num_CaMKII_hub_around_center / num_grid_around_center
	
	return rdf, p.rdf_bins
	
	
def get_rdfs_from_multiple_frames( dir_lammpstrj, filename_input, sampling_time_frames, center = None, multi_graph=None ):
	
	# Parameters
	rdf_grid_points   = get_lattice_grids()
	
	if center is None:
		center = get_average_center_of_mass( dir_lammpstrj, filename_input, sampling_time_frames )
	
	rdfs = { k: np.zeros( ( len(p.rdf_bins)-1, len(sampling_time_frames) ) ) for k in p.molecules_with_all.keys() }
	if multi_graph is not None:
		for targ in p.rdf_targets_multi_graph:
			rdfs[targ] = np.zeros( ( len(p.rdf_bins)-1, len(sampling_time_frames) ) )
			
	for i, id_frame in enumerate( sampling_time_frames ):
		types, positions, _, time_stamp = load_lammpstrj( dir_lammpstrj, filename_input, id_frame )
		print('Time_stamp (MC steps): ', time_stamp)
		current_rdfs = get_a_rdf(types, positions, rdf_grid_points, center, multi_graph=multi_graph )
		for k in rdfs.keys():
			rdfs[k][:,i] = current_rdfs[k]
		
	return rdfs, p.rdf_bins
	
	
def get_rdfs( dir_input, filename_input, target_frame, center=None, multi_graph=None ):
	
	# Target frames
	rdf_sampling_frames = list(range( target_frame - (p.rdf_num_sampling_frames - 1)*p.rdf_sampling_interval, \
			target_frame + p.rdf_sampling_interval,\
			p.rdf_sampling_interval))
	
	if np.any(np.array(rdf_sampling_frames) < 0):
		rdf_sampling_frames = [0]

	print('rdf_sampling_frames ', rdf_sampling_frames)
	rdfs, p.rdf_bins = get_rdfs_from_multiple_frames(dir_input, filename_input, \
			rdf_sampling_frames, center = center, multi_graph=multi_graph)
	
	return rdfs, p.rdf_bins, rdf_sampling_frames
	
	
def plot_a_rdf( ax, d, errorbar='shaded', legend=True, target_molecules = p.molecules_without_all.keys() , ylim = (-0.006,0.66), pmax = False  ):
	
	if pmax == True:
		ylim = (-5.0,150)
	
	r = d['rdf_bins'][1:-1] / np.sqrt(3) ###
	for k in target_molecules:
		rdf_mean  = np.mean( d['rdf'][k][1:], axis = 1 )
		rdf_std   = np.std(  d['rdf'][k][1:], axis = 1 )
		
		if pmax == True:
			rdfmax = np.max( rdf_mean )
			rdf_mean  = rdf_mean / rdfmax * 100
			rdf_std   = rdf_std / rdfmax * 100
		
		color     = c.cmap_universal_ratio[k]
		if errorbar   == 'shaded':
			ax.fill_between(r, rdf_mean-rdf_std, rdf_mean+rdf_std,color=c.cmap_universal_ratio_light[k])
			ax.plot(r, rdf_mean, color=color, label=k)
		elif errorbar == 'line':
			'''
			ax.hist(r[:-1], r, weights=rdf_mean[:-1], \
				color=color, \
				histtype="step", label=k)
			'''
			ax.step(r, rdf_mean, where='post', color=color, label=k)
			ax.errorbar(r+0.5, rdf_mean, yerr=rdf_std,\
				ecolor=color,\
				linestyle='',
				alpha = 0.4 )
		elif errorbar == 'final_frame_alone':
			ax.step( r, d['rdf'][k][1:,-1], color=color, label=k)
	if legend==True:
		ax.legend(frameon=False)

	'''
	ax.set_xlabel('Distance from \n center-of-mass (l.u.)')
	ax.set_xticks(range(0,50,10))
	ax.set_xlim(0,40)
	'''
	if pmax == True:
		ax.set_ylabel('(% max)')
	else:
		ax.set_ylabel('(beads / voxel)')
	ax.set_xlabel('Distance from \n center-of-mass (l.s.)')
	ax.set_xticks(range(0,25,5))
	ax.set_xlim(0,23)
	
	ax.set_ylim(*ylim)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	return


def plot_a_rdf_PSD95( ax, d, legend=True , ylim = (-0.006,0.66) ):
	r = d['rdf_PSD95_bins'][1:-1] / np.sqrt(3) ###
	target_molecules = ['CaMKII', 'GluN2B','PSD95', 'Shared PSD95', 'STG']
	for k in target_molecules:
		'''
		rdf_mean  = np.mean( d['rdf_PSD95'][k][1:], axis = 1 )
		rdf_std   = np.std(  d['rdf_PSD95'][k][1:], axis = 1 )
		color     = c.cmap_universal_ratio[k]

		ax.fill_between(r, rdf_mean-rdf_std, rdf_mean+rdf_std,color=c.cmap_universal_ratio_light[k])
		ax.plot(r, rdf_mean, color=color, label=k)
		'''
		color     = c.cmap_universal_ratio[k]
		#print(d['rdf_PSD95'][k].shape)
		#print(d['rdf_PSD95'][k][1:,-1])
		ax.plot(r, d['rdf_PSD95'][k][1:,-1], color=color, label=k)

	if legend==True:
		ax.legend(frameon=False)
	
	
	'''
	ax.set_xlabel('Distance from \n center-of-mass (l.u.)')
	ax.set_xticks(range(0,50,10))
	ax.set_xlim(0,40)
	'''
	ax.set_xlabel('Distance from \n center-of-mass (l.s.)')
	ax.set_xticks(range(0,25,5))
	ax.set_xlim(0,23)
	ax.set_ylabel('(beads / voxel)')
	
	ax.set_ylim(*ylim)

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	return

############### Profiles / Each panel plot

def arrange_graph_no_ticks(ax):
	#ax.xaxis.set_tick_params(labelbottom=False)
	#ax.yaxis.set_tick_params(labelleft=False)
	#ax.set_xticks([])
	#ax.set_yticks([])
	
	ax.set_axis_off()
	
def plot_colorbar(ax, cs):
	divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
	cax = divider.append_axes('right', '5%', pad='3%')
	cb = plt.colorbar(cs, cax=cax)
	ticks = copy.copy( cb.get_ticks() ).tolist()
	cb.set_ticks(ticks)
	cb.set_ticklabels(["{:.2f}".format(i) for i in ticks])
	cb.set_label('(beads / voxel)')
	return cb
	
	
def plot_scalebar(ax, col='k', linewidth=2):
	ax.plot([5,5+10*np.sqrt(3)],[5,5], '-', color=col, linewidth=linewidth) ###
	return
	
	
def make_a_panel_of_CaMKII_STG_condenstates(d, transp, slice, pre_rotated=False ):
	# Condensate region
	if pre_rotated == False:
		r_CaMKII   = d['region_condensate_in_grid_mesh']['CaMKII'].transpose(transp)[slice,:,:]
		r_STG      = d['region_condensate_in_grid_mesh']['STG'].transpose(transp)[slice,:,:]
	else:
		r_CaMKII   = d['rotated_region_condensate_in_grid_mesh']['CaMKII'].transpose(transp)[slice,:,:]
		r_STG      = d['rotated_region_condensate_in_grid_mesh']['STG'].transpose(transp)[slice,:,:]
	r_BOTH     = r_CaMKII & r_STG
	# white panel
	panel = np.ones( [r_CaMKII.shape[0],r_CaMKII.shape[1],3], dtype=np.uint8 )*255 
	col = c.cmap_universal_uint['STG']
	for k in range(3): panel[r_STG,k] = col[k]
	col = c.cmap_universal_uint['CaMKII']
	for k in range(3): panel[r_CaMKII,k] = col[k]
	col = [255, 241, 0] # Yellow
	for k in range(3): panel[r_BOTH,k] = col[k]
	return panel
	
	
def make_a_panel_of_CaMKII_All_condenstates(d, transp, slice, pre_rotated=False ):
	# Condensate region
	if pre_rotated == False:
		r_CaMKII   = d['region_condensate_in_grid_mesh']['CaMKII'].transpose(transp)[slice,:,:]
		r_All      = d['region_condensate_in_grid_mesh']['All'].transpose(transp)[slice,:,:]
	else:
		r_CaMKII   = d['rotated_region_condensate_in_grid_mesh']['CaMKII'].transpose(transp)[slice,:,:]
		r_All      = d['rotated_region_condensate_in_grid_mesh']['All'].transpose(transp)[slice,:,:]
	r_BOTH     = r_CaMKII & r_All
	# white panel
	panel = np.ones( [r_CaMKII.shape[0],r_CaMKII.shape[1],3], dtype=np.uint8 )*255 
	col = c.cmap_universal_uint['All']
	for k in range(3): panel[r_All,k] = col[k]
	col = c.cmap_universal_uint['CaMKII']
	for k in range(3): panel[r_CaMKII,k] = col[k]
	col = c.light_green_universal_uint # Lignt green
	for k in range(3): panel[r_BOTH,k] = col[k]
	return panel
	
	
def plot_regions_condenstate_from_a_direction(fig, num_rows, num_columns, row, column, d, transp = (0,1,2), title=True, scalebar=True, pre_rotated=False ):
	slice = int( p.space[0]/2 )
	panel = make_a_panel_of_CaMKII_All_condenstates(d, transp, slice, pre_rotated=pre_rotated )
	ax    = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
	ax.imshow( panel )
	if scalebar==True:
		plot_scalebar(ax)
	if title == True:
		ax.set_title('Green: CaMKII, \n Gray: All')
	arrange_graph_no_ticks(ax)
	return ax
	
	
def plot_regions_condenstate_from_a_direction_(fig, num_rows, num_columns, row, column, d, transp = (0,1,2), title=True, scalebar=True, pre_rotated=False ):
	slice = int( p.space[0]/2 )
	panel = make_a_panel_of_CaMKII_STG_condenstates(d, transp, slice, pre_rotated=pre_rotated )
	ax    = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
	ax.imshow( panel )
	if scalebar==True:
		plot_scalebar(ax)
	if title == True:
		ax.set_title('Green: CaMKII, \n Red: STG')
	arrange_graph_no_ticks(ax)
	return ax
	
	
def plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp = (0,1,2), title=True, colorbar=True, scalebar=True, vmin=0, vmax=None, pre_rotated=False ):
	slice = int( p.space[0]/2 )
	axes = []
	for i, (target, column) in enumerate(columns.items()):
		ax = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
		axes.append(ax)
		
		if pre_rotated == False:
			panel = d['concs_in_grid_mesh'][target].transpose(transp)[slice,:,:]
			cs = ax.imshow( panel , cmap=c.cmap[target], vmin=vmin, vmax=vmax )
			if i == 0 and scalebar==True:
				plot_scalebar(ax, col='w', linewidth=3)
		else:
			panel = d['rotated_concs_in_grid_mesh'][target].transpose(transp)[slice,:,:]
			cmap=c.cmap[target]
			cmap.set_bad(alpha = 0.0)
			mask = morphology.disk(p.center_np[0]-0.5)
			panel[np.logical_not(mask)] = float('nan')
			cs = ax.imshow( panel , cmap=cmap, vmin=vmin, vmax=vmax )
			
			#mask = np.ma.masked_where(mask == 0, mask)
			#ax.imshow(mask,alpha=1,cmap = 'Reds')
			
			if i == 0 and scalebar==True:
				plot_scalebar(ax, col='k', linewidth=3)

		if title is True:
			ax.set_title( target )
		elif title is False:
			pass
		else:
			ax.set_title( title )
		arrange_graph_no_ticks(ax)
		if colorbar == True:
			plot_colorbar(ax, cs)
	return axes
	
	
def plot_watershed_region_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp = (0,1,2), title=True, scalebar=True ):
	slice = int( p.space[0]/2 )
	axes = []
	for i, (target, column) in enumerate(columns.items()):
		panel = d['labels_watershed_in_grid_mesh'][target].transpose(transp)[slice,:,:]
		ax = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
		axes.append(ax)
		cs = ax.imshow( panel, cmap="binary", vmin = 0, vmax = 1.5 )
		if i == 0 and scalebar==True:
			plot_scalebar(ax)
		if title == True:
			ax.set_title('Watershed \n separated by '+ target )
		arrange_graph_no_ticks(ax)
	return axes
	
	
def arrange_graph_bar(ax, panel_dx, y0, panel_size_x, panel_size_y):
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax_h, ax_w = ax.bbox.height, ax.bbox.width
	# ax.set_aspect('auto') # "equal" , ax.set_aspect(1.0 / ax.get_data_ratio())
	#x, y, w, h = ax.get_position().bounds
	loc = ax.get_position()
	if y0 is None:
		y0 = loc.y0
	ax.set_position([loc.x0+panel_dx, y0, panel_size_x, panel_size_y])
	
	
def plot_conc_condensate_bar(ax, targ, ref, d):
	cc = d['conc_condensate']
	col = c.cmap_universal_ratio[targ]
	ax.bar(['In '+targ+'\n condensates', 'In '+ref+'\n condensates'], [cc[targ][targ], cc[ref][targ]], width=0.5, color=col)
	ax.set_title('Conc of {}'.format(targ))
	ax.set_ylabel('(beads / voxel)')
	ax.set_ylim(0,0.6)
	ax.tick_params(axis='x', rotation=45)
	
	
def plot_conc_ratio_condensate_bar(ax, targs, counterparts, d):
	cc = d['conc_condensate']
	conc_ratio = [ cc[t][t] / cc[c][t] for t, c in zip(targs, counterparts)]
	cols = [ c.cmap_universal_ratio[targ] for targ in targs ]
	ax.bar(targs, conc_ratio, width=0.5, color=cols)
	ax.set_title('Partition index')
	ax.set_ylabel('(target / counterpart)')
	ax.set_ylim(0,40)
	
	
def plot_concs_condensate_bar(ax,ref, d):
	data   = { k: d['conc_condensate'][ref][k]  for k in p.molecules_with_all.keys()}
	colormap_conc_bargraph =[c.cmap_universal_ratio[k] for k in data.keys()]
	data['Total'] = data.pop('All') ###
	ax.set_title('Concentration \n in {} condensate'.format( ref ))
	ax.bar(*zip(*data.items()), width=0.6, color=colormap_conc_bargraph )
	ax.set_ylim(0,1.0)
	ax.set_ylabel('(beads / voxel)')
	ax.tick_params(axis='x', rotation=45)
	
	
def plot_binding_energy_bar(ax,d, type_energy, target_condensate, ymax):
	# type_energy: 'energy_anisotropic', 'energy_anisotropic_self', 'energy_isotropic'
	# target_condensate: 'CaMKII', 'STG'
	
	vol   = np.sum( d['region_condensate_in_grid_mesh'][target_condensate] )
	data  = { k: -d[type_energy][target_condensate][k] / vol for k in p.molecules_with_all.keys()}
	colormap_conc_bargraph =[c.cmap_universal_ratio[k] for k in data.keys()]
	data['Total'] = data.pop('All') ###
	ax.set_title('{} per voxel \n in {} cond (Total: {:.3g})'.format( type_energy, target_condensate, data['Total'] ))
	ax.bar(*zip(*data.items()), width=0.6, color=colormap_conc_bargraph )
	ax.set_ylim(0,ymax)
	ax.set_ylabel('(1 / voxel)')
	ax.tick_params(axis='x', rotation=45, horizontalalignment="right")
	
	return data['Total']


def rotate(mesh_CaMKII, mesh_STG): 
	# Rotation
	# https://stackoverflow.com/questions/14607640/rotating-a-vector-in-3d-space
	CaMKII_dir = np.mean(mesh_CaMKII.vertices, axis=0)
	STG_dir    = np.mean(mesh_STG.vertices, axis=0)
	
	direction = STG_dir - CaMKII_dir
	direction = direction / np.linalg.norm(direction)
	
	x = direction[0]
	y = direction[1]
	z = direction[2]
	x2_y2= np.sqrt(x*x+y*y)
	theta_xy = np.arctan2(y, x)
	theta_xz = np.arctan2(x2_y2, z)
	r1 = np.array([[np.cos(theta_xy), np.sin(theta_xy), 0],[-np.sin(theta_xy), np.cos(theta_xy), 0],[0,0,1]])
	r2 = np.array([[np.cos(theta_xz), 0, -np.sin(theta_xz)],[0, 1, 0],[np.sin(theta_xz), 0, np.cos(theta_xz)]])
	rot_matrix = np.eye(4)
	rot_matrix[:3,:3] = np.dot(r2, r1)
	
	# mm = trimesh.transformations.random_rotation_matrix()
	mesh_CaMKII.apply_transform(rot_matrix)
	mesh_STG.apply_transform(rot_matrix)
	return rot_matrix
	
	
def select_plot(target, fig, num_rows, num_columns, row, column, d, title):
	
	scalebar = (row == 0) and (column == 1)
	value = True
	
	if target == 'region_condensates':
		plot_regions_condenstate_from_a_direction(fig, num_rows, num_columns, row, column, d, title=title, scalebar=scalebar )
	elif target == 'shared_PSD95':
		columns = {'Shared PSD95':column}
		plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, \
			title=title, colorbar=False, scalebar=scalebar, pre_rotated=True )
	elif target == 'unshared_PSD95':
		columns = {'Unshared PSD95':column}
		plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, \
			title=title, colorbar=False, scalebar=scalebar, pre_rotated=True )
	elif target == 'conc_CaMKII':
		columns = {'CaMKII':column}
		plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, \
			title=title, colorbar=False, scalebar=scalebar, pre_rotated=True )
	elif target == 'conc_unrotated_CaMKII':
		columns = {'CaMKII':column}
		plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, \
			title=title, colorbar=False, scalebar=scalebar, pre_rotated=False )
	elif target == 'conc_PSD95':
		columns = {'PSD95':column}
		plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, \
			title=title, colorbar=False, scalebar=scalebar, pre_rotated=True )
	elif target == 'conc_GluN2B':
		columns = {'GluN2B':column}
		plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, \
			title=title, colorbar=False, scalebar=scalebar, pre_rotated=True )
	elif target == 'conc_STG':
		columns = {'STG':column}
		plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, \
			title=title, colorbar=False, scalebar=scalebar, pre_rotated=True )
	elif target == 'watershed_CaMKII':
		columns = {'CaMKII':column}
		plot_watershed_region_from_a_direction(fig, num_rows, num_columns, row, columns, d, title=False, scalebar=scalebar )
	elif target == 'rdf':
		errorbar= 'shaded' # errorbar='shaded', 'line', or 'final_frame_alone'
		ax    = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
		plot_a_rdf( ax, d, errorbar=errorbar, legend=scalebar )
	elif target == 'rdf_CG':
		errorbar= 'shaded' # errorbar='shaded', 'line', or 'final_frame_alone'
		ax    = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
		plot_a_rdf( ax, d, errorbar=errorbar, legend=scalebar, target_molecules = ['All', 'CaMKII', 'GluN2B'] , ylim = (-0.01,1.0) )
	elif target == 'rdf_PSD95':
		errorbar= 'shaded' # errorbar='shaded', 'line', or 'final_frame_alone'
		ax    = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
		plot_a_rdf_PSD95( ax, d, legend=scalebar )
	elif target == 'concs_in_CaMKII':
		ax    = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
		t = 'CaMKII'
		plot_concs_condensate_bar(ax, t, d)
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
	elif target == 'concs_in_STG':
		ax    = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
		t = 'STG'
		plot_concs_condensate_bar(ax, t, d)
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
	elif 'energy' in target: # target == 'energy_anisotropic_STG'
		targets = target.split('_')
		type_energy       = targets[0] + '_' + targets[1]
		target_condensate = targets[-1]
		ax    = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
		ymax = 1.5
		value = plot_binding_energy_bar(ax, d, type_energy, target_condensate, ymax)
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
	else:
		raise ValueError("Target function has not been implemented: {}".format(target))

	return value



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



def flatten(sequence):
    result = []
    for item in sequence:
        if isinstance(item, (list, tuple, range, dict, set, frozenset)):
            result.extend(flatten(item))
        else:
            result.append(item)
    return result


def get_ratio_code(hex_code):
	hex_code   = hex_code.lstrip("#")
	ratio_code = [int(hex_code[i:i+2], 16)/255 for i in range(0, 6, 2)]
	return ratio_code


