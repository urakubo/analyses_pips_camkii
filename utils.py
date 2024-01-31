import os, glob, pickle, pprint, copy
import numpy as np
from skimage import morphology

from scipy import ndimage

import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1

from sklearn import mixture
from ovito.io import import_file

import colormap as c


plt.rcParams.update({
                    'pdf.fonttype' : 42,
                    'font.family' : 'sans-serif',
                    'font.sans-serif' : 'Arial',
                    'font.style' : 'normal'})

bsize    = 120
space    = [bsize, bsize, bsize]

edge0 =  list(range(-int(space[0]/2), int(space[0]/2), 1))
edge1 =  list(range(-int(space[1]/2), int(space[1]/2), 1))
edge2 =  list(range(-int(space[2]/2), int(space[2]/2), 1))


subunits = \
	{'GluN2Bc'	:{'id':3},\
	'CaMKII Hub'    :{'id':0},\
	'CaMKII Catalyst':{'id':5},\
	'STGc1'         :{'id':4},\
	'STGc2'         :{'id':2},\
	'PSD1'			:{'id':1}}


molecules_with_all = \
	{'GluN2B'	:{'s':['GluN2Bc']			,'c':'#ED0DD9'},\
	'CaMKII'	:{'s':['CaMKII Catalyst']	,'c':'#228B22'},\
	'STG'		:{'s':['STGc1','STGc2']		,'c':'r'},\
	'PSD95'		:{'s':['PSD1']				,'c':'#00FFFF'},\
	'All'		:{'s':['GluN2Bc','CaMKII Hub','CaMKII Catalyst','STGc1','STGc2','PSD1']			,'c':'k'}}

molecules_without_all = \
	{'GluN2B'	:{'s':['GluN2Bc']			,'c':'#ED0DD9'},\
	'CaMKII'	:{'s':['CaMKII Catalyst']	,'c':'#228B22'},\
	'STG'		:{'s':['STGc1','STGc2']		,'c':'r'},\
	'PSD95'		:{'s':['PSD1']				,'c':'#00FFFF'}}


for k, v in molecules_with_all.items():
	molecules_with_all[k]['id'] = [subunits[s]['id'] for s in v['s']]


for k, v in molecules_without_all.items():
	molecules_without_all[k]['id'] = [subunits[s]['id'] for s in v['s']]


	## It provides the same operation but offers shorter time: 
	## ids_loc = np.nonzero( utils.get_hist(locs_subunits) )
def get_locs_in_grid_coord(locs):
	locs_int = np.floor(locs).astype(int)
	loc0 = np.array( [edge0.index(i) for i in locs_int[:,0]] )
	loc1 = np.array( [edge1.index(i) for i in locs_int[:,1]] )
	loc2 = np.array( [edge2.index(i) for i in locs_int[:,2]] )
	locs_in_grid_space = np.vstack([loc0, loc1, loc2]).T
	return locs_in_grid_space


def get_periphery_in_grid_mesh():
	radius = (np.min(space) - 1) / 2
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
	
	
def load_data(dir_data, filename_data, id_frame):
	print('Load data.')
	data_all            = import_file(os.path.join(dir_data, filename_data), input_format= "lammps/dump" )
	data_target_frame   = data_all.compute(id_frame)
	types, positions, id_molecule = decode_data(data_target_frame)
	return types, positions, id_molecule
	
	
def decode_data(data_frame):
	type         = np.array( data_frame.particles['Particle Type'] )
	position     = np.array( data_frame.particles['Position'] )
	id_molecule  = np.array( data_frame.particles['Molecule Identifier'] )
	# data_target_frame.particles['bp']
	return type, position, id_molecule
	
	
def get_high(conc, th = 0.5):
	return (conc > np.max(conc)*th)

	
def get_hist(positions):
	edges0 =  list(range(-int(space[0]/2), int(space[0]/2) + 1, 1))
	edges1 =  list(range(-int(space[1]/2), int(space[1]/2) + 1, 1))
	edges2 =  list(range(-int(space[2]/2), int(space[2]/2) + 1, 1))
	H, (xedges, yedges, zedges) = np.histogramdd(positions, bins=(edges0, edges1, edges2))
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
	
	
def obtain_center_of_mass(position_ref):
	
	# print('position_ref: ', position_ref)
	x = (position_ref[:,0] / space[0]) * 2 * np.pi
	y = (position_ref[:,1] / space[1]) * 2 * np.pi
	z = (position_ref[:,2] / space[2]) * 2 * np.pi
	x0, x1 = np.cos(x), np.sin(x)
	y0, y1 = np.cos(y), np.sin(y)
	z0, z1 = np.cos(z), np.sin(z)

	x0, x1 = np.mean(x0), np.mean(x1)
	y0, y1 = np.mean(y0), np.mean(y1)
	z0, z1 = np.mean(z0), np.mean(z1)

	theta_x = np.arctan2(x1, x0)
	theta_y = np.arctan2(y1, y0)
	theta_z = np.arctan2(z1, z0)
	#theta_x = np.arctan2(-x1, -x0)
	#theta_y = np.arctan2(-y1, -y0)
	#theta_z = np.arctan2(-z1, -z0)

	center_of_mass = np.array( [theta_x * space[0], theta_y * space[1] , theta_z * space[2] ] ) + np.pi
	center_of_mass /= (2 * np.pi)
	print('Center_of_mass: ', center_of_mass )
	return center_of_mass

	# Reference
	# https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions


def get_center(types_, positions_, reference_molecule_for_centering):
	types = {}
	for k in molecules_with_all.keys():
		types[k] = [True if t in molecules_with_all[k]['id'] else False  for t in types_ ]
	p_for_centering = positions_[types[reference_molecule_for_centering],:]
	center     = obtain_center_of_mass(p_for_centering)
	return center



def centering(p, center): 
	p_centered = p - center
	for dim in [0,1,2]:
		over  = (p_centered[:,dim]  >  space[dim]/2)
		under = (p_centered[:,dim] <= -space[dim]/2)
		p_centered[over ,dim] -= space[dim]
		p_centered[under,dim] += space[dim]
	return p_centered


def decode_species(types, positions):
	types_binary = {k: [True if t in v['id'] else False for t in types] for k, v in molecules_with_all.items() }
	types_positions = {k: positions[types_binary[k],:] for k in molecules_with_all.keys() }
	return types_positions


############### Radial distrbution function

def get_lattice_grids():
	x, y, z = np.meshgrid(edge0, edge1, edge2)
	x, y, z = np.ravel(x), np.ravel(y), np.ravel(z)
	grids   = np.vstack((x,y,z)).T
	return grids


def get_rdf(types, positions, rdf_grid_points, rdf_bins, reference_molecule_for_centering):
	
	# Centering
	center              = get_center(types, positions, reference_molecule_for_centering)
	positions           = centering(positions, center)
	positions_grid_centered = centering(rdf_grid_points, center)
	
	# Decode species
	types_positions     = decode_species(types, positions)
	
	# Get distances from the center
	dists = {k: np.linalg.norm(v, axis=1) for k, v in types_positions.items()}
	dists_grid = np.linalg.norm(positions_grid_centered, axis=1)

	# Get radial distribution function (rdf)
	num_grid_around_center, _  = np.histogram(dists_grid , bins=rdf_bins)
	rdf = {}
	for k, v in dists.items():
		num_molecule_around_center, _  = np.histogram(v , bins=rdf_bins)
		rdf[k] = num_molecule_around_center / num_grid_around_center
	
	return rdf



############### Profiles / Each panel plot

def arrange_graph_no_ticks(ax):
	ax.xaxis.set_tick_params(labelbottom=False)
	ax.yaxis.set_tick_params(labelleft=False)
	ax.set_xticks([])
	ax.set_yticks([])
	
	
def show_colorbar(ax, cs):
	divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
	cax = divider.append_axes('right', '5%', pad='3%')
	cb = plt.colorbar(cs, cax=cax)
	ticks = copy.copy( cb.get_ticks() ).tolist()
	cb.set_ticks(ticks)
	cb.set_ticklabels(["{:.2f}".format(i) for i in ticks])
	return cb
	
	
def make_a_panel_of_CaMKII_STG_condenstates(d, transp, slice):
	# Condensate region
	r_CaMKII   = d['region_condensate']['CaMKII'].transpose(transp)[slice,:,:]
	r_STG      = d['region_condensate']['STG'].transpose(transp)[slice,:,:]
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
	
	
def plot_regions_condenstate_from_a_direction(fig, num_rows, num_columns, row, column, d, transp = (0,1,2), title=True):
	slice = int( space[0]/2 )
	panel = make_a_panel_of_CaMKII_STG_condenstates(d, transp, slice)
	ax    = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
	ax.imshow( panel )
	if title == True:
		ax.set_title('Green: CaMKII, \n Red: STG')
	arrange_graph_no_ticks(ax)
	
	
def plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp = (0,1,2), title=True, colorbar=True ):
	slice = int( space[0]/2 )
	for target, column in columns.items():
		panel = d['concs_in_grid_mesh'][target].transpose(transp)[slice,:,:]
		ax = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
		cs = ax.imshow( panel , cmap=c.cmap[target] )
		if title == True:
			ax.set_title('Smoothed '+ target )
		arrange_graph_no_ticks(ax)
		if colorbar == True:
			cb = show_colorbar(ax, cs)
			cb.set_label('(Beads / voxel)')
	
	
def plot_watershed_region_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp = (0,1,2), title=True ):
	slice = int( space[0]/2 )
	for target, column in columns.items():
		panel = d['labels_watershed_in_grid_mesh'][target].transpose(transp)[slice,:,:]
		ax = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
		cs = ax.imshow( panel, cmap="binary", vmin = 0, vmax = 1.5 )
		if title == True:
			ax.set_title('Watershed \n separated by '+ target )
		arrange_graph_no_ticks(ax)


