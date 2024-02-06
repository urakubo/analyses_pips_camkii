import os, glob, pickle, pprint, copy
import numpy as np
from skimage import morphology
from skimage import measure
from skimage.segmentation import watershed

from scipy import ndimage

import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1



from sklearn import mixture
#https://matsci.org/t/compatibility-issue-between-python-ovito-library-and-matplotlib/50794
os.environ['OVITO_GUI_MODE'] = '1'
from ovito.io import import_file
import pyvista
import trimesh
import colormap as c

plt.rcParams.update({
                    'pdf.fonttype' : 'truetype',
                    'svg.fonttype' : 'none',
                    'font.family' : 'sans-serif',
                    'font.sans-serif' : 'Arial',
                    'font.style' : 'normal'})

bsize     = 120
space     = [bsize, bsize, bsize]
space_np  = np.array(space)
center_np = space_np / 2

edge0 =  list(range(-int(space[0]/2), int(space[0]/2), 1))
edge1 =  list(range(-int(space[1]/2), int(space[1]/2), 1))
edge2 =  list(range(-int(space[2]/2), int(space[2]/2), 1))


reference_molecule_for_centering = 'All'

subunits = \
	{'GluN2Bc'	:{'id':3},\
	'CaMKII Hub'    :{'id':0},\
	'CaMKII Catalyst':{'id':5},\
	'STGc1'         :{'id':4},\
	'STGc2'         :{'id':2},\
	'PSD1'			:{'id':1}}


molecules_with_all = \
	{'CaMKII'	:{'s':['CaMKII Catalyst']	,'c':'#228B22'},\
	'GluN2B'	:{'s':['GluN2Bc']			,'c':'#ED0DD9'},\
	'STG'		:{'s':['STGc1','STGc2']		,'c':'r'},\
	'PSD95'		:{'s':['PSD1']				,'c':'#00FFFF'},\
	'All'		:{'s':['GluN2Bc','CaMKII Hub','CaMKII Catalyst','STGc1','STGc2','PSD1']			,'c':'k'}}

molecules_without_all = \
	{'CaMKII'	:{'s':['CaMKII Catalyst']	,'c':'#228B22'},\
	'GluN2B'	:{'s':['GluN2Bc']			,'c':'#ED0DD9'},\
	'STG'		:{'s':['STGc1','STGc2']		,'c':'r'},\
	'PSD95'		:{'s':['PSD1']				,'c':'#00FFFF'}}


for k, v in molecules_with_all.items():
	molecules_with_all[k]['id'] = [subunits[s]['id'] for s in v['s']]


for k, v in molecules_without_all.items():
	molecules_without_all[k]['id'] = [subunits[s]['id'] for s in v['s']]
	

# RDF parameters

rdf_bins          = np.arange(0, 35) # You can set "np.arange(0, 35, 2)"
rdf_num_sampling_frames = 5
rdf_sampling_interval   = 2


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
	
	
def load_data(dir_data, filename_data, id_frame):
	print('Load data.')
	data_all            = import_file(os.path.join(dir_data, filename_data), input_format= "lammps/dump" )
	data_target_frame   = data_all.compute(id_frame)
	types, positions, id_molecule = decode_data(data_target_frame)
	return types, positions, id_molecule
	
	
def load_lammpstrj(dir_data, filename_data, id_frame):
	print('Load data.')
	data_all            = import_file(os.path.join(dir_data, filename_data), input_format= "lammps/dump" )
	data_target_frame   = data_all.compute(id_frame)
	time_stamp 			= data_all.compute(id_frame).attributes['Timestep']
	types, positions, id_molecule = decode_data(data_target_frame)

	return types, positions, id_molecule, time_stamp
	
	
def decode_data(data_frame):
	type         = np.array( data_frame.particles['Particle Type'] )
	position     = np.array( data_frame.particles['Position'] )
	id_molecule  = np.array( data_frame.particles['Molecule Identifier'] )
	# data_target_frame.particles['bp']
	return type, position, id_molecule
	
def decode_species(types, positions):
	types_binary = {k: [True if t in v['id'] else False for t in types] for k, v in molecules_with_all.items() }
	types_positions = {k: positions[types_binary[k],:] for k in molecules_with_all.keys() }
	return types_positions
	
	
	
############### Edit data
	
	## It provides the same operation but offers shorter time: 
	## ids_loc = np.nonzero( get_hist(locs_subunits) )
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
	
	
	#
	# grid_coord: Molecular locations in the grid space [0, 1,..., 119]
	# (numpy uint) [(x0,y0,z0), (x1,y1,z1), ..., (xn,yn,zn)]
	#
	# real_coord: Molecular locations in the real space [-60, 60) 
	# (numpy float) [(x0,y0,z0), (x1,y1,z1), ..., (xn,yn,zn)]
	#
	# grid_mesh : molecular existance in a 3D numpy variable
	# (True/False in the 3D space (util.space[0], util.space[1], util.space[2]))
	#
	
def get_concs_and_condensates(types, positions, ids_molecule, sigma=2):
	
	# Parameters
	targs_molecule  = molecules_with_all.keys() # ['GluN2B', 'CaMKII', 'STG', 'PSD95', 'All']
	
	# Get the locations of target molecules separately.
	flag_type          = {t: [True if k in molecules_with_all[t]['id'] else False for k in types] for t in targs_molecule}
	locs_in_real_coord = {k: positions[flag_type[k],:] for k in targs_molecule } 
	locs_in_grid_mesh  = {k: get_hist(locs_in_real_coord[k]) for k in targs_molecule}
	
	# Get a peripheral region (a region outside of a sphere) 
	# and obtain the concentrations in this area (diluted region concentration for a baseline).
	region_periphery = get_periphery_in_grid_mesh()
	concs_periphery  = {t: np.sum(locs_in_grid_mesh[t] * region_periphery ) / np.sum( region_periphery ) for t in targs_molecule}
	
	# Get the condenate regions of targs_molecule.
	concs_in_grid_mesh = {t: ndimage.gaussian_filter(locs_in_grid_mesh[t], sigma = sigma) for t in targs_molecule}
	regions_condensate_in_grid_mesh = {t: get_high(concs_in_grid_mesh[t]-concs_periphery[t]) for t in targs_molecule}
	
	
	def get_concs_condensate(ref_molecule):
		return {t: np.sum(locs_in_grid_mesh[t] * regions_condensate_in_grid_mesh[ref_molecule])/ \
				np.sum( regions_condensate_in_grid_mesh[ref_molecule] ) for t in targs_molecule }
	concs_condensate = {t: get_concs_condensate(t) for t in targs_molecule}
	
	# The line below was prepared for visualization. So far, it does not contribute to the improvement.
	regions_condensate_in_grid_mesh_0_15 = {t: get_high(concs_in_grid_mesh[t]-concs_periphery[t], th=0.15) for t in targs_molecule}
	
	
	# Summary
	d = {
			'locs_in_grid_mesh':	locs_in_grid_mesh,
			'concs_in_grid_mesh':	concs_in_grid_mesh,
			'region_condensate_in_grid_mesh': regions_condensate_in_grid_mesh,
			'region_condensate_in_grid_mesh_0_15': regions_condensate_in_grid_mesh_0_15,
			'conc_periphery'	:	concs_periphery,
			'conc_condensate'	:	concs_condensate,
		}
	return d
	
	
############### Watershed segmentation
	
def watershed_segmentation( d ):
	markers = np.zeros(space, dtype = int)
	loc = int(space[0]/2), int(space[1]/2), int(space[2]/2)
	markers[ loc ]    = 1
	markers[10,10,10] = 2
	
	tot_volume = space[0] * space[1] * space[2]
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

def get_center_of_mass(types_, positions_, reference_molecule_for_centering = reference_molecule_for_centering):
	
	types = [True if t in molecules_with_all[reference_molecule_for_centering]['id'] else False for t in types_ ]
	position_ref = positions_[types,:]
	
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

	center_of_mass = np.array( [theta_x * space[0], theta_y * space[1] , theta_z * space[2] ] ) + np.pi
	center_of_mass /= (2 * np.pi)

	return center_of_mass

	# Reference
	# https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions


def centering(p, center): 
	p_centered = p - center
	for dim in [0,1,2]:
		over  = (p_centered[:,dim]  >  space[dim]/2)
		under = (p_centered[:,dim] <= -space[dim]/2)
		p_centered[over ,dim] -= space[dim]
		p_centered[under,dim] += space[dim]
	return p_centered


def get_average_center_of_mass( dir_lammpstrj, filename_input, sampling_time_frames ):
	centers = np.zeros( (len(sampling_time_frames), 3) )
	
	for i, target_frame in enumerate( sampling_time_frames ):
		types, positions_grid_coord,ids_molecule = load_data( dir_lammpstrj, filename_input, target_frame )
		positions_real_coord = centering(positions_grid_coord, center_np)
		centers[i,:] = get_center_of_mass(types, positions_real_coord) - center_np
	
	center = np.mean(centers, axis = 0)
	center -= space_np * (center >  center_np)
	center += space_np * (center <= -center_np)
	
	return center


############### Radial distrbution function

def get_lattice_grids():
	x, y, z = np.meshgrid(edge0, edge1, edge2)
	x, y, z = np.ravel(x), np.ravel(y), np.ravel(z)
	grids   = np.vstack((x,y,z)).T
	return grids


def get_a_rdf(types, positions, rdf_grid_points, rdf_bins, center = None):
	
	# Centering
	if center is None:
		center              = get_center_of_mass(types, positions)
	
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


def get_rdfs_from_multiple_frames( dir_lammpstrj, filename_input, sampling_time_frames, rdf_bins, rdf_grid_points, center = None ):

	if center is None:
		get_average_center_of_mass( dir_lammpstrj, filename_input, sampling_time_frames )
		
	rdfs = { k: np.zeros( ( len(rdf_bins)-1, len(sampling_time_frames) ) ) for k in molecules_with_all.keys() }
	for i, id_frame in enumerate( sampling_time_frames ):
		types, positions, _ = load_data( dir_lammpstrj, filename_input, id_frame )
		current_rdfs = get_a_rdf(types, positions, rdf_grid_points, rdf_bins, center )
		for k in rdfs.keys():
			rdfs[k][:,i] = current_rdfs[k]
	return rdfs

	
def get_rdfs( dir_input, filename_input, target_frame, center=None ):
	
	# Parameters
	rdf_grid_points   = get_lattice_grids()
	
	# Target frames
	rdf_sampling_frames = list(range( target_frame - (rdf_num_sampling_frames - 1)*rdf_sampling_interval, \
			target_frame + rdf_sampling_interval,\
			rdf_sampling_interval))
	
	if np.any(np.array(rdf_sampling_frames) < 0):
		rdf_sampling_frames = [0]

	print('rdf_sampling_frames ', rdf_sampling_frames)
	rdf = get_rdfs_from_multiple_frames(dir_input, filename_input, \
			rdf_sampling_frames, rdf_bins, rdf_grid_points, center = center)
	
	return rdf, rdf_bins, rdf_sampling_frames


def plot_a_rdf( ax, d, errorbar='shaded', legend=True, target_molecules = molecules_without_all.keys() , ylim = (-0.006,0.66) ):
	r = d['rdf_bins'][1:-1]
	for k in target_molecules:
		rdf_mean  = np.mean( d['rdf'][k][1:], axis = 1 )
		rdf_std   = np.std(  d['rdf'][k][1:], axis = 1 )
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

	ax.set_xlabel('Distance from \n center-of-mass (l.u.)')
	ax.set_ylabel('(beads / volume)')
	ax.set_xlim(0,30)
	ax.set_ylim(*ylim)
	ax.set_xticks(range(0,40,10))
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
	return cb
	
	
def plot_scalebar(ax, col='k', linewidth=2):
	ax.plot([5,25],[5,5], '-', color=col, linewidth=linewidth)
	return
	
	
def make_a_panel_of_CaMKII_STG_condenstates(d, transp, slice):
	# Condensate region
	r_CaMKII   = d['region_condensate_in_grid_mesh']['CaMKII'].transpose(transp)[slice,:,:]
	r_STG      = d['region_condensate_in_grid_mesh']['STG'].transpose(transp)[slice,:,:]
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
	
	
	
def make_a_panel_of_CaMKII_All_condenstates(d, transp, slice):
	# Condensate region
	r_CaMKII   = d['region_condensate_in_grid_mesh']['CaMKII'].transpose(transp)[slice,:,:]
	r_STG      = d['region_condensate_in_grid_mesh']['All'].transpose(transp)[slice,:,:]
	r_BOTH     = r_CaMKII & r_STG
	# white panel
	panel = np.ones( [r_CaMKII.shape[0],r_CaMKII.shape[1],3], dtype=np.uint8 )*255 
	col = c.cmap_universal_uint['All']
	for k in range(3): panel[r_STG,k] = col[k]
	col = c.cmap_universal_uint['CaMKII']
	for k in range(3): panel[r_CaMKII,k] = col[k]
	col = c.light_green_universal_uint # Lignt green
	for k in range(3): panel[r_BOTH,k] = col[k]
	return panel
	
	
def plot_regions_condenstate_from_a_direction(fig, num_rows, num_columns, row, column, d, transp = (0,1,2), title=True, scalebar=True):
	slice = int( space[0]/2 )
	panel = make_a_panel_of_CaMKII_All_condenstates(d, transp, slice)
	ax    = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
	ax.imshow( panel )
	if scalebar==True:
		plot_scalebar(ax)
	if title == True:
		ax.set_title('Green: CaMKII, \n Gray: All')
	arrange_graph_no_ticks(ax)
	return ax
	
	
def plot_regions_condenstate_from_a_direction_(fig, num_rows, num_columns, row, column, d, transp = (0,1,2), title=True, scalebar=True):
	slice = int( space[0]/2 )
	panel = make_a_panel_of_CaMKII_STG_condenstates(d, transp, slice)
	ax    = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
	ax.imshow( panel )
	if scalebar==True:
		plot_scalebar(ax)
	if title == True:
		ax.set_title('Green: CaMKII, \n Red: STG')
	arrange_graph_no_ticks(ax)
	return ax
	
	
def plot_concs_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp = (0,1,2), title=True, colorbar=True, scalebar=True, vmin=0, vmax=None ):
	slice = int( space[0]/2 )
	axes = []
	for i, (target, column) in enumerate(columns.items()):
		panel = d['concs_in_grid_mesh'][target].transpose(transp)[slice,:,:]
		ax = fig.add_subplot( num_rows, num_columns, row*num_columns+column )
		axes.append(ax)
		cs = ax.imshow( panel , cmap=c.cmap[target], vmin=vmin, vmax=vmax )
		if i == 0 and scalebar==True:
			plot_scalebar(ax, col='w', linewidth=3)
		if title == True:
			ax.set_title('Smoothed '+ target )
		arrange_graph_no_ticks(ax)
		if colorbar == True:
			cb = plot_colorbar(ax, cs)
			cb.set_label('(beads / voxel)')
	return axes
	
	
def plot_watershed_region_from_a_direction(fig, num_rows, num_columns, row, columns, d, transp = (0,1,2), title=True, scalebar=True ):
	slice = int( space[0]/2 )
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
	
	
def plot_concs_condensate_bar(ax, targ, ref, d):
	cc = d['conc_condensate']
	col = c.cmap_universal_ratio[targ]
	ax.bar(['In '+targ+'\n condensates', 'In '+ref+'\n condensates'], [cc[targ][targ], cc[ref][targ]], width=0.5, color=col)
	ax.set_title('Conc of {}'.format(targ))
	ax.set_ylabel('(beads / volume)')
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
	

############### 3D plot using pyvista


def square_zx():
	x = space[0]/2
	z = space[2]/2
	pointa = [-x,  0.0, z]
	pointb = [-x, 0.0, -z]
	pointc = [x , 0.0, -z]
	pointd = [x , 0.0,  z]
	return pyvista.Rectangle([pointa, pointb, pointc, pointd])


def square_xy():
	x = space[0]/2
	y = space[1]/2
	pointa = [-x,  y, 0.0]
	pointb = [-x, -y, 0.0]
	pointc = [x , -y, 0.0]
	pointd = [x ,  y, 0.0]
	return pyvista.Rectangle([pointa, pointb, pointc, pointd])


def square_yz():
	y = space[1]/2
	z = space[2]/2
	pointa = [0.0, -y,  z]
	pointb = [0.0, -y, -z]
	pointc = [0.0, y , -z]
	pointd = [0.0, y ,  z]
	return pyvista.Rectangle([pointa, pointb, pointc, pointd])


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
	
	
def generate_mesh(volume, num_smoothing = 1, flipx = False, flipy = False, flipz = False):
	v_march, f_march, normals, values = measure.marching_cubes(volume, 0.5, spacing=(1,1,1), gradient_direction='ascent')
	center = np.array(space)/2
	v_march = v_march - center
	
	if flipx == True:
		v_march[:,0] = -v_march[:,0]
	if flipy == True:
		v_march[:,1] = -v_march[:,1]
	if flipz == True:
		v_march[:,2] = -v_march[:,2]
	
	mesh = trimesh.Trimesh(vertices=v_march, faces=f_march)
	mesh = trimesh.smoothing.filter_humphrey(mesh, alpha = 1.0, beta=0.0, iterations=num_smoothing)
	return mesh
	
	
	
