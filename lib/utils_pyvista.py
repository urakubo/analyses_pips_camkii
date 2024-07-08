
import os, glob, pickle, pprint, copy
import numpy as np

import pyvista

import trimesh
from skimage import measure


import lib.utils as utils
import lib.parameters as p
import lib.colormap as c



def square_zx():
	x = p.space[0]/2
	z = p.space[2]/2
	pointa = [-x,  0.0, z]
	pointb = [-x, 0.0, -z]
	pointc = [x , 0.0, -z]
	pointd = [x , 0.0,  z]
	return pyvista.Rectangle([pointa, pointb, pointc])


def square_xy():
	x = p.space[0]/2
	y = p.space[1]/2
	pointa = [-x,  y, 0.0]
	pointb = [-x, -y, 0.0]
	pointc = [x , -y, 0.0]
	pointd = [x ,  y, 0.0]
	return pyvista.Rectangle([pointa, pointb, pointc])


def square_yz():
	y = p.space[1]/2
	z = p.space[2]/2
	pointa = [0.0, -y,  z]
	pointb = [0.0, -y, -z]
	pointc = [0.0, y , -z]
	pointd = [0.0, y ,  z]
	return pyvista.Rectangle([pointa, pointb, pointc])

def square_yz_mag(magnification):
	y = p.space[1]/2/magnification
	z = p.space[2]/2/magnification
	pointa = [0.0, -y,  z]
	pointb = [0.0, -y, -z]
	pointc = [0.0, y , -z]
	pointd = [0.0, y ,  z]
	return pyvista.Rectangle([pointa, pointb, pointc])




def plot_a_condensate_pyvista(d, pl, rotation=False): 
	
	flipz = False
	# Generate mesh
	r_CaMKII   = d['region_condensate_in_grid_mesh']['CaMKII'].astype(float)
	r_STG      = d['region_condensate_in_grid_mesh']['STG'].astype(float)
	
	if np.unique(r_CaMKII).shape[0] > 1:
		mesh_CaMKII = generate_mesh(r_CaMKII, flipz = flipz)
	if np.unique(r_STG).shape[0] > 1:
		mesh_STG    = generate_mesh(r_STG   , flipz = flipz)
	
	if rotation == True:
		utils.rotate(mesh_CaMKII, mesh_STG)
	
	# Add cube
	'''
	cube  = pyvista.Cube(center=(0,0,0), \
		x_length=utils.space[0], y_length=utils.space[1], z_length=utils.space[2])
	pl.add_mesh(cube, color='black', style='wireframe')
	'''
	if np.unique(r_CaMKII).shape[0] > 1:
		pl.add_mesh(mesh_CaMKII, color='green', show_edges=False,  opacity=0.4)
	if np.unique(r_STG).shape[0] > 1:
		pl.add_mesh(mesh_STG   , color='red', show_edges=False,  opacity=0.4)
	pl.set_background('white')
	
	
def plot_a_pre_rotated_condensate_pyvista(d, pl): 
	
	flipz = False
	# Generate mesh
	r_CaMKII   = d['rotated_region_condensate_in_grid_mesh']['CaMKII'].astype(float)
	if np.unique(r_CaMKII).shape[0] > 1:
		mesh_CaMKII = generate_mesh(r_CaMKII, flipz = flipz)
		pl.add_mesh(mesh_CaMKII, color='green', show_edges=False,  opacity=0.4)
	
	r_STG      = d['rotated_region_condensate_in_grid_mesh']['STG'].astype(float)
	if np.unique(r_STG).shape[0] > 1:
		mesh_STG    = generate_mesh(r_STG   , flipz = flipz)
		pl.add_mesh(mesh_STG   , color='red', show_edges=False,  opacity=0.4)
	
	pl.set_background('white')
	
	
def generate_mesh(volume, num_smoothing = 1, flipx = False, flipy = False, flipz = False):
	
	try:
		v_march, f_march, normals, values = measure.marching_cubes(volume, 0.5, spacing=(1,1,1), gradient_direction='ascent')
		center = np.array(p.space)/2
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
	
	except ValueError:
		print('Mesh was not generated.')
		return None
	
