
import os, glob, pickle, pprint, copy
import numpy as np

import utils
import pyvista

from skimage import measure
from skimage import morphology
import trimesh


def square_zx():
	x = utils.space[0]/2
	z = utils.space[2]/2
	pointa = [-x,  0.0, z]
	pointb = [-x, 0.0, -z]
	pointc = [x , 0.0, -z]
	pointd = [x , 0.0,  z]
	return pyvista.Rectangle([pointa, pointb, pointc, pointd])

def square_xy():
	x = utils.space[0]/2
	y = utils.space[1]/2
	pointa = [-x,  y, 0.0]
	pointb = [-x, -y, 0.0]
	pointc = [x , -y, 0.0]
	pointd = [x ,  y, 0.0]
	return pyvista.Rectangle([pointa, pointb, pointc, pointd])

def square_yz():
	y = utils.space[1]/2
	z = utils.space[2]/2
	pointa = [0.0, -y,  z]
	pointb = [0.0, -y, -z]
	pointc = [0.0, y , -z]
	pointd = [0.0, y ,  z]
	return pyvista.Rectangle([pointa, pointb, pointc, pointd])


def rorate(mesh_CaMKII, mesh_STG): 
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
	
	# From top-bottom to left-right
	#theta_xz = np.arctan2(-1, 0)
	#r3 = np.array([[np.cos(theta_xz), 0, -np.sin(theta_xz)],[0, 1, 0],[np.sin(theta_xz), 0, np.cos(theta_xz)]])
	#rot_matrix[:3,:3] = np.dot(r3, rot_matrix[:3,:3])
	#
	
	mm = trimesh.transformations.random_rotation_matrix()
	mesh_CaMKII.apply_transform(rot_matrix)
	mesh_STG.apply_transform(rot_matrix)

def generate_mesh(volume, num_smoothing = 1, flipx = False, flipy = False, flipz = False):
	v_march, f_march, normals, values = measure.marching_cubes(volume, 0.5, spacing=(1,1,1), gradient_direction='ascent')
	center = np.array(utils.space)/2
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
	
def plot_a_image(d, pl, rotation=True): 
	
	
	flipz = False
	# Generate mesh
	r_CaMKII   = d['region_condensate']['CaMKII'].astype(float)
	r_STG      = d['region_condensate']['STG'].astype(float)
	mesh_CaMKII = generate_mesh(r_CaMKII, flipz = flipz)
	mesh_STG    = generate_mesh(r_STG   , flipz = flipz)
	
	if rotation == True:
		rorate(mesh_CaMKII, mesh_STG)
	
	# Add cube
	'''
	cube  = pyvista.Cube(center=(0,0,0), \
		x_length=utils.space[0], y_length=utils.space[1], z_length=utils.space[2])
	pl.add_mesh(cube, color='black', style='wireframe')
	'''
	
	pl.add_mesh(mesh_CaMKII, color='green', show_edges=False,  opacity=0.4)
	pl.add_mesh(mesh_STG   , color='red', show_edges=False,  opacity=0.4)
	pl.set_background('white')
	
	
	
def save_a_plot(d, dir_img, prefix, suffix):
	pl = pyvista.Plotter(window_size=[300,800], shape=(3, 1), border=False)
	#pl.add_text( '{}_{}'.format(prefix, suffix), position='lower_left', color='k', font='arial', font_size=10)
	
	pl.subplot(0, 0)
	plot_a_image(d, pl, rotation=False)
	pl.add_mesh(square_yz(), color='black', style='wireframe')
	pl.view_yz()
	pl.camera.roll -= 90
	pl.camera.Zoom(1)
	
	pl.subplot(1, 0)
	plot_a_image(d, pl, rotation=False)
	pl.add_mesh(square_zx(), color='black', style='wireframe')
	pl.view_zx()
	pl.camera.roll -= 90
	pl.camera.Zoom(1)
	
	pl.subplot(2, 0)
	plot_a_image(d, pl, rotation=False)
	pl.add_mesh(square_xy(), color='black', style='wireframe')
	pl.view_xy()
	pl.camera.roll -= 90
	pl.camera.Zoom(1)
	
	#pl.show(interactive=True, auto_close=False) 
	pl.show(interactive=False, auto_close=True) # off_screen = True
	filename = os.path.join(dir_imgs, '{}_{}.png'.format(prefix, suffix))
	pl.screenshot(filename)
	
	
def save_plots_matrix(dir_data, dir_imgs, sigma): 
	STG    = [500, 1000, 2000, 3000, 4000]
	GluN2B = [500, 2000, 4000, 6000, 8000, 12000]
	suffix = 'sigma_{}'.format(sigma)
	
	pl = pyvista.Plotter(window_size=[1500,1500], shape=(len(GluN2B), len(STG)), border=False)
	
	for i, stg in enumerate(STG):
		for j, glun in enumerate(GluN2B):
			# Load data
			id = i + j * len(STG)
			prefix = str(id).zfill(3)
			d      = utils.load(dir_data, prefix, suffix)
			print('Target: {}, sigma: {}'.format(prefix, sigma))
			
			row    = i
			column = len(GluN2B)-j-1
			pl.subplot(column, row)
			
			plot_a_image(d, pl)
			pl.view_yz()
			pl.camera.roll -= 90
			
			#plot_a_image(d, pl, rotation=False)
			#pl.view_xy()
			#pl.camera.roll -= 90

			pl.camera.Zoom(1)
			
	pl.show(interactive=False, auto_close=True) # off_screen = True
	filename = os.path.join(dir_imgs, 'summary_{}.png'.format(suffix))
	pl.screenshot(filename)
	
	
if __name__ == '__main__':

	# Dataset 1
	filenames	= [	'PIPS',\
					'iPIPS',\
					'PartialE',\
					'Homo']	
	dir_data	= 'data'
	dir_imgs = 'imgs/region_condensates_pyvista'
	sigma    = 2 # 2, 3, or 4
	
	os.makedirs(dir_imgs, exist_ok=True)
	
	sigmas = [2,3,4]
	for sigma in sigmas:
		for filename in filenames:
			# Load data
			prefix = filename
			suffix = 'sigma_{}'.format(sigma)
			d      = utils.load(dir_data, prefix, suffix)	
			print('Target: {}, sigma: {}'.format(filename, sigma))
			save_a_plot(d, dir_imgs, prefix, suffix)
	
	
	'''
	# Dataset 2
	dir_data = 'data'
	dir_imgs = 'imgs/region_condensates_pyvista'
	sigma    = 2 # 2, 3, or 4
	save_plots_matrix(dir_data, dir_imgs, sigma)
	'''
	
	
	#r_CaMKII   = d['concs_in_grid_mesh']['CaMKII']
	#r_STG      = d['concs_in_grid_mesh']['STG']
	#pl.add_volume(r_CaMKII, clim=[0, 4.0], cmap="Greens", opacity='linear' )
	#pl.add_volume(r_STG   , clim=[0, 4.0], cmap="Reds"  , opacity='linear' ) # opacity='linear', "sigmoid"
	

