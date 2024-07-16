
import numpy as np
import networkx as nx

from ovito.pipeline import ModifierInterface, Pipeline, StaticSource
from ovito.data import DataCollection, ParticleType

import lib.utils as utils
import lib.parameters as p
import lib.utils_graph as utils_graph


def _get_beads_in_largest_cluster(data_all, time_frame):
	data   = data_all.compute(time_frame)
	types, positions, ids_molecule = utils.decode_data(data)
	bp 	= np.array( data.particles['bp'] ).astype('int')
	multi_graph = utils_graph.get_multi_graph(ids_molecule, types, bp, positions)
	clusters = sorted(nx.connected_components(multi_graph), key=len, reverse=True)
	g_largest_cluster = nx.MultiGraph( multi_graph.subgraph(clusters[0]) )
	ids_GluN2B = [v['ids_bead'][0].tolist() for n, v in g_largest_cluster.nodes.items() if v['species'] == 'GluN2B' ]
	ids_CaMKII = [v['ids_bead'][0].tolist() for n, v in g_largest_cluster.nodes.items() if v['species'] == 'CaMKII' ]
	ids_PSD95  = [v['ids_bead'][0].tolist() for n, v in g_largest_cluster.nodes.items() if v['species'] == 'PSD95' ]
	ids_STG    = [v['ids_bead'][0].tolist() for n, v in g_largest_cluster.nodes.items() if v['species'] == 'STG' ]
	ids_GluN2B = utils.flatten(ids_GluN2B)
	ids_CaMKII = utils.flatten(ids_CaMKII)
	ids_PSD95  = utils.flatten(ids_PSD95 )
	ids_STG    = utils.flatten(ids_STG   )
	return ids_GluN2B, ids_CaMKII, ids_PSD95, ids_STG


def get_beads_in_largest_cluster(data_all, time_frame):
	ids_GluN2B, ids_CaMKII, ids_PSD95, ids_STG = _get_beads_in_largest_cluster(data_all, time_frame)
	ids_bead_cluster = ids_GluN2B + ids_CaMKII + ids_PSD95 + ids_STG
	return ids_bead_cluster
	
def get_CaMKII_beads_in_largest_cluster(data_all, time_frame):
	ids_GluN2B, ids_CaMKII, ids_PSD95, ids_STG = _get_beads_in_largest_cluster(data_all, time_frame)
	return ids_CaMKII
	
def get_GluN2B_beads_in_largest_cluster(data_all, time_frame):
	ids_GluN2B, ids_CaMKII, ids_PSD95, ids_STG = _get_beads_in_largest_cluster(data_all, time_frame)
	return ids_GluN2B
	
def get_PSD95_beads_in_largest_cluster(data_all, time_frame):
	ids_GluN2B, ids_CaMKII, ids_PSD95, ids_STG = _get_beads_in_largest_cluster(data_all, time_frame)
	return ids_PSD95

def get_STG_beads_in_largest_cluster(data_all, time_frame):
	ids_GluN2B, ids_CaMKII, ids_PSD95, ids_STG = _get_beads_in_largest_cluster(data_all, time_frame)
	return ids_STG




class CenteringEachFrameModifier(ModifierInterface):
	def modify(self, data, frame, **kwargs):
		positions = data.particles.positions
		center    = utils.get_center_of_mass_simple(positions)
		data.particles_.positions_ = utils.centering_lattice_space(positions, center)
		
		
class CenteringModifier(ModifierInterface):
	def modify(self, data, frame, **kwargs):
		positions = data.particles.positions
		data.particles_.positions_ = utils.centering_lattice_space(positions, self.center)
	
	'''
	# How to use
	num_frames = data_all.source.num_frames
	types, positions, ids_molecule = utils.decode_data(data_all.compute(num_frames))
	center = utils.get_center_of_mass(types, positions)
	for dim in [0,1,2]:
		center[dim] += - p.space[dim] * (center[dim]  >=  p.space[dim]) + p.space[dim] * (center[dim]  <  0)
	'''
	
	
def delete_particle(frame, data, ids):
	transparency = np.ones((data.particles.count), dtype='int')
	transparency[ids] = 0

class DeleteParticle(ModifierInterface):
	def modify(self, data, frame, **kwargs):
		transparency = np.ones((data.particles.count), dtype='int')
		transparency[self.ids] = 0
		data.particles_.delete_elements(transparency)
		#print('transparency ', transparency)

class MakeTransparent(ModifierInterface):
	def modify(self, data, frame, **kwargs):
		transparency = np.ones((data.particles.count), dtype='float')*0.95
		transparency[self.ids] = 0.0
		data.particles_.create_property('Transparency', data=transparency)

class SetColorsById(ModifierInterface):
	def modify(self, data, frame, **kwargs):
		color_values = np.ones((data.particles.count, 3), dtype='float')
		for id, col in self.ids_colors.items():
			color_values[id,:] = col
		data.particles_.create_property('Color', data=color_values)

def compute_particle_transparency(frame, data):
    transparency = np.ones((data.particles.count))*0.3
    data.particles_.create_property('Transparency', data=transparency)


class FixBoundingbox(ModifierInterface):
	def modify(self, data, frame, **kwargs):
		# Bounding box
		# https://ovito.org/manual/python/introduction/examples/modifiers/shrink_wrap_box.html
		#   (x_max-x_min  0            0            x_min)
		#   (0            y_max-y_min  0            y_min)
		#   (0            0            z_max-z_min  z_min)
		coords_min = np.array([0,0,0])
		matrix = np.empty((3,4))
		matrix[:,:3] = np.diag( p.space_np - coords_min )
		matrix[:, 3] = coords_min

		# Assign the cell matrix - or create whole new SimulationCell object in
		# the DataCollection if there isn't one already.
		data.create_cell(matrix, (False, False, False))


def add_sphere_yellow(center, rm):
	
	print('center: ', center)
	#loc_ = (120,120,120) #12_002 colony 2
	#loc_ = (0,0,120) #12_006 colony 3
	data_photobleach = DataCollection()
	particles = data_photobleach.create_particles()
	#pos = [(center[0]+120, center[1]+120+rm, center[2]+120)]
	for dim in [0,1,2]:
		if center[dim]  >=  p.space[dim]:
			center[dim] -= p.space[dim]
		elif center[dim]  <  0:
			center[dim] += p.space[dim]
	print('center: ', center)
	pos = [(center[0], center[1]+rm, center[2])]
	pos_prop  = particles.create_property('Position', data=pos)
	type_prop = particles.create_property('Particle Type')
	type_prop.types.append(ParticleType(id = 1, name = 'Photobleach', color = (1.0,1.0,0.0), radius=rm))
	type_prop[0] = 1
	particles.create_property('Transparency', data=0.6)
	pipeline = Pipeline(source = StaticSource(data = data_photobleach))
	return pipeline
	

'''
line_vis = data_all.source.data.line.vis
line_vis.color = (0.0, 0.0, 0.0)
line_vis.width = 0.5
positions = np.array([[0,0,60],[0,120,60],[120,120,60],[120,0,60],[0,0,60] ])
lines_data = ovito.data.Lines
lines_data.create_line = positions
lines_data.color = (0.0, 0.0, 0.0)
lines_data.width = 0.5

line = GenerateTrajectoryLinesModifier(only_selected=False)
# line.generate()
line.vis.color = (0.0, 0.0, 0.0)
line.vis.width = 0.5
data_all.modifiers.append(line)
'''

	