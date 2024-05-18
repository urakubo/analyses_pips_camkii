
import numpy as np
from ovito.pipeline import ModifierInterface
from ovito.data import DataCollection
#from ovito.pipeline import *
#from ovito.qt_compat import QtCore
import lib.utils as utils
import lib.parameters as p


class CenteringModifier(ModifierInterface):
	def modify(self, data, frame, **kwargs):
		positions = data.particles.positions
		data.particles_.positions_ = utils.centering_lattice_space(positions, self.center)
		
		
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

	