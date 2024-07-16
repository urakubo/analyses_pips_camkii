
# https://yasakura.me/post/classification-and-regression-with-decision-tree/
# https://tech.nkhn37.net/scikit-learn-decision-tree-random-forest/
# https://scikit-learn.org/stable/auto_examples/tree/plot_unveil_tree_structure.html#sphx-glr-auto-examples-tree-plot-unveil-tree-structure-py

import os, glob, pickle, pprint, copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import lib.utils as utils
import lib.utils_fitting as utils_fitting

import lib.parameters as p
import lib.colormap as c


from skimage.measure import label
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import r2_score

from scipy.optimize import curve_fit
import warnings



from matplotlib.colors import ListedColormap
from sklearn.tree import DecisionTreeClassifier, export_graphviz
import pydotplus
import matplotlib.pyplot as plt
from pydotplus import graph_from_dot_data

from bin.conc_dependence1_calc_connectivity_plot_phase_diagram import \
	HandleConnectivityPhaseDiagramConcDependence, \
	HandleCondVolumePhaseDiagramConcDependence, \
	PlotPhaseDiagramConcDependence

from specification_datasets import SpecDatasets



plt.rcParams.update( p.rc_param )


	
def prepare_plot():
	# fig  = plt.figure(figsize=(5, 5))
	fig  = plt.figure(figsize=(4.22, 4.22))
	fig.subplots_adjust(left = 0.20, bottom=0.2)
	ax = fig.add_subplot( 1, 1, 1 )
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_box_aspect(1)
	return fig, ax
	
	
def export_tree_graph(tree):
	
	dot_data = export_graphviz(tree,
	  filled=True,
	  rounded=True,
	  class_names=['Setosa', 'Versicolor', 'Virginica'],
	  feature_names=['petal length','petal width'],
	  out_file=None)
	graph = graph_from_dot_data(dot_data)
	graph.progs = {'dot': u"C:\\Program Files\\Graphviz\\bin\\dot.exe"}
	graph.write_png('tree.png') 
	
	
def set_colors( data, target_type ):
	
	## Set Color
	m = 255
	gray = [230/m,230/m,230/m]
	cols = []
	for two_phase, partial_engulfment, phase_diagram in \
		zip( data['two_phase_condensate'], data['partial_engulfment'], data['phase_diagram'] ):
		if target_type == 'phase_diagram':
			if (two_phase == 1) and (partial_engulfment == 0):
				cols.append([255/m,191/m/2,192/m/2]) # 'r'
			elif (two_phase == 1) and (partial_engulfment == 1):
				cols.append('w')
			elif (two_phase == 0) and (phase_diagram == 3):
				cols.append([77/m,196/m, 255/m]) # 'b'
			else:
				cols.append(gray)
		elif target_type == 'two_phase':
			cols.append('k' if two_phase == 1 else 'w')
		elif target_type == 'pips_partial':
			cols.append('k' if partial_engulfment == 0 else 'w')
		else:
			sys.exit('Unknown target type: {}'.format(target_type) )
			
	return cols
	
	
	
class Fit2DPhaseDiagramConcDependence():
	
	'''
	type_color: 'phase_diagram', 'two_phase_condensate', 'pips_partial'
	'''
	
	def __init__( self ):
		
		
		## Load data
		variable_names = [ 	'Number of GluN2B bound to one CaMKII',\
							'Number of STG bound to one PSD-95' ,
							'Raio of PSD95 bound to both GluN2B and STG',
							'Number of GluN2B bound to one PSD-95' ]
		species        = [ 'CaMKII', 'PSD95' , 'PSD95' , 'PSD95']
		type_analyses  = [ 'average', 'average_STG', 'ratio', 'average_GluN2B']
		
		data = {}
		for s, t, variable_name in zip( species, type_analyses, variable_names ):
			pl = HandleConnectivityPhaseDiagramConcDependence(s, t)
			pl.conc_dependence_merged()
			pl.load_data()
			data[variable_name] = pl.data[:-1,:].flatten()
			# print(variable_name)
			# print('data[variable_name].shape ', data[variable_name].shape)
		
		
		variable_names = [	'two_phase_condensate',
							'partial_engulfment',
							'phase_diagram' ]
		objects = [	PlotPhaseDiagramConcDependence().phase_diagram_two_condensates,
					PlotPhaseDiagramConcDependence().phase_diagram_partial_engulfment,
					PlotPhaseDiagramConcDependence().phase_diagram]
		for obj, variable_name in zip( objects, variable_names ):
			data[variable_name] = np.fliplr(obj)[:-1,:].flatten()
		
		self.data = data


	def plot( self, \
		x1title = 'Number of STG bound to one PSD-95',\
		x2title = 'Number of GluN2B bound to one PSD-95', \
		target_type = 'phase_diagram'  ):
		#x2title = 'Raio of PSD95 bound to both GluN2B and STG'
		#x2title = 'Number of GluN2B bound to one CaMKII'
		
		x1data = self.data[x1title]
		x2data = self.data[x2title]
		min_x1, max_x1 = x1data.min() * 0.9, x1data.max() * 1.1
		min_x2, max_x2 = x2data.min() * 0.9, x2data.max() * 1.1
		cols = set_colors( self.data, target_type = 'phase_diagram' )
		#cols = set_colors( self.data, target_type = 'two_phase' )
		
		
		
		min_x1, min_x2 = 0, 0
		
		fig, ax = prepare_plot()
		ax.set_xlabel(x1title)
		ax.set_ylabel(x2title)
		for x1, x2, c in zip(x1data.tolist() , x2data.tolist() , cols):
			ax.plot(x1, x2,'o', \
				markersize = 5, \
				color = 'k', \
				markerfacecolor = c )
		
		ax.set_xlim([min_x1, max_x1])
		ax.set_ylim([min_x2, max_x2])
		
		#
		self.x1title = x1title
		self.x2title = x2title
		self.x1data  = x1data
		self.x2data  = x2data
		self.target_type = target_type
		self.min_x1, self.max_x1, self.min_x2, self.max_x2 = \
			min_x1, max_x1, min_x2, max_x2
		self.fig = fig
		self.ax  = ax
		
		#
		
	def overlay_decision_tree( self, max_depth=3 ):
		
		X = np.vstack([self.x1data,self.x2data])
		Y = self.data[self.target_type]
		
		clf = DecisionTreeClassifier(max_depth=max_depth)
		clf = clf.fit(X.T, Y.T)
		xx1, xx2 = np.meshgrid(
	        np.arange(self.min_x1, self.max_x1, 0.01), np.arange(self.min_x2, self.max_x2, 0.01)
	        )
		mesh_for_prediction = np.array([xx1.ravel(), xx2.ravel()]).T
		
		print("accuracy:{:.3f}".format(clf.score(X.T, Y.T)))
		export_tree_graph(clf)
		pred = clf.predict(mesh_for_prediction)
		
		cs = self.ax.contourf(xx1, xx2, pred.reshape(xx1.shape), cmap = 'binary', vmin=0, vmax=10, zorder=-1) # alpha=0.3
		self.ax.contour(cs, colors='k',linewidths=0.5)
		#self.ax.contour(cs, colors='gray',linestyles = 'dashed', linewidths=1.5)
		
	def save_plots( self ):
		
		t = SpecDatasets()
		t.conc_dependence_merged()
		dir_imgs = os.path.join(t.dir_imgs_root, 'fitting')
		os.makedirs(dir_imgs, exist_ok=True)
		savename = self.x1title.replace(' ', '_')+'_'+self.x2title.replace(' ', '_')
		
		self.fig.savefig( os.path.join(dir_imgs, savename + '.svg' ) )
		self.fig.savefig( os.path.join(dir_imgs, savename + '.png' ) , dpi=150)
		plt.show()
		plt.clf()
		plt.close(fig=self.fig)
		
		
		
if __name__ == '__main__':
	'''
	obj = Fit2DPhaseDiagramConcDependence()
	obj.plot(
		x1title = 'Number of STG bound to one PSD-95', \
		x2title = 'Number of GluN2B bound to one CaMKII', \
		target_type =  'two_phase_condensate' #
		)
	obj.overlay_decision_tree()
	obj.save_plots()
	'''
	
	obj = Fit2DPhaseDiagramConcDependence()
	obj.plot(
		x1title = 'Number of STG bound to one PSD-95', \
		x2title = 'Number of GluN2B bound to one PSD-95', \
		target_type =  'two_phase_condensate' #
		)
	obj.overlay_decision_tree()
	obj.save_plots()
	
	
