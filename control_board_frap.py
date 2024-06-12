from bin.valency_length_CG3_plot_FRAP_video_3d_ovito import MakeOvitoVideoFRAP
from bin.valency_length_FRAP1_simulate import SimulatePhotobleach
from bin.valency_length_FRAP2_fit_plot import PlotFRAPMatrixValencyLength, PlotFRAPSingleValencyLength


'''
centers  = [[0,0,0],[60,0,0],[0,60,0],[0,0,60],[60,60,0],[0,60,60],[60,0,60],[60,60,60]]
suffixes = ['FRAP{}'.format(i) for i in range(8)]
for center, suffix in zip(centers, suffixes):
	obj = SimulatePhotobleach()
	obj.C_valency_length_FRAP_Control()
	obj.center = center
	obj.suffix = suffix
	obj.repeat_runs()

from specification_datasets import SpecDatasets
import lib.utils as utils
import numpy as np



class CalcAverageFRAPValencyLength(SpecDatasets):
		
	def __init__( self ):
		#SpecDatasetsFRAP.__init__(self, target_dataset )
		self.center = None
		self.suffix_ = 'FRAP_'
		self.suffixes = ['FRAP{}'.format(i) for i in range(8)]
		
	def convert_save( self ):
		for valency in self.valencies_frap:
			for length in self.lengths_frap:
				prefix = self.filename_edited_matrix(valency, length)
				print( prefix )
				molecular_concentrations = []
				for suffix in self.suffixes:
					d = utils.load(self.dir_edited_data, prefix, suffix)
					molecular_concentrations.append( d['molecular_concentration_in_target_area'] )
				
				d_ = {}
				d_['time_steps'] = d['time_steps']
				d_['molecular_concentration_in_target_area'] = np.mean( molecular_concentrations ,  axis=0 )
				d_['std_molecular_concentration_in_target_area']     = np.std(  molecular_concentrations ,  axis=0 )
				
				utils.save(self.dir_edited_data, prefix, self.suffix_, d_)

obj = CalcAverageFRAPValencyLength()
obj.C_valency_length_FRAP_Control()
obj.convert_save()

# Single profiles of FRAP
target_molecule = 'CaMKII' # 'GluN2B', 'CaMKII'
valency         = 12 # 2, 6
length          = 6 # 2, 6

obj = PlotFRAPSingleValencyLength(target_molecule)
obj.C_valency_length_FRAP_Control()
obj.suffix_ = 'FRAP_'
obj.plot(valency, length)
obj.save()

'''

obj = PlotFRAPMatrixValencyLength()
obj.C_valency_length_FRAP_Control()
obj.suffix_ = 'FRAP_'
obj.run(frap= True)
obj.save()




# Two FRAP files are merged and saved.
'''
prefix = 'FRAP'
suffix = 'matrix'
dir_edited_data = os.path.join('data4', 'small_colony2')
d_colony2 = utils.load(dir_edited_data, prefix, suffix)
dir_edited_data = os.path.join('data4', 'small_colony3')
d_colony3 = utils.load(dir_edited_data, prefix, suffix)

for valency_length in ['12_000','10_000','08_000','12_001','10_001','12_002']:
	d_colony3[valency_length] = d_colony2[valency_length]

dir_edited_data = os.path.join('data4', 'small_colony2')
prefix = 'FRAP_merged'
utils.save(dir_edited_data, prefix, suffix, d_colony3)
'''
