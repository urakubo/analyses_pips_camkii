'''
from bin.valency_length_FRAP1_simulate import SimulatePhotobleach
from bin.valency_length_FRAP2_calc_average import CalcAverageFRAPValencyLength
from bin.valency_length_FRAP3_fit_plot import PlotFRAPMatrixValencyLength, PlotFRAPSingleValencyLength
'''

from bin.valency_length1_plot_phase_diagram import \
	PlotConnectivityValencyLength, \
	PlotPhaseDiagramValencyLength, \
	PlotRelaxzationTimeValencyLength, \
	PlotPropertiesValencyLength


from specification_datasets import SpecDatasets
import lib.utils as utils
import os, pprint
import numpy as np
os.environ['OVITO_GUI_MODE'] = '1' # Request a session with OpenGL support



pl = PlotRelaxzationTimeValencyLength()
pl.C_valency_length_FRAP_Control_fine_sampling(frap = True)
pl.prefix_loadname = 'FRAP_taus_CaMKII_ratio'
pl.basename = 'FRAP_taus_CaMKII_ratio'
pl.levels   =  np.linspace(0,2,8)
pl.ticks_level = [0, 1, 2]
pl.plot_mixture()
pl.save()


'''
## Get ratio
obj1 = SpecDatasets()
obj1.CG_valency_length_only_local_move(frap= True)
prefix = 'FRAP_taus_CaMKII_merged'
suffix = 'matrix'
d1 = utils.load( obj1.dir_edited_data, prefix, suffix )

obj2 = SpecDatasets()
obj2.C_valency_length_FRAP_Control_fine_sampling(frap = True)
prefix = 'FRAP_taus_CaMKII'
suffix = 'matrix'
d2 = utils.load( obj2.dir_edited_data, prefix, suffix )

dd = d1 / d2

prefix = 'FRAP_taus_CaMKII_ratio'
pprint.pprint(dd)
utils.save( obj2.dir_edited_data, prefix, suffix, dd )
##


obj = PlotFRAPSingleValencyLength()
obj.C_valency_length_FRAP_Control_fine_sampling(frap = True)
obj.plot_dissociated_CaMKII()
obj.save()

import numpy as np
pl = PlotRelaxzationTimeValencyLength()
pl.C_valency_length_FRAP_Control_fine_sampling(frap = True)
pl.prefix_loadname = 'FRAP_taus_CaMKII'
pl.basename = 'FRAP_taus_CaMKII'
pl.levels   =  np.linspace(-4,1,8)
pl.ticks_level = [-4, -3, -2, -1, 0, 1]
pl.plot_mixture()
pl.save()

centers  = [[0,0,0],[60,0,0],[0,60,0],[0,0,60],[60,60,0],[0,60,60],[60,0,60],[60,60,60]]
suffixes = ['FRAP{}'.format(i) for i in range(8)]

for center, suffix in zip(centers, suffixes):
	obj = SimulatePhotobleach()
	obj.C_valency_length_FRAP_Control_fine_sampling(frap = True)
	#obj.C_valency_length_FRAP_Control(frap = True)
	#obj.num_skip_frames = 10
	obj.center = center
	obj.suffix = suffix
	obj.repeat_runs()


suffixes = ['FRAP{}'.format(i) for i in range(8)]
suffix   = 'FRAP_'
obj = CalcAverageFRAPValencyLength()
obj.suffixes = suffixes
obj.suffix_unified = suffix
obj.C_valency_length_FRAP_Control_fine_sampling(frap = True)
#obj.C_valency_length_FRAP_Control(frap = True)
obj.convert_save()


suffix   = 'FRAP_'
obj = PlotFRAPMatrixValencyLength()
obj.C_valency_length_FRAP_Control_fine_sampling(frap= True)
#obj.C_valency_length_FRAP_Control(frap = True)
obj.suffix = suffix
fitting_tau = obj.run()
obj.save()
obj.save_taus_fitting()



obj = SimulatePhotobleach()
obj.C_valency_length_FRAP_Control_fine_sampling()
obj.center = [0,60,60]
obj.suffix = 'FRAP2'
obj.repeat_runs()


center = [0,0,0]
suffix = 'FRAP0'
obj = SimulatePhotobleach()
obj.C_valency_length_FRAP_Control_fine_sampling(frap = True)
obj.center = center
obj.suffix = suffix
obj.repeat_runs()


obj = PlotFRAPMatrixValencyLength()
obj.C_valency_length_FRAP_Control_fine_sampling(frap= True)
obj.suffix = 'FRAP0'
obj.run()
obj.save()


suffix   = 'FRAP_'
obj = PlotFRAPMatrixValencyLength()
#obj.C_valency_length_FRAP_Control_fine_sampling(frap= True)
obj.C_valency_length_FRAP_Control(frap = True)
obj.suffix = suffix
fitting_tau = obj.run()
obj.save()
obj.save_taus_fitting()



obj.C_valency_length_FRAP_Control(frap= True)
prefix = 'FRAP_taus_CaMKII_merged'
pprint.pprint(taus_std)
utils.save( obj.dir_edited_data, prefix, suffix, taus_std )

import numpy as np
pl = PlotRelaxzationTimeValencyLength()
pl.C_valency_length_FRAP_Control(frap= True)
pl.prefix_loadname = 'FRAP_taus_CaMKII_merged'
pl.basename = 'FRAP_taus_CaMKII_merged'
pl.levels   = np.linspace(-4,1,8)
pl.ticks_level = [-4, -3, -2, -1, 0, 1]
pl.plot_mixture()
pl.save()
'''



