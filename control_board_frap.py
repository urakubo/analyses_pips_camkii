from bin.valency_length_FRAP1_simulate import SimulatePhotobleach
from bin.valency_length_FRAP2_calc_average import CalcAverageFRAPValencyLength
from bin.valency_length_FRAP3_fit_plot import PlotFRAPMatrixValencyLength, PlotFRAPSingleValencyLength


from bin.valency_length1_plot_phase_diagram import \
	PlotConnectivityValencyLength, \
	PlotPhaseDiagramValencyLength, \
	PlotRelaxzationTimeValencyLength, \
	PlotPropertiesValencyLength

from specification_datasets import SpecDatasets
import lib.utils as utils
import pprint
import numpy as np


# Single profiles of FRAP
target_molecule = 'CaMKII' # 'GluN2B', 'CaMKII'
valency         = 12 
length          = 5 # 2, 5, 6

obj = PlotFRAPSingleValencyLength(target_molecule)
#obj.CG_valency_length_only_local_move_fine_sampling(frap = True)
obj.CG_valency_length_only_local_move(frap= True)
obj.suffix_ = 'FRAP_'
obj.plot(valency, length)
obj.save()

'''



## CaMKII FRAP CG_valency_length_only_local_move

obj = SimulatePhotobleach()
obj.CG_valency_length_only_local_move(frap = True)
#obj.CG_valency_length_only_local_move_fine_sampling(frap = True)
#obj.num_skip_frames = 10
obj.repeat_runs()



obj = PlotFRAPMatrixValencyLength( target_molecule = 'CaMKII' ) # CaMKII, GluN2B
obj.CG_valency_length_only_local_move(frap = True)
#obj.CG_valency_length_only_local_move_fine_sampling(frap = True)
obj.run()
obj.save()
obj.save_taus_fitting()



import numpy as np
pl = PlotRelaxzationTimeValencyLength()
pl.CG_valency_length_only_local_move(frap= True)
pl.prefix_loadname = 'FRAP_taus_CaMKII_merged'
pl.basename = 'FRAP_taus_CaMKII_merged'
pl.levels   = np.linspace(-4,1,8)
pl.ticks_level = [-4, -3, -2, -1, 0, 1]
pl.plot_mixture()
pl.save()


## GluN2B FRAP CG_valency_length_only_local_move

obj = PlotFRAPMatrixValencyLength( target_molecule = 'GluN2B' ) # CaMKII, GluN2B
obj.CG_valency_length_only_local_move(frap = True)
#obj.CG_valency_length_only_local_move_fine_sampling(frap = True)
obj.run()
obj.save()
obj.save_taus_fitting()

obj = SimulatePhotobleach()
#obj.CG_valency_length_only_local_move(frap = True)
obj.CG_valency_length_only_local_move_fine_sampling(frap = True)
obj.num_skip_frames = 10
obj.repeat_runs()


# Single profiles of FRAP
target_molecule = 'GluN2B' # 'GluN2B', 'CaMKII'
valency         = 12 # 2, 6
length          = 6 # 2, 6


obj = PlotFRAPMatrixValencyLength( target_molecule = 'GluN2B' ) # CaMKII, GluN2B
obj.CG_valency_length_only_local_move(frap = True)
#obj.CG_valency_length_only_local_move_fine_sampling(frap = True)
obj.run()
obj.save()
obj.save_taus_fitting()


# Single profiles of FRAP
target_molecule = 'CaMKII' # 'GluN2B', 'CaMKII'
valency         = 12 # 2, 6
length          = 6 # 2, 6

obj = PlotFRAPSingleValencyLength(target_molecule)
obj.CG_valency_length_only_local_move_fine_sampling(frap = True)
#obj.CG_valency_length_only_local_move(frap= True)
obj.suffix_ = 'FRAP_'
obj.plot(valency, length)
obj.save()

obj = PlotFRAPMatrixValencyLength( target_molecule = 'CaMKII' ) # CaMKII, GluN2B
#obj.CG_valency_length_only_local_move(frap = True)
obj.CG_valency_length_only_local_move_fine_sampling(frap = True)
obj.run()
obj.save()
obj.save_taus_fitting()


# Merge two taus
prefix = 'FRAP_taus_CaMKII'
suffix = 'matrix'
obj = SpecDatasets()
obj.CG_valency_length_only_local_move(frap= True)
taus_std = utils.load(obj.dir_edited_data, prefix, suffix)
obj.CG_valency_length_only_local_move_fine_sampling(frap= True)
taus_fine = utils.load(obj.dir_edited_data, prefix, suffix)
for i, column in enumerate( obj.lengths ):
	for j, row in enumerate( obj.valencies ):
		filename = obj.filename_edited_matrix(row,column)
		print(filename)
		if ('_006' in filename) or ('04_' in filename):
		#if  ('04_' in filename):
			taus_std[j,i] =taus_fine[j,i]
	
obj.CG_valency_length_only_local_move(frap= True)
prefix = 'FRAP_taus_CaMKII_merged'
pprint.pprint(taus_std)
utils.save( obj.dir_edited_data, prefix, suffix, taus_std )


import numpy as np
pl = PlotRelaxzationTimeValencyLength()
pl.CG_valency_length_only_local_move(frap= True)
pl.prefix_loadname = 'FRAP_taus_CaMKII_merged'
pl.basename = 'FRAP_taus_CaMKII_merged'
pl.levels   = np.linspace(-4,1,8)
pl.ticks_level = [-4, -3, -2, -1, 0, 1]
pl.plot_mixture()
pl.save()

# Single profiles of FRAP
target_molecule = 'CaMKII' # 'GluN2B', 'CaMKII'
valency         = 12 # 2, 6
length          = 6 # 2, 6

obj = PlotFRAPSingleValencyLength(target_molecule)
obj.CG_valency_length_only_local_move_fine_sampling(frap = True)
#obj.CG_valency_length_only_local_move(frap= True)
obj.suffix_ = 'FRAP_'
obj.plot(valency, length)
obj.save()


'''



