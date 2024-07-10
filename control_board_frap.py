from bin.valency_length_CG3_plot_FRAP_video_3d_ovito import MakeOvitoVideoFRAP
from bin.valency_length_FRAP1_simulate import SimulatePhotobleach
from bin.valency_length_FRAP2_calc_average import CalcAverageFRAPValencyLength
from bin.valency_length_FRAP3_fit_plot import PlotFRAPMatrixValencyLength, PlotFRAPSingleValencyLength


from bin.valency_length1_plot_phase_diagram import \
	PlotConnectivityValencyLength, \
	PlotPhaseDiagramValencyLength, \
	PlotRelaxzationTimeValencyLength, \
	PlotPropertiesValencyLength


from bin.valency_length_CG3_plot_FRAP_video_3d_ovito import MakeOvitoVideoFRAP

from specification_datasets import SpecDatasets
import lib.utils as utils
import pprint
import numpy as np

obj = PlotFRAPSingleValencyLength()
obj.C_valency_length_FRAP_Control_fine_sampling(frap = True)
obj.plot_dissociated_CaMKII()
obj.save()


'''

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

# Two FRAP files are merged and saved.prefix = 'FRAP'



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



target_molecule = 'CaMKII' # 'CaMKII', 'GluN2B', 'Both'
obj = MakeOvitoVideoFRAP()
#obj.valency_length_small_colony2()
#obj.inspect()
obj.valency_length_small_colony3()
i = 7*5+6 # val_12\R2_006
i = 7*5+2 # val_12\R2_006
obj.run(i, target_molecule)
obj.make_a_video(i)


## CaMKII FRAP CG_valency_length_only_local_move

obj = SimulatePhotobleach()
#obj.CG_valency_length_only_local_move(frap = True)
#obj.CG_valency_length_only_local_move_fine_sampling(frap = True)

#obj.num_skip_frames = 10
obj.repeat_runs()

obj = PlotFRAPMatrixValencyLength( target_molecule = 'CaMKII' ) # CaMKII, GluN2B
obj.CG_valency_length_only_local_move(frap = True)
#obj.CG_valency_length_only_local_move_fine_sampling(frap = True)
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
		#if ('_006' in filename) or ('04_' in filename):
		if  ('04_' in filename):
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



