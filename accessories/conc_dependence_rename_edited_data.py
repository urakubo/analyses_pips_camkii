import os, glob, pickle, pprint, copy
import shutil
import numpy as np
import lib.utils as utils


STGs    = [108,216,432,864,1728,2592,3456,4320,5184]
id_STG1 = {0: 0, 1: 1, 2: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8}

GluN2Bs = [270,540,1080,2160,4320,6480,8640,12960,17280]
iGluN2Bs= [0  ,1   ,2   ,3   ,4   ,5   ,6   ,7   ,8]


dir_target1 = 'conc_dependence'
dir_target2 = 'conc_dependence_0.33'
dir_targeto = 'conc_dependence_merged'
suffix = 'sigma_2'
suffix = 'connectivity_graph'

dir_edited_data1 = os.path.join('data4', dir_target1)
dir_edited_data2 = os.path.join('data4', dir_target2)
dir_edited_datao = os.path.join('data4', dir_targeto)
os.makedirs(dir_edited_datao, exist_ok=True)


for id_glun2b in range(9):
	for id_stg in range(10):
		if id_stg != 3:
			id     = id_STG1[id_stg] + id_glun2b * len(id_STG1)
			prefix = str(id).zfill(3)
			original_file = os.path.join(dir_edited_data1,'{}_{}.pickle'.format(prefix, suffix) )
		else:
			prefix = str(id_glun2b).zfill(3)
			original_file = os.path.join(dir_edited_data2,'{}_{}.pickle'.format(prefix, suffix) )
		
		prefix = str(id_stg).zfill(2)+'_'+str(id_glun2b).zfill(2)
		print(prefix)
		destination_file = os.path.join(dir_edited_datao,'{}_{}.pickle'.format(prefix, suffix) )
		shutil.copy(original_file, destination_file)


