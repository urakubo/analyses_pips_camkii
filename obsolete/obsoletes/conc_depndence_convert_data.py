
import os, glob, pickle, pprint, copy
import numpy as np

import lib.utils as utils
import lib.parameters as p
import pyvista


STGs    = [108,216,432,864,1728,2592,3456,4320,5184]
id_STG1 = {0: 0, 1: 1, 2: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8}

GluN2Bs = [270,540,1080,2160,4320,6480,8640,12960,17280]
iGluN2Bs= [0  ,1   ,2   ,3   ,4   ,5   ,6   ,7   ,8]
volume  = np.prod(space_np)

STGs    = [ s / volume for s in STGs    ]
GluN2Bs = [ n / volume for n in GluN2Bs ]

dir_target1       = 'conc_dependence'
dir_target2       = 'conc_dependence_0.33'
dir_target_outout = 'conc_dependence_merged'
suffix            = 'sigma_2'

dir_edited_data1	= os.path.join('data4', dir_target1)
dir_edited_data2	= os.path.join('data4', dir_target2)
dir_edited_datao	= os.path.join('data4', dir_target_output)
os.makedirs(dir_edited_datao, exist_ok=True)


for id_glun2b in range(9):
	for id_stg in range(9):
		if id_stg != 3:
			id     = id_STG1[id_stg] + id_glun2b * len(id_STG1)
			prefix = str(id).zfill(3)
			d      = utils.load(dir_edited_data1, prefix, suffix)
		else:
			prefix = str(id_glun2b).zfill(3)
			d      = utils.load(dir_edited_data2, prefix, suffix)
		
		prefix = str(id_stg).zfill(2)+'_'+str(id_glun2b).zfill(2)
		utils.save(dir_edited_datao, prefix, suffix)


