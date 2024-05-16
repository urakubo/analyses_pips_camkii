import os, glob, pickle, pprint, copy
import numpy as np



bsize     = 120
space     = [bsize, bsize, bsize]
space_np  = np.array(space)
center_np = space_np / 2

edge0 =  list(range(-int(space[0]/2), int(space[0]/2), 1))
edge1 =  list(range(-int(space[1]/2), int(space[1]/2), 1))
edge2 =  list(range(-int(space[2]/2), int(space[2]/2), 1))


edges0 =  list(range(-int(space[0]/2), int(space[0]/2) + 1, 1))
edges1 =  list(range(-int(space[1]/2), int(space[1]/2) + 1, 1))
edges2 =  list(range(-int(space[2]/2), int(space[2]/2) + 1, 1))


reference_molecule_for_centering = 'All'

subunits = \
	{'GluN2B binding site'	:{'id':3},\
	'CaMKII hub'    :{'id':0},\
	'CaMKII binding site':{'id':5},\
	'STG binding site' :{'id':4},\
	'STG hub'         :{'id':2},\
	'PSD'			:{'id':1}}

molecules_with_all = \
	{'CaMKII'	:{'s':['CaMKII binding site', 'CaMKII hub']	,'c':'#228B22'},\
	'GluN2B'	:{'s':['GluN2B binding site']			,'c':'#ED0DD9'},\
	'PSD95'		:{'s':['PSD']				,'c':'#00FFFF'},\
	'STG'		:{'s':['STG binding site','STG hub']	,'c':'r'},\
	'All'		:{'s':['GluN2B binding site','CaMKII hub','CaMKII binding site','STG binding site','STG hub','PSD']	,'c':'k'}}

molecules_without_all = \
	{'CaMKII'	:{'s':['CaMKII binding site', 'CaMKII hub']	,'c':'#228B22'},\
	'GluN2B'	:{'s':['GluN2B binding site']			,'c':'#ED0DD9'},\
	'PSD95'		:{'s':['PSD']							,'c':'#00FFFF'},\
	'STG'		:{'s':['STG binding site','STG hub']	,'c':'r'}}

for k, v in molecules_with_all.items():
	molecules_with_all[k]['id'] = [subunits[s]['id'] for s in v['s']]


for k, v in molecules_without_all.items():
	molecules_without_all[k]['id'] = [subunits[s]['id'] for s in v['s']]
	

# RDF parameters

rdf_bins = np.arange(0, 40) # You can set "np.arange(0, 35, 2)"
rdf_num_sampling_frames = 5
rdf_sampling_interval   = 5 # 2


rdf_targets_multi_graph = ['Shared PSD95']

# Matplotlib
rc_param = {'pdf.fonttype' : 'truetype',
	'svg.fonttype' : 'none',
	'font.family' : 'sans-serif',
	'font.sans-serif' : 'Arial',
	'font.style' : 'normal',
	'legend.frameon': False}

m = [-1, 0, 1]
neighbors26 = np.stack(np.meshgrid(m, m, m), axis=-1).reshape(-1, 3)
neighbors26 = np.delete(neighbors26, 13, axis=0)


# Conc dependence

STGs    = [108,216,432,576,864,1728,2592,3456,4320,5184]
GluN2Bs = [270,540,1080,2160,4320,6480,8640,12960,17280]
volume  = np.prod(space_np)

STGs    = [ s / volume for s in STGs    ]
GluN2Bs = [ n / volume for n in GluN2Bs ]


# valency-length

valencies = list(range(2,14,2)) 
lengths   = [1, 2, 3, 4, 5, 6, 9]

valencies2 = list(range(4,14,2)) 
lengths2   = [1, 2, 3, 4, 5, 6, 9]

fnames_valency = {v: str(v).zfill(2) for v in valencies }
fnames_length  = {l: str(i).zfill(3) for i,l in enumerate(lengths) }

fnames_valency2 = {v: str(v).zfill(2) for v in valencies2 }
fnames_length2  = {l: str(i).zfill(3) for i,l in enumerate(lengths2) }


