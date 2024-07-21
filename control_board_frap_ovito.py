from bin.valency_length_FRAP_make_video_ovito import MakeOvitoVideoFRAP


'''
# Dispersed
obj = MakeOvitoVideoFRAP()
obj.C_valency_length_FRAP_Control_fine_sampling(frap = True)
obj.inspect()
target_molecule = 'Both_condensate_diluent' # 'CaMKII', 'GluN2B', 'Both', 'Both_condensate_diluent'
i = 29 # 25: 12_002, 29: 12_006
obj.num_skip_frames_for_sampling = 10 # 10, 100
obj.run(i, target_molecule)

obj = MakeOvitoVideoFRAP()
#obj.CG_valency_length_only_local_move_fine_sampling(frap = True)
i = 28 # 25: 12_002, 28: 12_005, 29: 12_006
target_molecule = 'CaMKII' # 'CaMKII', 'GluN2B', 'Both'
obj.num_skip_frames_for_sampling = 100
obj.run(i, target_molecule)


'''

obj = MakeOvitoVideoFRAP()
obj.CG_valency_length_only_local_move(frap = True)
i = 28 # 25: 12_002, 28: 12_005
target_molecule = 'CaMKII' 
obj.num_skip_frames_for_sampling = 1 # 1 or 2
obj.run(i, target_molecule)


