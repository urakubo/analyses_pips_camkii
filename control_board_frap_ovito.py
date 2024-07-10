from bin.valency_length_FRAP_make_video_ovito import MakeOvitoVideoFRAP



# Small colony 2
obj = MakeOvitoVideoFRAP()
obj.C_valency_length_FRAP_Control_fine_sampling(frap = True)
obj.inspect()
i = 29 # 25: 12_002, 29: 12_006
target_molecule = 'Both_condensate_diluent' # 'CaMKII', 'GluN2B', 'Both', 'Both_condensate_diluent'
obj.num_skip_frames_for_sampling = 100
obj.run(i, target_molecule)
#obj.make_a_video(i)
