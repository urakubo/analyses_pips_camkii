
from bin.connectivity_make_video_ovito_modularity import MakeOvitoVideoModularity
from bin.connectivity_calc1_graph_time_development import MakeConnectivityGraphTimeDev
from bin.connectivity_plot_matrix_dendrogram import PlotConnectivityMatrixDendrogram


from bin.connectivity_plot_profile_3d_ovito import Plot3dOvitoConnectivity

from bin.valency_length1_plot_phase_diagram import \
	PlotConnectivityValencyLength, \
	PlotPhaseDiagramValencyLength, \
	PlotRelaxzationTimeValencyLength, \
	PlotPropertiesValencyLength

# from bin.connectivity_calc2_modularity import CalcModularityClustring




obj = Plot3dOvitoConnectivity()
obj.CG_valency_length_only_local_move(frap= True)
#obj.inspect()
obj.plot_an_image(28) # 25, 29



'''
obj = PlotConnectivityMatrixDendrogram()
obj.CG_valency_length_only_local_move(frap= True)
#obj.inspect()
obj.plot_connectivity_save_graph(28) # 25, 29


property = 'clustering_log' # 'density', 'modularity', 'clustering', 'clustering_log', 'FRAP'
obj = PlotPropertiesValencyLength(property)
obj.CG_valency_length_only_local_move(frap= True)
obj.plot2()
obj.save()



obj = MakeConnectivityGraphTimeDev()
obj.valency_length_small_colony2()
obj.inspect()

# repeat_for_valency_length(filenames_lammpstrj, filenames_edited)

i = 7*5+1 # val_12\R2_002
i = 7*4+0 # val_10\R2_000
i = 7*5+0 # val_12\R2_000
i = 7*4+1 # val_10\R2_001
#i = 7*3+0 # val_08\R2_000
obj.repeat_for_time_development( i )
#obj.repeat_for_valency_length()


obj = MakeOvitoVideoModularity()
obj.valency_length_small_colony2()
obj.inspect()
#i = 7*5+2 # val_12\R2_002
i = 7*5+6 # val_12\R2_006
obj.run( i = i )
obj.make_a_video( i )




obj = CalcModularityClustringC()
obj.CG_valency_length_only_local_move(frap = True)
obj.run()
obj.save()

'''



