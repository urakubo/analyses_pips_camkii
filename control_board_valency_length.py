from bin.both2_plot_matrix_3d_pyvista import PlotMatrixPyvista
from bin.both1_plot_matrix import PlotConcMatrix, PlotConnectivityMatrix

from bin.valency_length1_plot_phase_diagram import \
	PlotConnectivityValencyLength, \
	PlotPhaseDiagramValencyLength, \
	PlotRelaxzationTimeValencyLength, \
	PlotPropertiesValencyLength


from bin.valency_length_CG3_plot_FRAP_video_3d_ovito import MakeOvitoVideoFRAP


'''


#species, type_analysis = 'CaMKII', 'conc_in_CaMKII_condensate'
#species, type_analysis = 'GluN2B', 'conc_in_CaMKII_condensate'
species, type_analysis = 'All', 'conc_in_CaMKII_condensate'
obj = PlotConnectivityValencyLength(species, type_analysis)
obj.CG_valency_length_only_local_move(frap = True)
obj.run()
obj.save()


obj = PlotMatrixPyvista()
#obj.CG_valency_length(sub = True)
obj.CG_valency_length()
obj.run(non_rotated = True)

obj = PlotMatrixPyvista()
#obj.valency_length(sub = True)
obj.valency_length()
obj.run(non_rotated = True)

species, type_analysis = 'CaMKII', 'average'
#species, type_analysis = 'PSD95' , 'average'
#species, type_analysis = 'PSD95' , 'ratio'
pl = PlotConnectivityValencyLength(species, type_analysis)
pl.valency_length()
pl.run()
pl.save()


#species, type_analysis = 'CaMKII', 'average'
#species, type_analysis = 'PSD95' , 'average'
species, type_analysis = 'PSD95' , 'ratio'
pl = PlotConnectivityMatrix(species, type_analysis)
pl.valency_length()
pl.run()
pl.save()



pl = PlotPhaseDiagramValencyLength()
pl.valency_length() # Only for the save directory?
pl.plot()
pl.save()


obj = PlotMatrixPyvista()
obj.valency_length(sub = True)
#obj.valency_length()
obj.run()


'''

