

from bin.both1_plot_matrix import PlotConcMatrix, PlotConnectivityMatrix
from bin.both2_plot_matrix_3d_pyvista import PlotMatrixPyvista

from bin.conc_dependence1_calc_connectivity_plot_phase_diagram import \
	HandleConnectivityPhaseDiagramConcDependence, \
	HandleCondVolumePhaseDiagramConcDependence, \
	PlotPhaseDiagramConcDependence



species       = 'CaMKII' # 'STG','GluN2B', 'PSD95','CaMKII'
type_analysis = 'distribution'
# 'average' and 'distribution' for all,
# species: 'GluN2B', type_analysis 'CaMKII' or 'PSD95'
# species: 'PSD95' , type_analysis 'ratio'

obj = PlotConnectivityMatrix(species, type_analysis)
obj.valency_length()
values = obj.run()
obj.save()


'''


sub = True
obj = PlotMatrixValencyLengthPyvista()
obj.valency_length()
obj.run(sub = sub)

sub = True
obj = PlotMatrixConcPyvista()
obj.conc_dependence_merged()
obj.run(sub = sub)



# 'region_condensates', 'conc_CaMKII', 'conc_PSD95', 'conc_STG', 'conc_GluN2B', 'rdf',  'rdf_PSD95'
# 'concs_in_CaMKII', 'concs_in_STG',
# 'shared_PSD95', 'unshared_PSD95', 'conc_unrotated_CaMKII'

target = 'conc_CaMKII'
obj = PlotConcMatrix(target)
#obj.valency_length()
obj.conc_dependence_merged()
obj.run()
obj.save()



pl = PlotPhaseDiagramConcDependence()
pl.conc_dependence_merged()
pl.plot()
pl.save_plots()



species, type_analysis = 'CaMKII', 'average'
#species, type_analysis = 'PSD95' , 'average'
#species, type_analysis = 'PSD95' , 'ratio'
# species, type_analysis = 'PSD95' , 'average_GluN2B'
pl = HandleConnectivityPhaseDiagramConcDependence(species, type_analysis)
pl.conc_dependence_merged()
#pl.edit_data_save_them()
pl.load_data()
pl.plot_data()
pl.save_plots()



species = 'CaMKII' # 'CaMKII', 'STG'
pl = HandleCondVolumePhaseDiagramConcDependence(species)
pl.conc_dependence_merged()
#pl.edit_data_save_them()
pl.load_data()
pl.plot_data()
pl.save_plots()


'''


