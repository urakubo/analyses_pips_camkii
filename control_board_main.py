

from bin.all1_edit_data import EditData
from bin.all2_edit_connectivity_graph import EditConnectivityGraph
from bin.all3_plot_profile import PlotProfiles
from bin.all4_plot_profile_3d_ovito import Plot3dOvito
from bin.all5_plot_profile_3d_pyvista import Plot3dPyvista
from bin.all6_plot_profile_shared_PSD95 import PlotProfilesSharedPSD95
from bin.all7_make_video_3d_ovito import MakeOvitoVideo

from bin.both1_plot_matrix import PlotConcMatrix, PlotConnectivityMatrix
from bin.both2_plot_matrix_3d_pyvista import PlotMatrixPyvista

from bin.conc_dependence1_calc_connectivity_plot_phase_diagram import \
	HandleConnectivityPhaseDiagramConcDependence, \
	HandleCondVolumePhaseDiagramConcDependence, \
	PlotPhaseDiagramConcDependence


obj = Plot3dPyvista()
#obj.boundary_conditions2() #  conc_dependence_merged(), valency_length(), CG_valency_length(), boundary_conditions2()
obj.boundary_conditions2()
obj.run()



'''


obj = EditConnectivityGraph()
#obj.conc_dependence()
obj.conc_dependence_033()
obj.run()

obj = PlotProfiles()
obj.CG_valency_length_only_local_move(frap = True)
obj.run()


obj = Plot3dOvito()
#obj.C_valency_length_FRAP_Control(frap = True)
obj.CG_valency_length_only_local_move(frap = True)
obj.repeat_run()


# 'region_condensates', 'conc_CaMKII', 'conc_PSD95', 'conc_STG', 'conc_GluN2B', 'rdf',  'rdf_PSD95'
# 'concs_in_CaMKII', 'concs_in_STG',
# 'shared_PSD95', 'unshared_PSD95', 'conc_unrotated_CaMKII'

obj = Plot3dPyvista()
#obj.boundary_conditions2() #  conc_dependence_merged(), valency_length(), CG_valency_length(), boundary_conditions2()
obj.CG_valency_length_only_local_move(frap = True)
obj.run()


target = 'rdf_CG'
obj = PlotConcMatrix(target=target)
#obj.conc_dependence_merged()
obj.CG_valency_length_only_local_move(frap = True)
obj.run()
obj.save()


obj = EditData()
obj.CG_valency_length_only_local_move(frap = True)
obj.run()

obj = EditConnectivityGraph()
obj.CG_valency_length_only_local_move(frap= True)
obj.run()

obj = PlotProfiles()
obj.CaMKII_blocked4()
obj.run()


obj = Plot3dOvito()
obj.boundary_conditions2() #  conc_dependence(), valency_length(), CG_valency_length(), boundary_conditions2()
obj.repeat_run()

obj = Plot3dPyvista()
obj.boundary_conditions2() #  conc_dependence_merged(), valency_length(), CG_valency_length(), boundary_conditions2()
obj.run()

obj = PlotProfilesSharedPSD95()
obj.valency_length() #  conc_dependence(), valency_length(), CG_valency_length(), boundary_conditions2()
obj.run()

sub = True
obj = PlotMatrixPyvista()
obj.conc_dependence_merged(sub = sub)
obj.run()


# 'region_condensates', 'conc_CaMKII', 'conc_PSD95', 'conc_STG', 'conc_GluN2B', 'rdf',  'rdf_PSD95'
# 'concs_in_CaMKII', 'concs_in_STG',
# 'shared_PSD95', 'unshared_PSD95', 'conc_unrotated_CaMKII'

target = 'conc_unrotated_CaMKII'
obj = PlotConcMatrix(target)
obj.conc_dependence_merged()
#obj.valency_length()
obj.run()
obj.save()


species       = 'CaMKII' # 'STG','GluN2B', 'PSD95','CaMKII'
type_analysis = 'distribution'
# 'average' and 'distribution' for all,
# species: 'GluN2B', type_analysis 'CaMKII' or 'PSD95'
# species: 'PSD95' , type_analysis 'ratio'

obj = PlotConnectivityMatrix(species, type_analysis)
obj.valency_length()
values = obj.run()
obj.save()



obj = MakeOvitoVideo()
obj.boundary_conditions2()
obj.inspect()
obj.run(3)
obj.make_a_video(3)


'''


