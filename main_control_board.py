

from bin.all1_edit_data import EditData
from bin.all2_edit_connectivity_graph import EditConnectivityGraph
from bin.all3_plot_profile import PlotProfiles
from bin.all4_plot_profile_3d_ovito import Plot3dOvito
from bin.all5_plot_profile_3d_pyvista import Plot3dPyvista
from bin.all6_plot_profile_shared_PSD95 import PlotProfilesSharedPSD95
from bin.all7_make_video_3d_ovito import MakeOvitoVideo

from bin.valency_length_CG3_plot_FRAP_video_3d_ovito import MakeOvitoVideoFRAP
from bin.valency_length_FRAP1_simulate import SimulatePhotobleach
from bin.valency_length_FRAP2_fit_plot import PlotFRAPMatrixValencyLength, PlotFRAPSingleValencyLength


from bin.both1_plot_matrix import PlotConcMatrixValencyLength, PlotConcMatrixConcDependence
from bin.both2_plot_matrix_3d_pyvista import PlotMatrixValencyLengthPyvista, PlotMatrixConcPyvista

from bin.surface_tension1_calc_radius_CaMKII_condensate import HandleRDPCaMKII
from bin.surface_tension2_plot import PlotSurfaceTension
from bin.surface_tension3_plot_phase_diagram import PlotSurfaceTensionPhaseDiagramValencyLength

from bin.conc_dependence1_calc_connectivity_plot_phase_diagram import HandleConnectivityPhaseDiagramConcDependence, HandleCondVolumePhaseDiagramConcDependence, PlotPhaseDiagramConcDependence

from bin.valency_length1_plot_phase_diagram import PlotConnectivityValencyLength, PlotPhaseDiagramValencyLength, PlotRelaxzationTimeForMixtureValencyLength, PlotPropertiesValencyLength


#species, type_analysis = 'CaMKII', 'average'
#species, type_analysis = 'PSD95' , 'average'
species, type_analysis = 'PSD95' , 'ratio'
pl = PlotConnectivityValencyLength(species, type_analysis)
pl.valency_length()
pl.run()
pl.save()


'''



centers  = [[0,0,0],[60,0,0],[0,60,0],[0,0,60],[60,60,0],[0,60,60],[60,0,60],[60,60,60]]
suffixes = ['FRAP{}'.format(i) for i in range(8)]
for center, suffix in zip(centers, suffixes):
	obj = SimulatePhotobleach()
	obj.C_valency_length_FRAP_Control()
	obj.center = center
	obj.suffix = suffix
	obj.repeat_runs()




pl = PlotPhaseDiagramValencyLength()
pl.valency_length() # Only for the save directory?
pl.plot()
pl.save()


#species, type_analysis = 'CaMKII', 'average'
#species, type_analysis = 'PSD95' , 'average'
species, type_analysis = 'PSD95' , 'ratio'
pl = PlotConnectivityValencyLength(species, type_analysis)
pl.valency_length()
pl.run()
pl.save()

pl = PlotRelaxzationTimeForMixtureValencyLength()
pl.valency_length_small_colony2()
pl.run_mixture()
pl.save()


property = 'modularity' # 'density', 'modularity', 'clustering', 'FRAP'
pl = PlotPropertiesValencyLength(property)
pl.valency_length_small_colony2()
pl.plot()
pl.save()


property = 'FRAP_GluN2B' # 'density', 'modularity', 'clustering', 'FRAP', 'FRAP_GluN2B'
pl = PlotPropertiesValencyLength(property)
pl.valency_length_small_colony3()
pl.plot()
pl.save()


obj = PlotSurfaceTensionPhaseDiagramValencyLength()
obj.CG_valency_length()
obj.reflect_spec()
#graph.run_calc()
#graph.save_data()
obj.load_data()
obj.plot_phase_diagrams()
obj.plot_logistic_regression()


obj = PlotSurfaceTension()
obj.CG_valency_length()
obj.reflect_spec()
# obj.inspect_targets()
obj.multiple_run_plot()
obj.show_polar_graphs()


obj = HandleRDPCaMKII()
obj.CG_valency_length()
obj.run()
obj.save_figs()
obj.save_hmws()


sub = True
obj = PlotMatrixValencyLengthPyvista()
obj.CG_valency_length()
obj.run(sub = sub)

sub = True
obj = PlotMatrixConcPyvista()
obj.conc_dependence_merged()
obj.run(sub = sub)



# 'region_condensates', 'conc_CaMKII', 'conc_PSD95', 'conc_STG', 'conc_GluN2B', 'rdf',  'rdf_PSD95'
# 'concs_in_CaMKII', 'concs_in_STG',
# 'shared_PSD95', 'unshared_PSD95', 'conc_unrotated_CaMKII'

target = 'conc_unrotated_CaMKII'
obj = PlotConcMatrixConcDependence(target)
obj.conc_dependence_merged()
obj.run()
obj.save()


# 'region_condensates', 'conc_CaMKII', 'conc_PSD95', 'conc_STG', 'conc_GluN2B', 'rdf',  'rdf_PSD95'
# 'concs_in_CaMKII', 'concs_in_STG',
# 'shared_PSD95', 'unshared_PSD95', 'conc_unrotated_CaMKII'

target = 'conc_unrotated_CaMKII'
obj = PlotConcMatrixValencyLength(target)
obj.CG_valency_length()
obj.run()
obj.save()


centers  = [[0,0,0],[60,0,0],[0,60,0],[0,0,60],[60,60,0],[0,60,60],[60,0,60],[60,60,60]]
suffixes = ['FRAP{}'.format(i) for i in range(8)]
for center, suffix in zip(centers, suffixes):
	obj = SimulatePhotobleach()
	obj.C_valency_length_FRAP_Control()
	obj.center = center
	obj.suffix = suffix
	obj.repeat_runs()



obj = SimulatePhotobleach()
obj.C_valency_length_FRAP_Control()
obj.center = [0,60,60]
obj.suffix = 'FRAP2'
obj.repeat_runs()


obj = PlotFRAPMatrixValencyLength()
obj.C_valency_length_FRAP_Control()
obj.run(frap= True)
obj.save()

obj = PlotFRAPMatrixValencyLength()
obj.CG_valency_length_only_local_move()
#obj.valency_length_small_colony3()
obj.run(frap= True)
obj.save()



obj = SimulatePhotobleach()
obj.C_valency_length_FRAP_Control()
obj.repeat_runs()

obj = SimulatePhotobleach()
obj.CG_valency_length_only_local_move()
#obj.valency_length_small_colony3()
obj.repeat_runs()


obj = PlotFRAPMatrixValencyLength( target_molecule = 'CaMKII' ) # CaMKII, GluN2B
obj.CG_valency_length_only_local_move()
#obj.valency_length_small_colony3()
obj.run(frap = True)
obj.save()


obj = EditData()
obj.CG_valency_length() #  conc_dependence(), valency_length(), valency_length_CG(), inhibitor()
obj.run()

obj = EditConnectivityGraph()
obj.CG_valency_length() #  conc_dependence(), valency_length(), CG_valency_length(), boundary_conditions2
obj.run()

obj = PlotProfiles()
obj.inhibitor()  #  conc_dependence(), valency_length(), CG_valency_length(), boundary_conditions2()
obj.run()

# I manually ran each one of them,
# because I do not know how to fully reset the ovito visualization system.
obj = Plot3dOvito()
obj.boundary_conditions2() #  conc_dependence(), valency_length(), CG_valency_length(), boundary_conditions2()
obj.inspect()
obj.run(1)

obj = Plot3dPyvista()
obj.boundary_conditions2() #  conc_dependence_merged(), valency_length(), CG_valency_length(), boundary_conditions2()
obj.run()

obj = PlotProfilesSharedPSD95()
obj.valency_length() #  conc_dependence(), valency_length(), CG_valency_length(), boundary_conditions2()
obj.run()


obj = MakeOvitoVideo()
obj.boundary_conditions2()
obj.inspect()
obj.run(3)
obj.make_a_video(3)

target_molecule = 'CaMKII' # 'CaMKII', 'GluN2B', 'Both'
obj = MakeOvitoVideoFRAP()
#obj.valency_length_small_colony2()
#obj.inspect()
obj.valency_length_small_colony3()
i = 7*5+6 # val_12\R2_006
i = 7*5+2 # val_12\R2_006
obj.run(i, target_molecule)
obj.make_a_video(i)



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


