
from bin.surface_tension0_plot_pyvista_CaMKII import PlotPyvistaCaMKII

from bin.surface_tension1_calc_radius_CaMKII_condensate import HandleRDPCaMKII
from bin.surface_tension2_plot import PlotSurfaceTension
from bin.surface_tension3_plot_phase_diagram import PlotSurfaceTensionPhaseDiagramValencyLength


from bin.connectivity_plot_profile_3d_ovito import Plot3dOvitoConnectivity



from bin.both1_plot_matrix import PlotConcMatrix, PlotConnectivityMatrix
from bin.all1_edit_data import EditData
from bin.all3_plot_profile import PlotProfiles



## Example plot for demo in Pyvista
num_samples = 8
prefix, random_seed = '12_005', 0
num_samples = 10
prefix, random_seed = '12_002', 1

obj = PlotPyvistaCaMKII()
obj.CG_valency_length()
obj.plot_save( prefix = prefix, random_seed = random_seed, num_samples = num_samples )
##


'''


# 'region_condensates', 'conc_CaMKII', 'conc_PSD95', 'conc_STG', 'conc_GluN2B', 'rdf',  'rdf_PSD95'
# 'concs_in_CaMKII', 'concs_in_STG',
# 'shared_PSD95', 'unshared_PSD95', 'conc_unrotated_CaMKII'

target = 'rdf'
obj = PlotConcMatrix(target)
obj.CG_valency_length(sub2 = True)
obj.run()
obj.save()


from bin.all4_plot_profile_3d_ovito import Plot3dOvito

# Compare two representative case: spatial arrangement.
obj = PlotSurfaceTension()
obj.CG_valency_length()
obj.apply_specification()
#obj.inspect_targets()
#obj.multiple_run_plot(['12_002','12_006'])
#obj.multiple_run_plot(['12_002','12_005'])
#obj.show_polar_graphs(targets = ['12_001','12_002','12_003','12_004','12_005','12_006'])
#obj.show_polar_graphs(targets = ['12_002','12_005'], mode = 'angle_from_condensate_center')
obj.show_polar_graphs(targets = ['12_002','12_005'], mode = 'angle_from_hub')


obj = Plot3dOvito()
obj.CG_valency_length()
obj.inspect()
#obj.repeat_run()
obj.run(34)

obj = Plot3dOvitoConnectivity()
obj.CG_valency_length()
obj.inspect()
#obj.plot_an_image(37, mode = 'CaMKII_hub_beads') # 37, 41
obj.plot_an_image(34) # 31, 34 # , mode = 'CaMKII_hub_beads'


obj = EditData()
obj.CG_valency_length()
obj.run()



graph = PlotSurfaceTensionPhaseDiagramValencyLength()
graph.CG_valency_length()
graph.reflect_spec()
#graph.run_calc()
#graph.save_data()
graph.load_data()
graph.plot_phase_diagrams()
graph.plot_logistic_regression()


# 'angle_from_condensate_center'  'angle_from_hub' # targets = ['12_006']
# , mode_surrogate = True





# Compare two representative case: spatial arrangement.
obj = PlotSurfaceTension()
obj.CG_valency_length()
obj.apply_specification()
#obj.inspect_targets()
#obj.multiple_run_plot(['12_001','12_002','12_003','12_004','12_005','12_006'])
#obj.show_polar_graphs(targets = ['12_001','12_002','12_003','12_004','12_005','12_006'])
obj.show_polar_graphs(targets = ['12_005'], mode = 'angle_from_condensate_center')

# obj.show_polar_graphs(targets = ['12_005'], mode = 'angle_from_hub')


# 'angle_from_condensate_center'  'angle_from_hub' # targets = ['12_006']
# , mode_surrogate = True



obj = PlotProfiles()
obj.CG_valency_length()
obj.run()


graph = PlotSurfaceTensionPhaseDiagramValencyLength()
graph.CG_valency_length()
graph.reflect_spec()
#graph.run_calc()
#graph.save_data()
graph.load_data()
#graph.plot_phase_diagrams()
graph.plot_logistic_regression()


## Radial distribution profile and Radii
obj = HandleRDPCaMKII()
obj.CG_valency_length()
obj.run(separate_hub_binding = True)
obj.save_figs()
#obj.save_hmws()
##




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



## Radial distribution profile and Radii
obj = HandleRDPCaMKII()
obj.CG_valency_length()
obj.run()
obj.save_figs()
obj.save_hmws()
##


'''

