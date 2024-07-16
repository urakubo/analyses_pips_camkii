
from bin.surface_tension0_plot_pyvista_CaMKII_tmp import PlotPyvistaCaMKII

from bin.surface_tension1_calc_radius_CaMKII_condensate import HandleRDPCaMKII
from bin.surface_tension2_plot import PlotSurfaceTension
from bin.surface_tension3_plot_phase_diagram import PlotSurfaceTensionPhaseDiagramValencyLength


from bin.connectivity_plot_profile_3d_ovito import Plot3dOvitoConnectivity


## Example plot

num_samples = 10
prefix, random_seed = '12_002', 1

num_samples = 6
prefix, random_seed = '12_006', 8


obj = PlotPyvistaCaMKII()
obj.CG_valency_length()
obj.plot_save( prefix = prefix, random_seed = random_seed, num_samples = num_samples )
##



'''

obj = Plot3dOvitoConnectivity()
obj.CG_valency_length()
#obj.inspect()
#obj.plot_an_image(37, mode = 'CaMKII_hub_beads') # 37, 41
obj.plot_an_image(38, mode = 'CaMKII_hub_beads') # 37, 41



## Radial distribution profile and Radii
obj = HandleRDPCaMKII()
obj.CG_valency_length()
obj.run(separate_hub_binding = True)
obj.save_figs()
#obj.save_hmws()
##


# Compare two representative case: spatial arrangement.
obj = PlotSurfaceTension()
obj.CG_valency_length()
obj.apply_specification()
#obj.inspect_targets()
obj.multiple_run_plot(['12_001','12_002','12_003','12_004','12_005','12_006'])
#obj.show_polar_graphs(filenames_prefix = ['12_001','12_002','12_003','12_004','12_005','12_006'])
#obj.show_polar_graphs(filenames_prefix = ['12_006'])



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

