
from bin.surface_tension0_plot_pyvista_CaMKII import PlotPyvistaCaMKII


## Target file definition

num_samples = 30
prefix, random_seed = '12_002', 1
prefix, random_seed = '12_006', 0
#prefix, random_seed = '12_005', 2


obj = PlotPyvistaCaMKII()
obj.CG_valency_length()
obj.plot_save( prefix = prefix, random_seed = random_seed, num_samples = num_samples )




