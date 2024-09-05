
from bin.all1_edit_data import EditData
from bin.all3_plot_profile import PlotProfiles
from bin.all4_plot_profile_3d_ovito import Plot3dOvito
from bin.all5_plot_profile_3d_pyvista import Plot3dPyvista




obj = EditData()
obj.examples()
obj.run()

obj = PlotProfiles()
obj.examples()
obj.run()

obj = Plot3dOvito()
obj.examples()
obj.run()


obj = Plot3dPyvista()
obj.examples()
obj.run()


