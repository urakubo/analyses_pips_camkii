
from bin.all1_edit_data import EditData
from bin.all2_edit_connectivity_graph import EditConnectivityGraph
from bin.all3_plot_profile import PlotProfiles
from bin.all4_plot_profile_3d_ovito import Plot3dOvito
from bin.all5_plot_profile_3d_pyvista import Plot3dPyvista
from bin.all6_plot_profile_shared_PSD95 import PlotProfilesSharedPSD95
from bin.all7_plot_video_3d_ovito import MakeOvitoVideo

from bin.valency_length_CG3_plot_FRAP_video_3d_ovito import MakeOvitoVideoFRAP


'''
obj = EditData()
obj.inhibitor() #  conc_dependence(), valency_length(), valency_length_CG(), inhibitor()
obj.run()

obj = EditConnectivityGraph()
obj.boundary_conditions2() #  conc_dependence(), valency_length(), valency_length_CG()
obj.run()

obj = PlotProfiles()
obj.inhibitor()  #  conc_dependence(), valency_length(), valency_length_CG(), boundary_conditions2()
obj.run()

# I manually ran each one of them,
# because I do not know how to fully reset the ovito visualization system.
obj = Plot3dOvito()
obj.boundary_conditions2() #  conc_dependence(), valency_length(), valency_length_CG(), boundary_conditions2()
obj.inspect()
obj.run(1)

obj = Plot3dPyvista()
obj.boundary_conditions2() #  conc_dependence_merged(), valency_length(), valency_length_CG(), boundary_conditions2()
obj.run()

obj = PlotProfilesSharedPSD95()
obj.valency_length() #  conc_dependence(), valency_length(), valency_length_CG(), boundary_conditions2()
obj.run()


obj = MakeOvitoVideo()
obj.boundary_conditions2()
obj.inspect()
obj.run(3)
obj.make_a_video(3)
'''


target_molecule = 'CaMKII' # 'CaMKII', 'GluN2B', 'Both'
obj = MakeOvitoVideoFRAP()
#obj.valency_length_small_colony2()
#obj.inspect()
obj.valency_length_small_colony3()
i = 7*5+6 # val_12\R2_006
i = 7*5+2 # val_12\R2_006
obj.run(i, target_molecule)
obj.make_a_video(i)


