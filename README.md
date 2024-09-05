![System requirements](https://img.shields.io/badge/python-3.8-red.svg)
![System requirements](https://img.shields.io/badge/platform-win%2064,%20linux%2064-green.svg)
[![License: GPL v3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# Analysis code for "the role of protein shape in multiphasic separation within condensates"
The GitHub repository contains analysis programs for the study "the role of protein shape in multiphasic separation within condensates" by Pandey et al. [1].
All programs were written in Python3.8 (Windows) and designed for the analyses of output from LASSI simulation engine [2].

[1] https://www.biorxiv.org/content/10.1101/2024.08.26.606306v1

[2] https://github.com/Pappulab/LASSI


| directory | contents |
| -------- | -------- |
| **`bin`** |Executable programs. |
| **`lib`**| Shared libraries. |
| **`workspace`**| Accessary programs. |


### specification_datasets.py

Analyses are conducted through two steps: the conversion of lammpstrj files into intermediate data files and the visualization based on the intermediate data. Their directories and filenames are specified as instance variables of a superclass "SpecDatasets" in "specification_datasets.py." In the methods of "SpecDatasets", the instance variables below should be defined. They are referred by executable programs (bin/).

| instance variable | variable type | content |
| -------- | -------- | -------- |
| **`self.dir_lammpstrj`** | str | Directory for target lammpstrj files. |
| **`self.dir_edited_data`**| str |Directory for intermediate data files. |
| **`self.dir_imgs_root`**| str | Directory for image files. |
| **`self.filenames_lammpstrj`**| list/tuple | Filenames of target lammpstrj files. |
| **`self.filenames_edited`**| list/tuple | Filenames of intermediate data files. |
| **` self.filename_lammpstrj_matrix `**| func(v, l)  | Filename of a target lammpstrj file specified by valency/GluN2B conc (v) and length/STG conc (l). It is used to draw matrix graphs and phase planes. |
| **`self.filename_edited_matrix`**| func(v, l) | Filename of the intermediate data file specified by valency/GluN2B conc (v) and length/STG conc (l). It is used to draw matrix graphs and phase planes. |



### Some .py files

| .py file | functions |
| -------- | -------- |
| **`lib/paramters.py`** | It defines basic parameters such as “ID versus molecular name“, “the size of lattice space” and so on.  The defined variables (e.g., space) are referred such as "p.space" in executable programs (bin/). |
| **`lib/colormap.py`**| It defines colors. The color universal design was utilized for color blindness. |
| **` control_board_*.py `**| Workspaces. They would be edited depending on further simulation and analyses. |


## Example

#### Programs and data

| diretory / file | contents |
| -------- | -------- |
| **`specification_datasets.py`** | The method "examples" is referred. |
| **`control_board_example.py`** | Control board that calls executables. |
| **`example_data/lammpstrj`**| CG.zip, PS.zip, CGPS.zip. Please unzip them. |
| **`example_data/edited`**| Empty. |
| **`example_data/imgs`**| Empty. |

#### control_board_example.py

It calls the following classes.

- EditData
- PlotProfiles
- Plot3dOvito
- Plot3dPyvista

#### EditData

The 'PlotProfiles' class converts the lammpstrj files ('CG.lammpstrj', 'PS.lammpstrj', and 'CGPS.lammpstrj' in example_data/lammpstrj) into intermediate data files ('CG_.pickle', 'PS_.pickle', and 'CGPS_.pickle'; example_data/edited). Each pickle file has a dict variable with the following keys:

| key | value type | description |
| -------- | -------- | -------- |
| **`mc_step`** | int | Sampled MC step. Final MC step in general. |
| **`sampling_frame`** | int | Sampled time frame. Final time frame in general. |
| **`region_condensate_in_grid_mesh`** | dict[Z] | Condensate region of X. Each condensate region is defined by the region over the half maximal levels of blurred Y. |
| **`conc_condensate`** | dict[Y][X] | Conc of X in the condensate region Y. |
| **`conc_periphery`** | dict[X] | Conc of X in the periphery region. The periphery region is defined by the region over d/2 distant from the center in the lattice space [d, d, d]. |
| **` concs_in_grid_mesh`** | dict[X] | Conc of X in the 3d space (3d np.array, float). Grid locations of beads were blurred by the gaussian (sd: p.sigma in lattice unit). |
| **` locs_in_grid_mesh`** | dict[X] | Conc of X in the 3d space (3d np.array, int). Grid locations of beads. |
| **`rdf_bins`**| 1d np.array, int | Bins for radial distribution profile (RDP) (in lattice unit).  |
| **`rdf_sampling_frames`**| list, int | Sampled frames for RDP |
| **`rdf `**| dict[X] | Radial distribution profile X (2d np.array, float)  |
| **`dir_lammpstrj `**| str | dir_lammpstrj. |
| **`filename_lammpstrj`**| str | filename_lammpstrj. |
| **`rotated_region_condensate_in_grid_mesh`** | dict[X] | Same as **` region_condensate_in_grid_mesh `** but rotated. |
| **`rotated_concs_in_grid_mesh`** | dict[X] | Same as **`concs_in_grid_mesh`** but rotated. |

Here, X, Y ∈ ['All', 'CaMKII', 'GluN2B', 'STG', 'PSD95'], Z ∈ ['All', 'CaMKII', 'GluN2B', 'STG', 'PSD95', 'dilute'].

#### PlotProfiles

It produces graphs that show the intensity levels at the section that divided the center of the condensate as well as the RDPs from condensate center-of-mass ('CG.svg', 'PS.svg', and 'CGPS.svg' in example_data/imgs/profiles).


#### Plot3dOvito

It visualizes protein bead distributions in the 3D space
('CG.svg', 'PS.svg', and 'CGPS.svg' in example_data/imgs/profiles_3d_ovito).


#### Plot3dPyvista

It visualizes the 3D shapes of the condensates using their rendered volumes
('CG.svg', 'PS.svg', and 'CGPS.svg' in example_data/imgs/profiles_3d_pyvista).

