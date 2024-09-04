![System requirements](https://img.shields.io/badge/python-3.8-red.svg)
![System requirements](https://img.shields.io/badge/platform-win%2064,%20linux%2064-green.svg)
[![License: GPL v3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# Analysis code for “the role of protein shape in multiphasic separation within condensates”
The GitHub space contains analysis programs for the study “the role of protein shape in multiphasic separation within condensates” by Vikas Pandey,  Tomohisa Hosokawa,  Yasunori Hayashi,  Hidetoshi Urakubo [1].
All programs were written in Python3.8 (Windows) and designed for the analyses of output from LASSI simulation engine [2].

[1] https://www.biorxiv.org/content/10.1101/2024.08.26.606306v1

[2] https://github.com/Pappulab/LASSI


| directory | contents |
| -------- | -------- |
| **`bin`** |Executable programs. |
| **`lib`**| Shared libraries. |
| **`workspace`**| Accessary programs. |


### specification_datasets.py

Analyses are conducted through two steps: the conversion of lammpstrj files into intermediate data files and the visualization based on the intermediate data. The "SpecDatasets" is a superclass of executable programs (bin/). In the methods of "SpecDatasets" class, the following instance variables must be defined:

| instance variable | variable type | content |
| -------- | -------- | -------- |
| **`self.dir_lammpstrj`** | str | Directory for target lammpstrj files. |
| **`self.dir_edited_data`**| str |Directory for intermediate data files. |
| **`self.dir_imgs_root`**| str | Directory for image files. |
| **`self.filenames_lammpstrj`**| list/tuple | Filenames of target lammpstrj files. |
| **`self.filenames_edited`**| list/tuple | Filenames of intermediate data files. |
| **` self.filename_lammpstrj_matrix `**| func(v, l)  | Filename of a target lammpstrj file specified by valency/GluN2B conc (v) and length/STG conc (l). It is used to draw matrix graphs and phase planes. |
| **`self.filename_edited_matrix`**| func(v, l) | Filename of the intermediate data file specified by valency/GluN2B conc (v) and length/STG conc (l). It is used to draw matrix graphs and phase planes. |

### Some key .py files

| .py file | functions |
| -------- | -------- |
| **`lib/paramters.py`** | It defines basic parameters such as “ID versus molecular name“, “the size of lattice space” and so on.  The defined variables are referred such as "p.space". |
| **`lib/colormap.py`**| It defines colors. The color universal design was utilized for color-blindness. |
| **` control_board_*.py `**| Workspaces. They would be edited depending on further simulation and analyses. |


## control_board_example.py


