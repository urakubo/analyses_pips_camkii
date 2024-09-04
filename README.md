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

Analyses were conducted through two steps: the conversion of lammpstrj files into intermediate data files and the visualization based on the intermediate data. In the methods of “SpecDatasets” class, the following instance variables should be specified:

| instance variable | variable type | content |
| -------- | -------- | -------- |
| **`self.dir_lammpstrj`** | str | Directory for target lammpstrj files. |
| **`self.dir_edited_data`**| str |Directory for intermediate data files. |
| **`self.dir_imgs_root`**| str | Directory for image files. |
| **`self.filenames_lammpstrj`**| list/tuple | Filenames of target lammpstrj files. |
| **`self.filenames_edited`**| list/tuple | Filenames of intermediate data files. |



### specification_datasets.py
Specifies the directories and filenames of lammpstrj and intermediate data files.  


Under documentation...
