# PeptoGrid

PeptoGrid is a set of tools implementing rescoring function for AutoDock Vina scoring predictions

## Related Publication

[Zalevsky, Zlobin, Gedzun, Reshetnikov, Lovat, Malyshev, Doronin, Babkin, Golovin, "PeptoGrid—Rescoring Function for AutoDock Vina to Identify New Bioactive Molecules from Short Peptide Libraries"](https://doi.org/10.3390/molecules24020277)

## Install 

### Prerequisites

PeptoGrid written in Python2.7 and depends on several other projects:

[NumPy](http://www.numpy.org/) for calculations

[H5PY](https://www.h5py.org/) for efficient data storage

[ODDT](https://github.com/oddt/oddt) for processing structural data

[OpenBabel](http://openbabel.org/wiki/Main_Page) as a backend for ODDT

[PyMOL](https://pymol.org/2/) for visualization

Python dependencies can be installed with pip:

> pip install numpy h5py oddt

OpenBabel should be installed with your package manager or with any other suitable [way](http://openbabel.org/wiki/Category:Installationhttp://openbabel.org/wiki/Category:Installation)

For Ubuntu/Debian users:

> apt-get install libopenbabel4v5 python-openbabel

PyMOL can be download from official [site](https://pymol.org/2/) or installed with package manager.

For Ubuntu/Debian users:
> apt-get install pymol

### Installation

PeptoGrid itself doesn't require any installation. Just clone or download repository and it's ready to go

> git clone https://github.com/aozalevsky/peptogrid.git

## Usage

No preprocessing (splitting or conversion to PDB) of AutoDock Vina results is required. PeptoGrid works with raw pdbqt outputs.

To build atom frequency grid:
> python peptogrid/grid.py -m grid.hdf5 -f *out.pdbqt -c vina.cfg -s 0.1

To rescore poses:
> python peptogrid/score.py -m grid.hdf5 -f *out.pdbqt -c vina.cfg -s 0.1 -o score_table.csv

To visualize in PyMOL:
* run PyMOL
* in command line type in following command:
> run peptogrid/load_grid.py
* load calculated grid:
> load_grid grid.hdf5

For details about volumetric data in PyMOL refer to [PyMOL Wiki](https://pymolwiki.org/index.php/Volume)
