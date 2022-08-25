# SeparatedTopologies
Exploratory work on free energy calculations with GROMACS using "separated topologies" approach (Rocklin et al., 2013).

## Installation
```bash
# Create conda environment using environment.yml file
conda env create -f environment.yml

#Clone the github repository
git clone https://github.com/MobleyLab/SeparatedTopologies.git
```
The package currently requires an OpenEye license for oechem and oespruce.

## Usage
Example python script for setting up topology and coordinate files
 - for SepTop in the complex phase: `prep_SepTop_complex.py`
 - for SepTop in the solvent phase: `prep_SepTop_solvent.py`
 - for absolute hydration free energy calculations: `prep_solvent_absolute.py`

Example MDP files can be found in the folder `MDP_files`
