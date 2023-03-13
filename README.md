# hyperspectral
Finding the limit of the number fluorescent molecules we can distinguish

## Imaging model
We illuminate the four fluorophores with four wavelengths.
We generate the Fisher information matrix and calculate figures of merit to access how well each fluorophore can be imaged.

## Overview
The main scripts are `Model.ipynb` and `demo.ipynb`. View these scripts in action with by running
```
streamlit run User_Interface.py
```

### Directory structure
fluoro_spectra
- fluorophore excitation and emission spectra data in CSV files
- downloaded from [FPbase](https://www.fpbase.org/)

fluorophore-spectra
- repeat of `fluoro_spectra` with an added test file

FOM
- calculated figures of merit for each fluorophore
- generated from `Model.ipynb`

Figures 
- variation of plots generated from `plotting/*` and `Model.ipynb`

playground
- scripts for temporary testing and comparison

quantum-efficiency-data
- detector quantum efficiency data in CSV files
- ORCA-Fusion BT camera
- script generating plots of the data

### Dependencies
- numpy
- pandas
- matplotlib
- notebook/ipython
- streamlit