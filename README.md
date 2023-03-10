# hyperspectral
Finding the limit of the number fluorescent molecules we can distinguish

## Imaging model
We illuminate the four fluorophores with four wavelengths.
We generate the Fisher information matrix and calculate figures of merit to access how well each fluorophore can be imaged.

## Developers guide
The main script is `Model.ipynb`.

### Directory structure
fluoro_spectra
- fluorophore excitation and emission spectra data in CSV files
<!-- - downloaded from FPbase -->
FOM
- calculated figures of merit for each fluorophore
- generated from `Model.ipynb`

Figures 
- variation of plots generated from `plotting/*` and `Model.ipynb`

playground
- scripts for temporary testing and comparison

### Dependencies
- numpy
- pandas
- matplotlib
- notebook/ipython