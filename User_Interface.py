import streamlit as st

st.title("Hyperspectral microscopy")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from data_objects import Fluorophore, get_fluorophore_list, normalize
from imaging_model import fast_form_A
from information_matrix import read_qe, fast_form_q_vec, FIM

wavelength_range = st.slider(
    "Wavelength range (nm)", min_value=400, max_value=900, value=(400, 900), step=10
)
bin_width = st.slider(
    "Bin width (nm), must divide evenly into the wavelength range", 1, 20, 10
)

fluorophore_list = get_fluorophore_list(wavelength_range=wavelength_range)

tab1, tab2 = st.tabs(["FIM", "Fluorophores"])

names = map(Fluorophore.get_name, fluorophore_list)
name_list = list(names)
dict = {name_list[i]: fluorophore_list[i] for i in range(len(fluorophore_list))}
fluoro_name = tab2.selectbox("Choose a fluorophore", name_list, index=0)
fluorophore = dict[fluoro_name]
tab2.write("brightness = {}".format(fluorophore.brightness))

tab2.header("Fluorophore spectra")
fig, axs = plt.subplots(1, 1)
idx = np.arange(*wavelength_range)
normalized_excitation = normalize(fluorophore.spectra.excitation.copy())
axs.plot(idx, fluorophore.spectra.emission, color="green", label="emission")
axs.plot(idx, normalized_excitation, color="blue", label="excitation")
axs.fill_between(idx, fluorophore.spectra.emission, color="green", alpha=0.3)
axs.fill_between(idx, normalized_excitation, color="blue", alpha=0.3)
axs.set_xlabel("Wavelength (nm)")
axs.set_ylabel("")
axs.set_title("Normalized fluorophore spectra")
fig.legend()
tab2.pyplot(fig)


illumination_wavelengths = np.array(
    [405, 488, 561, 637]
)  # numpy array with the wavelength of each illumination
k = np.array(
    [1, 1, 1, 1]
)  # numpy array with (photon flux)*(voxel volume) for each illumination wavelength
bin_wavelength_range = wavelength_range  # length 2 ordered int tuple with first and last wavelengths detected
# bin_width = 10  # int denoting size of each wavelength bin

tab1.subheader("Experimental parameters")

# TODO: make this 2 column method work, issue may be due to the tabs
col1, col2 = tab1.columns(2)
with col1:
    tab1.markdown("#### Illumination")
    tab1.markdown("Wavelength of each illumination")
    df_illum = pd.DataFrame(illumination_wavelengths)
    df_illum_edited = tab1.experimental_data_editor(df_illum)
    illumination_wavelengths = np.array(df_illum_edited.squeeze())
    tab1.markdown("(photon flux)*(voxel volume) for each illumination")
    df_k_edited = tab1.experimental_data_editor(pd.DataFrame(k))
    k = np.array(df_k_edited.squeeze())

params = (
    illumination_wavelengths,
    k,
    bin_wavelength_range,
    bin_width,
    fluorophore_list,
)

A = fast_form_A(*params)
qe = read_qe(bin_wavelength_range)
q_vec = fast_form_q_vec(
    qe, bin_wavelength_range, bin_width, len(illumination_wavelengths)
)
qe = read_qe(bin_wavelength_range)
q_vec = fast_form_q_vec(
    qe, bin_wavelength_range, bin_width, len(illumination_wavelengths)
)

with col2:
    tab1.markdown("#### Sample")
    tab1.markdown("Concentration of fluorophores")
    x_vec = [1.0, 1.0, 1.0, 1.0]
    df_x_vec_edited = tab1.experimental_data_editor(pd.DataFrame(x_vec))
    x_vec = np.array(df_x_vec_edited.squeeze())
    variance = tab1.slider("Electronic Noise Variance", 1, 5, 2)

F = FIM(A, x_vec, q_vec, variance)
F_inv = np.linalg.inv(F)
CRLB = np.diagonal(F_inv)
results = pd.DataFrame(
    {
        "Fluorophore": list(map(Fluorophore.get_name, fluorophore_list)),
        "CRLB": CRLB,
        "FOM": x_vec / np.sqrt(CRLB),
    }
)

tab1.subheader("Detectibility results")
tab1.dataframe(results)
