import numpy as np
import pandas as pd

# reads quantum efficiency vector from data sheet
def read_qe(wavelength_range):
    q_df = pd.read_csv("quantum-efficiency-data/qe-sampled-points-interpolated.csv")
    lambda_min, lambda_max = wavelength_range
    q_temp = q_df.loc[
        (lambda_min <= q_df["Wavelength"]) & (q_df["Wavelength"] < lambda_max),
        "Quantum Efficiency",
    ]
    qe = pd.Series(q_temp.to_list(), index=range(*wavelength_range))
    return qe


# reformats quantum efficiency data into a vector where the value at each index i is the quantum efficiency at illumination i
def form_q_vec(illumination_df, qe):
    q_vec = qe[
        (illumination_df.bin_wavelength_min + illumination_df.bin_wavelength_max) // 2
    ]
    q_vec = np.array(q_vec.array)
    return q_vec


def fast_form_q_vec(
    qe,  # quantum efficiency
    bin_wavelength_range,  # length 2 ordered int tuple with first and last wavelengths detected
    bin_width,  # int denoting size of each wavelength bin
    N_ex,  # number of illuminations performed
):
    N_em = (bin_wavelength_range[1] - bin_wavelength_range[0]) // bin_width
    idx = np.arange(N_ex * N_em)
    bins_list = np.arange(*bin_wavelength_range, bin_width)
    bins_list_centered = bins_list + bin_width // 2
    q_vec = np.array(qe.loc[bins_list_centered])[idx % N_em]
    return q_vec


"""Computes Fisher information matrix of model, where
        A is the imaging matrix, which determines the mean number of photons that each pixel recieves
        x_vec is the vector of fluorophore concentrations
        q_vec is a series with quantum efficiencies at each illumination in A
    
    Uses the matrix formula to avoid for loops
"""


def FIM(A, x_vec, variance, q_vec = 1):
    y_vec = q_vec * (A @ x_vec)
    variance_vec = y_vec + variance
    inverse_variance_vec = 1 / variance_vec
    diag_vec = q_vec * q_vec * inverse_variance_vec * (1 + inverse_variance_vec / 2)
    F = (A.T * diag_vec) @ A
    return F


def FIM_poisson(A, x_vec):
    y = A @ x_vec
    y_inv = 1 / y
    y_inv[y_inv == np.inf] = 0
    return (A.T * y_inv) @ A
