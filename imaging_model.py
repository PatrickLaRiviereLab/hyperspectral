# functions for computing entries of imaging matrix (matrix A) and the matrix itself
# organizes data for illumination i into the rows of a dataframe
# depends on Fluorophore object in data_objects.py, even though data_objects is not imported

import numpy as np
import pandas as pd

""" Return a DataFrame where row i contains data to calculate matrix entry A_ij
"""


def form_illumination_df(
    illumination_data,  # dataframe with columns:
    # "wavelength" for wavelength of illumination
    # "k" for the value (photon flux)*(volume of voxel)
    bin_wavelength_range,  # pair with first and last detected wavelengths
    bin_width,  # width of each bin
):
    illumination_df = pd.DataFrame(
        {
            "illumination_wavelength": [],
            "bin_wavelength_min": [],
            "bin_wavelength_max": [],
            "k": [],
        }
    )
    bins_list = range(*bin_wavelength_range, bin_width)
    counter = 0
    for illumination_index in illumination_data.index:
        illumination_wavelength = illumination_data.wavelength[illumination_index]
        k = illumination_data.k[illumination_index]
        for bin_wavelength_min in bins_list:
            bin_wavelength_max = bin_wavelength_min + bin_width

            illumination_df.loc[counter] = {
                "illumination_wavelength": illumination_wavelength,
                "bin_wavelength_min": bin_wavelength_min,
                "bin_wavelength_max": bin_wavelength_max,
                "k": k,
            }
            counter += 1

    return illumination_df


def calc_emission_area(fluorophore, lambda_min, lambda_max):
    area = np.trapz(fluorophore.spectra.emission.loc[lambda_min : lambda_max - 1])
    return area


def calc_A_entry(illumination_df_row, fluorophore):
    k = illumination_df_row.k
    b = fluorophore.brightness
    excitation = fluorophore.spectra.excitation[
        illumination_df_row.illumination_wavelength
    ]
    emission = calc_emission_area(
        fluorophore,
        illumination_df_row.bin_wavelength_min,
        illumination_df_row.bin_wavelength_max,
    )

    return k * b * excitation * emission


def form_A(illumination_df, fluorophore_list):
    n = len(illumination_df.index)
    m = len(fluorophore_list)
    A = np.zeros((n, m))
    for i in range(n):
        for j in range(m):
            A[i, j] = calc_A_entry(illumination_df.loc[i], fluorophore_list[j])

    return A


# computes imaging matrix directly from parameters as quickly as possible
def fast_form_A(
    illumination_wavelengths,  # numpy array with the wavelength of each illumination
    k,  # numpy array with (photon flux)*(voxel volume) for each illumination wavelength
    bin_wavelength_range,  # length 2 ordered int tuple with first and last wavelengths detected
    bin_width,  # int denoting size of each wavelength bin
    fluorophore_list,  # list of fluorophores in image
):
    try:
        assert illumination_wavelengths.size == k.size
    except:
        raise ValueError(
            "arguments 'illumination_wavelengths' and 'k' must be the same size"
        )

    try:
        assert bin_wavelength_range[0] == int(bin_wavelength_range[0])
        assert bin_wavelength_range[1] == int(bin_wavelength_range[1])
        assert bin_width == int(bin_width)
    except:
        raise ValueError("wavelengths must be given as integers")

    try:
        assert (bin_wavelength_range[1] - bin_wavelength_range[0]) % bin_width == 0
    except:
        raise ValueError("size of wavelength range must be divisible by bin width")

    N_ex = illumination_wavelengths.size
    N_em = (bin_wavelength_range[1] - bin_wavelength_range[0]) // bin_width
    idx = np.arange(N_ex * N_em)
    excitation_index = idx // N_em
    emission_index = idx % N_em
    k = k[excitation_index]

    columns = []
    for fluorophore in fluorophore_list:
        column_excitation = np.array(
            fluorophore.spectra.excitation[illumination_wavelengths]
        )[excitation_index]
        fluorophore_emission_array = np.array(
            fluorophore.spectra.emission.loc[
                bin_wavelength_range[0] : bin_wavelength_range[1] - 1
            ]
        )
        column_emission = np.trapz(fluorophore_emission_array.reshape(N_em, bin_width))[
            emission_index
        ]

        column = k * fluorophore.brightness * column_excitation * column_emission
        columns.append(column)

    A = np.array(columns).T

    return A


# rescales the imaging matrix so that the number of photons detected are the same 
def fast_form_A_from_photons(
        desired_photons,  # float with the desired number of photons from each illumination
        illumination_wavelengths,  # numpy array with the wavelength of each illumination
        bin_wavelength_range,  # length 2 ordered int tuple with first and last wavelengths detected
        bin_width,  # int denoting size of each wavelength bin
        fluorophore_list,  # list of fluorophores in image
    ):
    N_ex = illumination_wavelengths.size
    N_em = (bin_wavelength_range[1] - bin_wavelength_range[0]) // bin_width
    A = fast_form_A(illumination_wavelengths, np.ones(illumination_wavelengths.shape),
                    bin_wavelength_range, bin_width, fluorophore_list)
    row_sums = np.sum(A, axis=1)
    illumination_sums = row_sums.reshape(N_ex, N_em).sum(axis=1)
    k = desired_photons / illumination_sums
    new_A = fast_form_A(illumination_wavelengths, k,
                    bin_wavelength_range, bin_width, fluorophore_list)
    return new_A