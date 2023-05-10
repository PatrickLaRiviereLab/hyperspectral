# defines "Fluorophore" object to handle fluorophore data
# creates list of Fluorophore objects for use in other code

import numpy as np
import pandas as pd

fluorophore_string_list = ["mCherry", "mEmerald", "mNeptune2.5", "mTagBFP2"]
brightness_list = [15.85, 39.1, 22.8, 32.38]

# returns a normalized version of the series
def normalize(series):
    area = np.trapz(series, dx=1)
    normalized = series / area
    return normalized


# sets the wavelength column of a fluorophore dataframe to be range(a,b)
def set_wavelength_range(df, a, b):
    df_new = pd.DataFrame(
        {
            "wavelength": range(a, b),
            "excitation": [0.0] * (b - a),
            "emission": [0.0] * (b - a),
        },
        index=range(a, b),
    )
    wavelength_list = list(df.wavelength)
    for wavelength in range(a, b):
        if wavelength in wavelength_list:
            boolean_index = df.wavelength == wavelength
            excitation = df.loc[boolean_index, "excitation"].iloc[0]
            emission = df.loc[boolean_index, "emission"].iloc[0]
            df_new.loc[wavelength:wavelength, ("excitation", "emission")] = [
                excitation,
                emission,
            ]
    return df_new


# reads and processes fluorophore data file
def load_fluorophore(
    name,
    filename=None,
    excitation_column_name=None,
    emission_column_name=None,
    index_range=(400, 900),
):
    if filename == None:
        filename = (
            "fluorophore-spectra/" + name.replace(".", "p") + "_fpbase_spectra.csv"
        )

    df = pd.read_csv(filename)
    df.fillna(0, inplace=True)

    if excitation_column_name == None:
        excitation_column_name = name + " ex"
    if emission_column_name == None:
        emission_column_name = name + " em"

    df[emission_column_name] = normalize(df[emission_column_name])
    #df[excitation_column_name] = normalize(df[excitation_column_name])

    df_columns = df.rename(
        columns={name + " ex": "excitation", name + " em": "emission"}
    )
    df_new = set_wavelength_range(df_columns, index_range[0], index_range[1])

    return df_new


# bundles together information about a single fluorophore
class Fluorophore:
    def __init__(
        self,
        name,
        brightness,
        filename=None,
        excitation_column_name=None,
        emission_column_name=None,
        index_range=(400, 900),
    ):
        self.spectra = load_fluorophore(
            name,
            filename=filename,
            excitation_column_name=excitation_column_name,
            emission_column_name=emission_column_name,
            index_range=index_range,
        )
        self.name = name
        self.brightness = brightness

    def get_spectra(self):
        return self.spectra

    def get_name(self):
        return self.name


# returns list of all spectra
def get_spectra_list(wavelength_range=(400, 900)):
    spectra_list = []
    for fluorophore_name in fluorophore_string_list:
        # print(fluorophore_name)
        spectra_list.append(
            load_fluorophore(fluorophore_name, index_range=wavelength_range)
        )
    return spectra_list


# returns list of all fluorophore objects
## preferable to above function (if you want to call both, just call this one and use map(Fluorophore.get_spectra,fluorophore_list)
## to compute spectra_list in order to reduce redundant computation)
def get_fluorophore_list(wavelength_range=(400, 900)):
    fluorophore_list = []
    for fluorophore_name, fluorophore_brightness in zip(
        fluorophore_string_list, brightness_list
    ):
        fluorophore_list.append(
            Fluorophore(
                fluorophore_name, fluorophore_brightness, index_range=wavelength_range
            )
        )
    return fluorophore_list
