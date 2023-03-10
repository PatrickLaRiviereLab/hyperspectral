import pandas as pd
import numpy as np


def generate_normed_spectra():
    df_mEmerald = pd.read_csv("fluoro_spectra/mEmerald_fpbase_spectra.csv")
    df_mEmerald.fillna(0, inplace=True)
    df_mTagBFP2 = pd.read_csv("fluoro_spectra/mTagBFP2_fpbase_spectra.csv")
    df_mTagBFP2.fillna(0, inplace=True)
    df_mCherry = pd.read_csv("fluoro_spectra/mCherry_fpbase_spectra.csv")
    df_mCherry.fillna(0, inplace=True)
    df_mNeptune2p5 = pd.read_csv("fluoro_spectra/mNeptune2p5_fbpase_spectra.csv")
    df_mNeptune2p5.fillna(0, inplace=True)

    mEmerald_columns = df_mEmerald.rename(
        columns={"mEmerald ex": "excitation", "mEmerald em": "emission"}
    )
    mTagBFP2_columns = df_mTagBFP2.rename(
        columns={"mTagBFP2 ex": "excitation", "mTagBFP2 em": "emission"}
    )
    mCherry_columns = df_mCherry.rename(
        columns={"mCherry ex": "excitation", "mCherry em": "emission"}
    )
    mNeptune2p5_columns = df_mNeptune2p5.rename(
        columns={"mNeptune2.5 ex": "excitation", "mNeptune2.5 em": "emission"}
    )

    spectra_list = [
        mEmerald_columns,
        mTagBFP2_columns,
        mCherry_columns,
        mNeptune2p5_columns,
    ]
    for j in range(0, 4):
        # print("j = ", j)
        excitation_column = spectra_list[j]["excitation"]
        area = np.trapz(excitation_column, dx=1)
        # print("Initial excitation area: ", area)
        normalized = excitation_column / area
        excitation_column = normalized
        spectra_list[j]["excitation"] = normalized
        new_area = np.trapz(excitation_column, dx=1)
        # print("New excitation area: ", new_area)
        column = spectra_list[j]["emission"]
        area_emission = np.trapz(column, dx=1)
        # print("Initial emission area: ", area_emission)
        normalized_em = column / area_emission
        column = normalized_em
        spectra_list[j]["emission"] = normalized_em
        new_area_em = np.trapz(column, dx=1)
        # print("New emission area: ", new_area_em)
    # return spectra_list
    return [mEmerald_columns, mTagBFP2_columns, mCherry_columns, mNeptune2p5_columns]


def append_zeros_to_spectra(
    mEmerald_columns, mTagBFP2_columns, mCherry_columns, mNeptune2p5_columns
):
    wavelength = 701
    for index in range(401, 601):
        mEmerald_columns.loc[len(mEmerald_columns.index)] = [wavelength, 0.0, 0.0]
        wavelength = wavelength + 1

    wavelength = 651
    for index in range(341, 591):
        mTagBFP2_columns.loc[len(mTagBFP2_columns.index)] = [wavelength, 0.0, 0.0]
        wavelength = wavelength + 1

    wavelength = 851
    for index in range(621, 671):
        mNeptune2p5_columns.loc[len(mNeptune2p5_columns.index)] = [wavelength, 0.0, 0.0]
        wavelength = wavelength + 1
    spectra_list = [
        mEmerald_columns,
        mTagBFP2_columns,
        mCherry_columns,
        mNeptune2p5_columns,
    ]
    return spectra_list
