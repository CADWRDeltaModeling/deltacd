# Utility functions for DeltaCD
#
# This module includes tools to read and convert old file formats.

import numpy as np
import pandas as pd
import xarray as xr


def split_landuse(path_landuse, path_landuse_split_ratios):
    """ Split the land use information using the split ratios

        Parameters
        ----------
        path_landuse: str
        path_landuse_split_ratios: str

        Return
        ------
        pandas.DataFrame
    """
    df_landuse = pd.read_csv(path_landuse)
    df_split_ratios = pd.read_csv(path_landuse_split_ratios)

    col_crops_landuse = ~df_landuse.columns.isin(["DATE", "TYPE", "area_id"])
    col_crops_ratios = ~df_split_ratios.columns.isin(["area_id", "target_area_id", "zone", "ex"])

    df_landuse_split = df_landuse.copy()
    dfs_landuse_new = []
    area_id_max = df_landuse_split["area_id"].max()
    # Iterate through the split ratios. Each row will be a new area
    for _, row in df_split_ratios.iterrows():
        area_id = row["area_id"]
        mask_area_id = df_landuse_split["area_id"] == area_id
        ratio = row[col_crops_ratios]
        landuse_split = df_landuse.loc[mask_area_id, col_crops_landuse] * ratio.values
        df_landuse_split.loc[mask_area_id,
            col_crops_landuse] = df_landuse_split.loc[mask_area_id, col_crops_landuse] - landuse_split
        # Create land use data for the new area
        df_new_landuse = df_landuse[mask_area_id].copy()
        df_new_landuse.loc[:, col_crops_landuse] = landuse_split
        area_id_max += 1
        df_new_landuse['area_id'] = area_id_max
        dfs_landuse_new.append(df_new_landuse)
    df_landuse_new = pd.concat([df_landuse_split] + dfs_landuse_new, ignore_index=True)
    return df_landuse_new


def split_irrigation_efficiencies(path_irrigation_efficiencies,
                                  path_landuse_split_ratios):
    """ Split irrigation efficiencies using the splitting information

        Parameters
        ----------
        path_landuse: str
        path_landuse_split_ratios: str

        Return
        ------
        pandas.DataFrame
    """
    df_eta = pd.read_csv(path_irrigation_efficiencies)
    df_split_ratios = pd.read_csv(path_landuse_split_ratios)

    area_ids = df_split_ratios['area_id'].values
    df_eta_new = df_eta.loc[area_ids, :].reset_index()
    df_eta_new["area_id"] = np.arange(len(area_ids)) + len(df_eta) + 1
    # print(df_eta_new)
    df_eta_all = pd.concat([df_eta, df_eta_new.set_index("area_id")])
    return df_eta_all


def split_leach_water(path_leachwater, path_landuse_split_ratios):
    """ Split leach water information

        Parameters
        ----------
        path_leachwater: str
        path_landuse_split_ratios: str

        Return
        ------
        pandas.DataFrame
    """
    df_leach = pd.read_csv(path_leachwater)
    df_split_ratios = pd.read_csv(path_landuse_split_ratios)

    df_leach_new = df_leach.set_index('area_id').loc[df_split_ratios['area_id']].reset_index()
    # NOTE Set the values zero following the original DCD
    df_leach_new.iloc[:, 1:] = 0.
    df_leach_new["area_id"] = np.arange(len(df_leach_new)) + len(df_leach) + 1
    df_leach_all = pd.concat((df_leach, df_leach_new), ignore_index=True)
    return df_leach_all


def split_groundwater(path_groundwater, path_landuse_split_ratios):
    """ Split groundwater information

    P   arameters
        ----------
        path_groundwater: str
        path_landuse_split_ratios: str

        Return
        ------
        pandas.DataFrame
    """
    df_gw = pd.read_csv(path_groundwater)
    df_split_ratios = pd.read_csv(path_landuse_split_ratios)
    df_gw_new = df_gw.loc[:, [str(i) for i in df_split_ratios['area_id']]]
    df_gw_new.columns = np.arange(df_gw_new.shape[1]) + df_gw.shape[1]
    return pd.concat((df_gw, df_gw_new), axis=1)


def split_distribution(path_distribution, path_landuse_split_ratios):
    """ Split node distribution info

        Parameters
        ----------
        path_distribution: str
        path_landuse_split_ratios: str

        Return
        ------
        pandas.DataFrame
    """
    df_gw = pd.read_csv(path_distribution)
    df_split_ratios = pd.read_csv(path_landuse_split_ratios)
    return


def read_old_leach_water_file(path_leachwater):
    """ Read leach water information in the original DCD format
    """
    df_leach = pd.read_csv(path_leachwater, delim_whitespace=True,
                          header=None)
    df_leach.rename(columns={0: 'area_id'}, inplace=True)
    df_leach.rename(columns={i + 1: (i + 9) % 12 + 1 for i in range(12)}, inplace=True)
    return df_leach


def read_old_irrigation_efficiencies(path_irrigation_efficiencies) -> pd.DataFrame:
    """ Read irrigation efficiencies in the old DETAW-DCD format file and
        return a DataFrame

        Parameters
        ----------
        path_irrigation_efficiencies: str
            file path of an irrigation efficiencies in the old DETAW-DCD style.

        Return
        ------
        pandas.DataFrame
            irrigation efficiencies per area (island)
    """
    df_eta = pd.read_csv(path_irrigation_efficiencies,
                         header=None).reset_index().rename(
        columns={'index': 'area_id', 0: 'irrigation_efficiency'})
    df_eta.loc[:, 'area_id'] += 1
    df_eta.set_index('area_id', inplace=True)
    return df_eta


def read_subarea_info_file(path_file: str) -> pd.DataFrame:
    """ Read the subarea info text file and return it in a DataFrame

        The file has fixed format and not have headers.

        Parameters
        ----------
        path_file
            file path to read

        Returns
        -------
        pandas.DataFrame
    """
    colnames = ["area", "name", "uplow", "acreage", "docregion"]
    colspecs = ((0, 3), (5, 56), (57, 58), (60, 75), (81, 90))

    df = pd.read_fwf(path_file,
                     colspecs=colspecs,
                     skiprows=6,
                     names=colnames)
    return df
