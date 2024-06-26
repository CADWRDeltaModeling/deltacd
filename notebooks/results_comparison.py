import pandas as pd
import geopandas as gpd
# from deltacd.utils.dss_to_df import dss_to_df
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import pylab
import pyhecdss
import os
import argparse
import yaml
import numpy as np

def create_argparser():
    """ Create an argument parser.

        Returns
        -------
        argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("model_name", type=str,
                        help="A str that refers to the model name e.g. dsm2")
    parser.add_argument("input_yaml", type=str,
                        help="An input YAML file to provide path to old and new filenames")

    return parser

def plot_diff(diff_desc_df, title, fig_fname=None):
    plt.figure(figsize=(10,4))
    plt.title(title)
    plt.bar(diff_desc_df.index,diff_desc_df['mean'])
    plt.ylabel('difference (cfs)')
    if fig_fname[:3] != 'cal':
        plt.xticks(diff_desc_df.index[::5],  rotation='vertical')
    plt.savefig(fig_fname)
    plt.show()


def compare_results(filenames):

    path_new_filename = filenames.get('path_new_filename')
    path_old_div_csv = filenames.get('path_old_div_csv')
    path_old_drain_csv = filenames.get('path_old_drain_csv')
    path_old_seep_csv = filenames.get('path_old_seep_csv')
    nodes_file = filenames.get('nodes_filename')
    outfile_prefix = filenames.get('outfile_prefix')

    # Read the new output file
    ds = xr.open_dataset(path_new_filename)
    dcd_new_df = ds.to_dataframe()

    # Make seep df
    new_drain_df = pd.pivot_table(dcd_new_df, values = 'drainage', index = 'time', columns = 'node')
    new_drain_df.index.name=None
    new_drain_df.index.name = "datetime"
    new_drain_df.columns.name = None

    # Make seep df
    new_seep_df = pd.pivot_table(dcd_new_df, values = 'seepage', index = 'time', columns = 'node')
    new_seep_df.index.name=None
    new_seep_df.index.name = "datetime"
    new_seep_df.columns.name = None

    # Make div df
    new_div_df = pd.pivot_table(dcd_new_df, values = 'diversion', index = 'time', columns = 'node')
    new_div_df.index.name=None
    new_div_df.index.name = "datetime"
    new_div_df.columns.name = None

    if outfile_prefix == 'calsim':
        new_div_m_df = new_div_df.groupby(pd.Grouper(freq='ME')).mean()
        new_drain_m_df = new_drain_df.groupby(pd.Grouper(freq='ME')).mean()
        new_seep_m_df = new_seep_df.groupby(pd.Grouper(freq='ME')).mean()

    # Read dsm2 nodes to sort by longitude. West to east.
    gdf = gpd.read_file(nodes_file).to_crs(epsg=26910)
    # Extract longitude from the geometry
    gdf['longitude'] = gdf.geometry.x

    # Sort by longitude
    gdf_sorted = gdf.sort_values('longitude')

    # Exclude 'datetime' and identify columns to be reordered
    csv_columns = new_div_df.columns.tolist()
    columns_to_sort = [col for col in csv_columns if col != 'datetime']

    # Check if these columns match with GeoDataFrame ids
    sorted_columns = gdf_sorted['id'].astype(str).tolist()

    # Ensure all sorted columns are present in the CSV
    sorted_columns = [col for col in sorted_columns if col in columns_to_sort]

    # Read old div output
    old_div_df = pd.read_csv(path_old_div_csv,parse_dates=[0])
    old_div_df = old_div_df.set_index(old_div_df.columns[0])

    # Subtract new and old div
    if outfile_prefix == 'calsim':
        div_diff_df = new_div_m_df - old_div_df
    else:
        div_diff_df = new_div_df - old_div_df

    # Compute descriptive stats sorted by longitude
    if outfile_prefix == 'calsim':
        diff_desc = div_diff_df.replace([np.inf, -np.inf],np.nan).describe()
        diff_desc = diff_desc.T
    else:
        diff_desc = div_diff_df[sorted_columns].replace([np.inf, -np.inf],np.nan).describe()
        diff_desc = diff_desc.T

    # Plot the diff
    title = 'DIV-FLOW difference between deltaCD and dcdv1.3 outputs'
    fig_fname = outfile_prefix + '_div_diff.png'
    plot_diff(diff_desc_df=diff_desc, title=title, fig_fname=fig_fname)

    # Read old drain output
    old_drain_df = pd.read_csv(path_old_drain_csv,parse_dates=[0])
    old_drain_df = old_drain_df.set_index(old_drain_df.columns[0])

    # Subtract new and old drain
    if outfile_prefix == 'calsim':
        drain_diff_df = new_drain_m_df - old_drain_df
    else:
        drain_diff_df = new_drain_df - old_drain_df

    # Compute descriptive stats sorted by longitude
    if outfile_prefix == 'calsim':
        diff_desc = drain_diff_df.replace([np.inf, -np.inf],np.nan).describe()
        diff_desc = diff_desc.T
    else:
        diff_desc = drain_diff_df[sorted_columns].replace([np.inf, -np.inf],np.nan).describe()
        diff_desc = diff_desc.T

    # Plot the diff
    title = 'DRAIN-FLOW difference between deltaCD and dcdv1.3 outputs'
    fig_fname = outfile_prefix + '_drain_diff.png'
    plot_diff(diff_desc_df=diff_desc, title=title, fig_fname=fig_fname)

    # Read old seep output
    old_seep_df = pd.read_csv(path_old_seep_csv,parse_dates=[0])
    old_seep_df = old_seep_df.set_index(old_seep_df.columns[0])

    # Subtract new and old seep
    if outfile_prefix == 'calsim':
        seep_diff_df = new_seep_m_df - old_seep_df
    else:
        seep_diff_df = new_seep_df - old_seep_df

    # Compute descriptive stats sorted by longitude
    if outfile_prefix == 'calsim':
        diff_desc = seep_diff_df.replace([np.inf, -np.inf],np.nan).describe()
        diff_desc = diff_desc.T
    else:
        diff_desc = seep_diff_df[sorted_columns].replace([np.inf, -np.inf],np.nan).describe()
        diff_desc = diff_desc.T

    # Plot the diff
    title = 'SEEP-FLOW difference between deltaCD and dcdv1.3 outputs'
    fig_fname = outfile_prefix + '_seep_diff.png'
    plot_diff(diff_desc_df=diff_desc, title=title, fig_fname=fig_fname)

    ds.close()

if __name__ == "__main__":
    parser = create_argparser()
    args = parser.parse_args()

    model_name = args.model_name
    fname_main_yaml = args.input_yaml

    with open(fname_main_yaml, 'r') as file_in:
        filenames_to_compare = yaml.safe_load(file_in)

    compare_results(filenames_to_compare.get(model_name))

