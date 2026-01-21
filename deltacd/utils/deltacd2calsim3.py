import xarray as xr
import argparse
import pandas as pd
import pyhecdss
import yaml

def create_argparser():
    """ Create an argument parser.

        Returns
        -------
        argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("config_file", type=str,
                        help="A str path that refers to the YAML config file")

    return parser

def daily_to_monthly(daily_netcdf, monthly_dss, config):
    """ Convert daily netCDF data to monthly netCDF data.

        Parameters
        ----------
        daily_netcdf : str
            Path to the daily netCDF file.
        monthly_dss : str
            Path to the output monthly DSS file.
        config : dict
            Configuration dictionary from YAML.
    """
    # Open the daily netCDF file
    ds_daily = xr.open_dataset(daily_netcdf)
    ds_dims = ds_daily.dims
    df_daily=ds_daily.to_dataframe().reset_index()

    rename_dict = config['rename_dict']
    df_daily = df_daily.rename(columns=rename_dict)
    print("Converting daily data to monthly data...")

    rename_vars = list(rename_dict.values())
    # Keep only the dimension columns and renamed variables. Limit to time and only variables that yaml specified.
    cols_to_keep = list(set(['time'] + list(ds_daily.sizes.keys()) + rename_vars))
    cols_to_keep = [col for col in cols_to_keep if col in df_daily.columns]
    df_daily = df_daily[cols_to_keep]

    with pyhecdss.DSSFile(str(monthly_dss), create_new=True) as dcddss:
        # Get ptype and units from config
        ptype = config['ptype']
        units = config['units']
        # assign the parts A, E, F of the dss file from config
        A = config['A']
        E = config['E']
        F = config['F']

        # Loop through unique nodes or area_ids. This was added to handle both node-based and area-based data.
        if 'node' in ds_dims:
            for node in df_daily["node"].unique():
                df_daily_node = df_daily[df_daily["node"] == node].drop(columns=["node"])
                df_daily_node.set_index("time", inplace=True)
                df_monthly_node = df_daily_node.resample('ME').mean().to_period('M')
                for var in rename_vars:
                    if var in df_monthly_node.columns:
                        path = f"/{A}/{node}/{var}//{E}/{F}/"
                        dcddss.write_rts(path, df_monthly_node[[var]], units, ptype)
        elif 'area_id' in ds_dims:
            for area_id in df_daily["area_id"].unique():
                # Skip area_ids that contain underscore
                if '_' in str(area_id):
                    continue
                df_daily_area = df_daily[df_daily["area_id"] == area_id].drop(columns=["area_id"])
                df_daily_area.set_index("time", inplace=True)
                df_monthly_area = df_daily_area.resample('ME').mean().to_period('M')
                for var in rename_vars:
                    if var in df_monthly_area.columns:
                        path = f"/{A}/{area_id}/{var}//{E}/{F}/"
                        dcddss.write_rts(path, df_monthly_area[[var]], units, ptype)
    ds_daily.close()

def main():
    parser = create_argparser()
    args = parser.parse_args()

    with open(args.config_file, 'r') as f:
        config = yaml.safe_load(f)

    daily_to_monthly(config['daily_netcdf'], config['monthly_dss'], config)

if __name__ == "__main__":
    main()