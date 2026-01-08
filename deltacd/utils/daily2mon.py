import xarray as xr
import argparse
import pandas as pd
import pyhecdss

def create_argparser():
    """ Create an argument parser.

        Returns
        -------
        argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("daily_netcdf", type=str,
                        help="A str path that refers to the daily netCDF file")
    parser.add_argument("monthly_netcdf", type=str,
                        help="A str path that refers to the monthly netCDF file")

    return parser

def daily_to_monthly(daily_netcdf, monthly_dss):
    """ Convert daily netCDF data to monthly netCDF data.

        Parameters
        ----------
        daily_netcdf : str
            Path to the daily netCDF file.
        monthly_netcdf : str
            Path to the output monthly netCDF file.
    """
    # Open the daily netCDF file
    ds_daily = xr.open_dataset(daily_netcdf)
    df_daily=ds_daily.to_dataframe().reset_index()
    df_daily = df_daily.rename({"DIVERSION": "DIV-FLOW",
                              "DRAINAGE":"DRAIN-FLOW",
                              "SEEPAGE":"SEEP-FLOW"})
    with pyhecdss.DSSFile(str(monthly_dss), create_new=True) as dcddss:
        ptype = "PER-AVER"
        units = "CFS"
        A = "DELTACD"
        E = "1MON"
        F = "DWR-BDO"

        for node in df_daily["node"].unique():
            df_daily_node = df_daily[df_daily["node"] == node].drop(columns=["node"])
            df_daily_node.set_index("time", inplace=True)
            df_monthly_node = df_daily_node.resample('ME').mean().to_period('M')
            for var in df_monthly_node.columns:
                path = f"/{A}/{node}/{var}//{E}/{F}/"
                dcddss.write_rts(path, df_monthly_node, units, ptype)
    ds_daily.close()

def main():
    parser = create_argparser()
    args = parser.parse_args()

    daily_netcdf = args.daily_netcdf
    monthly_netcdf = args.monthly_netcdf

    daily_to_monthly(daily_netcdf, monthly_netcdf)

if __name__ == "__main__":
    main()