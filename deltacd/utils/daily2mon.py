import xarray as xr
import argparse

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

def daily_to_monthly(daily_netcdf, monthly_netcdf):
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

    # Resample the data to monthly frequency
    ds_monthly = ds_daily.resample(time="ME").mean()

    # Save the monthly netCDF file
    ds_monthly.to_netcdf(monthly_netcdf)

def main():
    parser = create_argparser()
    args = parser.parse_args()

    daily_netcdf = args.daily_netcdf
    monthly_netcdf = args.monthly_netcdf

    daily_to_monthly(daily_netcdf, monthly_netcdf)

if __name__ == "__main__":
    main()