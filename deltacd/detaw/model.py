"""High-level `xarray` orchestrator and CLI entry point for DETAW."""

from __future__ import annotations

import argparse
import logging
import timeit
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from deltacd.detaw.config import DetawConfig
from deltacd.detaw.et import ETcProvider, create_etc_provider
from deltacd.detaw.etaw import compute_smoothed_etaw
from deltacd.detaw.io import (
    read_crop_behavior_csv,
    read_crop_params,
    read_landuse,
    read_subarea_params,
    write_detaw_netcdf,
    write_weather_netcdf,
)
from deltacd.detaw.meteorology import MeteorologyProvider, create_meteorology_provider
from deltacd.detaw.params import ACRE_FT_PER_HA_MM_LEGACY
from deltacd.detaw.water_balance import run_soil_water_balance

logger = logging.getLogger(__name__)


class DetawModel:
    """Modular `xarray`-based DETAW simulation engine.

    Supports both full batch execution from YAML (`DetawModel.from_yaml(...).run()`)
    and interactive notebook / research workflows with custom meteorology or
    external ET (e.g., OpenET) providers.
    """

    def __init__(
        self,
        config: DetawConfig,
        meteorology_provider: MeteorologyProvider | None = None,
        etc_provider: ETcProvider | None = None,
    ):
        self.config = config
        if config.crop_behavior_path:
            self.extra_behaviors = (
                *config.extra_crop_behaviors,
                *read_crop_behavior_csv(config.crop_behavior_path),
            )
        else:
            self.extra_behaviors = config.extra_crop_behaviors

        self.meteorology_provider = meteorology_provider or create_meteorology_provider(config)
        self.etc_provider = etc_provider or create_etc_provider(
            config, extra_behaviors=self.extra_behaviors
        )

    @classmethod
    def from_yaml(cls, yaml_path: str | Path) -> "DetawModel":
        """Instantiate ``DetawModel`` from a YAML configuration file."""
        return cls(DetawConfig.from_yaml(yaml_path))

    def load_inputs(self) -> tuple[xr.Dataset, xr.Dataset, xr.Dataset, xr.Dataset]:
        """Load and coordinate-align landuse, crop parameters, subarea params, and meteorology."""
        water_years = list(range(self.config.start_water_year, self.config.end_water_year + 1))
        ds_landuse = read_landuse(
            self.config.landuse_path,
            water_years=water_years,
            extra_behaviors=self.extra_behaviors,
            column_map=self.config.column_map,
        )
        ds_crop_params = read_crop_params(
            self.config.critical_crop_path,
            self.config.noncritical_crop_path,
            target_crops=ds_landuse["crop"].values.tolist(),
            extra_behaviors=self.extra_behaviors,
        )
        ds_subarea = read_subarea_params(
            self.config.et_correction_path,
            subareas=ds_landuse["subarea"].values.tolist(),
            station_aliases=self.config.meteorology.station_aliases,
            column_map=self.config.column_map,
        )
        ds_meteo = self.meteorology_provider.load(
            self.config.start_water_year,
            self.config.end_water_year,
            ds_subarea,
        )
        return ds_landuse, ds_crop_params, ds_subarea, ds_meteo

    def run(self, write_outputs: bool = True) -> xr.Dataset:
        """Run the full DETAW simulation and return an ``xarray.Dataset``.

        Parameters
        ----------
        write_outputs : bool, default True
            If True, write ``precip_output``, ``et_output``, and ``detaw_output``
            NetCDF files specified in ``self.config``.

        Returns
        -------
        xarray.Dataset
            Dataset on ``(time, subarea, crop)`` with ``et_c``, ``s_e``,
            ``precip``, ``et_aw``, ``d_sw``, ``e_r``, and ``net_app``.
        """
        t0 = timeit.default_timer()
        ds_landuse, ds_crop_params, ds_subarea, ds_meteo_full = self.load_inputs()

        sim_start = pd.Timestamp(f"{self.config.start_water_year - 1}-10-01")
        sim_end = pd.Timestamp(f"{self.config.end_water_year}-09-30")
        time_index = pd.date_range(sim_start, sim_end, freq="D", name="time")

        if write_outputs:
            # Slice out any prior Sep 30 spin-up day so weather NetCDFs start on Oct 1
            write_weather_netcdf(
                ds_meteo_full.sel(time=time_index),
                self.config.precip_output_path,
                self.config.et_output_path,
            )

        ds_et = self.etc_provider.compute(
            time_index=time_index,
            ds_meteo_full=ds_meteo_full,
            ds_crop_params=ds_crop_params,
            da_year_type=ds_landuse["year_type"],
        )

        ds_wb = run_soil_water_balance(
            ds_et=ds_et,
            ds_crop_params=ds_crop_params,
            ds_subarea=ds_subarea,
            da_year_type=ds_landuse["year_type"],
            pad_year_type=self.config.options.default_pad_year_type,
            season_mode=self.config.et.season_mode,
            offseason_swd_cap=self.config.options.offseason_swd_cap,
            emulate_legacy_bugs=self.config.options.emulate_legacy_bugs,
        )

        if "external_et_aw_mm" in ds_et:
            da_etaw_mm = ds_et["external_et_aw_mm"]
        else:
            da_etaw_mm = compute_smoothed_etaw(
                ds_wb=ds_wb,
                ds_crop_params=ds_crop_params,
                cutoff_inclusive=self.config.options.cutoff_inclusive,
                emulate_legacy_bugs=self.config.options.emulate_legacy_bugs,
            )

        # Scale from depth (mm/day) to volume (acre-foot/day) when daily_output_unit == 1
        if self.config.daily_output_unit == 1:
            day_wy = xr.DataArray(
                np.where(time_index.month >= 10, time_index.year + 1, time_index.year),
                dims=("time",),
                coords={"time": time_index},
            )
            area_ha_daily = (
                ds_landuse["area_ha"]
                .sel(water_year=day_wy)
                .drop_vars("water_year")
                .transpose("time", "subarea", "crop")
            )
            scale = area_ha_daily * ACRE_FT_PER_HA_MM_LEGACY
            units_str = "acre-foot day-1"
        else:
            scale = 1.0
            units_str = "mm day-1"

        pcp_3d_mm = ds_wb["precip_mm"].broadcast_like(ds_wb["et_c_mm"])

        ds_result = xr.Dataset(
            data_vars={
                "et_c": (ds_wb["et_c_mm"] * scale).assign_attrs(
                    long_name="crop evapotranspiration", units=units_str
                ),
                "s_e": (ds_wb["s_e_mm"] * scale).assign_attrs(
                    long_name="Effective seepage", units=units_str
                ),
                "precip": (pcp_3d_mm * scale).assign_attrs(
                    long_name="precipitation", units=units_str
                ),
                "et_aw": (da_etaw_mm * scale).assign_attrs(
                    long_name="evapotranspiration of applied water", units=units_str
                ),
                "d_sw": (ds_wb["d_sw_mm"] * scale).assign_attrs(
                    long_name="change in soil water content", units=units_str
                ),
                "e_r": (ds_wb["e_r_mm"] * scale).assign_attrs(
                    long_name="effective rainfall", units=units_str
                ),
                "net_app": (ds_wb["net_app_mm"] * scale).assign_attrs(
                    long_name="net irrigation application", units=units_str
                ),
            },
            coords={
                "time": time_index,
                "subarea": ds_landuse["subarea"],
                "crop": ds_landuse["crop"],
                "crop_code": ds_landuse["crop_code"],
            },
        )

        if not self.config.daily_output:
            for v in ds_result.data_vars:
                ds_result[v] = xr.zeros_like(ds_result[v])

        if write_outputs:
            write_detaw_netcdf(
                ds_result,
                self.config.detaw_output_path,
                include_net_app=self.config.options.include_net_app_output,
            )

        elapsed = timeit.default_timer() - t0
        logger.info("DETAW simulation completed in %.3f seconds", elapsed)
        return ds_result


def run_detaw(yaml_path: str | Path, write_outputs: bool = True) -> xr.Dataset:
    """Load YAML configuration, run DETAW, and return the output ``xr.Dataset``."""
    model = DetawModel.from_yaml(yaml_path)
    return model.run(write_outputs=write_outputs)


def detaw(fname_main_yaml: str) -> None:
    """Backward-compatible entry point running DETAW from a YAML configuration path."""
    logging.info("Reading the main YAML input: %s", fname_main_yaml)
    run_detaw(fname_main_yaml, write_outputs=True)


def create_argparser() -> argparse.ArgumentParser:
    """Create CLI argument parser for ``detaw``."""
    parser = argparse.ArgumentParser(
        description="Run Delta Evapotranspiration of Applied Water (DETAW)"
    )
    parser.add_argument(
        "input_yaml",
        type=str,
        help="Path to the DETAW YAML configuration file",
    )
    return parser


def main() -> None:
    """CLI entry point for ``detaw``."""
    parser = create_argparser()
    args = parser.parse_args()
    detaw(args.input_yaml)
