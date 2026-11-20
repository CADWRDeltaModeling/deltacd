"""Delta Evapotranspiration of Applied Water (DETAW) subpackage.

Provides an extensible, coordinate-aligned `xarray` & `pandas` implementation of
DETAW (**Snyder et al., 2006**) while preserving backward compatibility with
``deltacd.detaw.detaw(yaml_path)`` and the ``detaw`` CLI command.
"""

from __future__ import annotations

from deltacd.detaw.config import DetawConfig, ETConfig, MeteorologyConfig, PhysicsOptions
from deltacd.detaw.et import (
    ETcProvider,
    ExternalETProvider,
    KcETcProvider,
    build_daily_ikc_and_kc,
    compute_bare_soil_okc,
)
from deltacd.detaw.etaw import compute_smoothed_etaw
from deltacd.detaw.io import (
    read_crop_behavior_csv,
    read_crop_params,
    read_landuse,
    read_station_precipitation,
    read_subarea_params,
    read_weather_station_csv,
    write_detaw_netcdf,
    write_weather_netcdf,
)
from deltacd.detaw.meteorology import (
    DatasetMeteorologyProvider,
    MeteorologyProvider,
    StationMeteorologyProvider,
    asce_penman_monteith,
    extraterrestrial_radiation,
    hargreaves_samani,
)
from deltacd.detaw.model import DetawModel, create_argparser, detaw, main, run_detaw
from deltacd.detaw.params import (
    ACRE_FT_PER_HA_MM_EXACT,
    ACRE_FT_PER_HA_MM_LEGACY,
    DEFAULT_CROP_BEHAVIORS,
    CropBehavior,
    resolve_crop_behavior,
)
from deltacd.detaw.water_balance import run_soil_water_balance
from deltacd.utils.paths import convert_to_absolute_path

__all__ = [
    "ACRE_FT_PER_HA_MM_EXACT",
    "ACRE_FT_PER_HA_MM_LEGACY",
    "CropBehavior",
    "DEFAULT_CROP_BEHAVIORS",
    "DatasetMeteorologyProvider",
    "DetawConfig",
    "DetawModel",
    "ETConfig",
    "ETcProvider",
    "ExternalETProvider",
    "KcETcProvider",
    "MeteorologyConfig",
    "MeteorologyProvider",
    "PhysicsOptions",
    "StationMeteorologyProvider",
    "asce_penman_monteith",
    "build_daily_ikc_and_kc",
    "compute_bare_soil_okc",
    "compute_smoothed_etaw",
    "convert_to_absolute_path",
    "create_argparser",
    "detaw",
    "extraterrestrial_radiation",
    "hargreaves_samani",
    "main",
    "read_crop_behavior_csv",
    "read_crop_params",
    "read_landuse",
    "read_station_precipitation",
    "read_subarea_params",
    "read_weather_station_csv",
    "resolve_crop_behavior",
    "run_detaw",
    "run_soil_water_balance",
    "write_detaw_netcdf",
    "write_weather_netcdf",
]
