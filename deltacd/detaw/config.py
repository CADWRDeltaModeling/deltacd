"""Typed configuration schema and YAML loader for DETAW."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from deltacd.detaw.params import CropBehavior
from deltacd.utils.paths import convert_to_absolute_path


@dataclass
class MeteorologyConfig:
    """Configuration for meteorology & reference evapotranspiration (ETo)."""

    provider: str = "station"
    eto_method: str = "hargreaves_samani"
    input_pcp: str | None = None
    input_temperature: str | None = None
    dataset_path: str | None = None
    latitude_deg: float = 38.5
    elevation_m: float = 12.2
    temp_qc_mode: str = "clamp"  # "clamp" (legacy baseline) or "swap"
    station_aliases: dict[str, str] = field(default_factory=dict)


@dataclass
class ETConfig:
    """Configuration for crop evapotranspiration (ETc) / external ET."""

    provider: str = "kc_eto"  # "kc_eto" (default) or "openet" / "external"
    okc_mode: str = "legacy_cumulative"  # "legacy_cumulative" or "daily_increment"
    beta: float = 2.6
    wet_threshold_mode: str = "mean_eto"  # "mean_eto" (legacy) or "two_eto"
    season_mode: str = "legacy_wrapped"  # "legacy_wrapped" or "continuous"
    type1_initial_kc_mode: str = "legacy"  # "legacy" or "irrigation_frequency"
    external_path: str | None = None
    external_variable: str = "et_c"
    external_units: str = "mm/d"
    temporal_resolution: str = "daily"  # "daily" or "monthly"
    disaggregation: str = "et0_scaled"  # "uniform" or "et0_scaled"


@dataclass
class PhysicsOptions:
    """Configurable hydrological options and legacy-emulation flags."""

    offseason_swd_cap: bool = False
    cutoff_inclusive: bool = True
    include_net_app_output: bool = False
    default_pad_year_type: str = "AN"
    # Diagnostic toggles used strictly in equivalence tests to isolate Category A fixes:
    emulate_legacy_bugs: bool = False


@dataclass
class DetawConfig:
    """Complete resolved configuration for a DETAW simulation."""

    start_water_year: int
    end_water_year: int
    landuse_path: str
    et_correction_path: str
    critical_crop_path: str
    noncritical_crop_path: str
    detaw_output_path: str
    precip_output_path: str | None = None
    et_output_path: str | None = None
    crop_behavior_path: str | None = None
    daily_output: bool = True
    daily_output_unit: int = 1  # 1 = acre-foot day-1, 0 = mm day-1
    target_model: str | None = None
    meteorology: MeteorologyConfig = field(default_factory=MeteorologyConfig)
    et: ETConfig = field(default_factory=ETConfig)
    options: PhysicsOptions = field(default_factory=PhysicsOptions)
    column_map: dict[str, str] = field(default_factory=dict)
    extra_crop_behaviors: tuple[CropBehavior, ...] = field(default_factory=tuple)
    base_dir: Path = field(default_factory=Path.cwd)

    @property
    def model_start_year(self) -> int:
        """Calendar year containing October 1 of ``start_water_year``."""
        return int(self.start_water_year) - 1

    @classmethod
    def from_yaml(cls, yaml_path: str | Path) -> "DetawConfig":
        """Load and validate a ``DetawConfig`` from a YAML file."""
        path = Path(yaml_path).resolve()
        with open(path, "r") as f:
            raw = yaml.safe_load(f)
        if not isinstance(raw, dict) or "detaw" not in raw:
            raise ValueError(f"Expected top-level 'detaw:' mapping in {path}")
        return cls.from_dict(raw["detaw"], base_dir=path.parent)

    @classmethod
    def from_dict(cls, params: dict[str, Any], base_dir: str | Path | None = None) -> "DetawConfig":
        """Construct ``DetawConfig`` from a parsed ``detaw`` dictionary."""
        base = Path(base_dir).resolve() if base_dir is not None else Path.cwd()

        def _abs(val: str | None) -> str | None:
            if val is None:
                return None
            return convert_to_absolute_path(val, base)

        start_wy = int(params["start_water_year"])
        end_wy = int(params["end_water_year"])
        if end_wy < start_wy:
            raise ValueError(f"end_water_year ({end_wy}) must be >= start_water_year ({start_wy})")

        met_block = dict(params.get("meteorology") or {})
        met_cfg = MeteorologyConfig(
            provider=str(met_block.get("provider", "station")).lower(),
            eto_method=str(met_block.get("eto_method", "hargreaves_samani")).lower(),
            input_pcp=_abs(met_block.get("input_pcp", params.get("input_pcp"))),
            input_temperature=_abs(
                met_block.get("input_temperature", params.get("input_temperature"))
            ),
            dataset_path=_abs(met_block.get("dataset_path")),
            latitude_deg=float(met_block.get("latitude_deg", 38.5)),
            elevation_m=float(met_block.get("elevation_m", 12.2)),
            temp_qc_mode=str(met_block.get("temp_qc_mode", "clamp")).lower(),
            station_aliases=dict(met_block.get("station_aliases") or {}),
        )

        et_block = dict(params.get("et") or {})
        et_cfg = ETConfig(
            provider=str(et_block.get("provider", "kc_eto")).lower(),
            okc_mode=str(et_block.get("okc_mode", "legacy_cumulative")).lower(),
            beta=float(et_block.get("beta", 2.6)),
            wet_threshold_mode=str(et_block.get("wet_threshold_mode", "mean_eto")).lower(),
            season_mode=str(et_block.get("season_mode", "legacy_wrapped")).lower(),
            type1_initial_kc_mode=str(et_block.get("type1_initial_kc_mode", "legacy")).lower(),
            external_path=_abs(et_block.get("path") or et_block.get("external_path")),
            external_variable=str(
                et_block.get("variable", et_block.get("external_variable", "et_c"))
            ),
            external_units=str(et_block.get("units", et_block.get("external_units", "mm/d"))),
            temporal_resolution=str(et_block.get("temporal_resolution", "daily")).lower(),
            disaggregation=str(et_block.get("disaggregation", "et0_scaled")).lower(),
        )

        opt_block = dict(params.get("options") or {})
        opt_cfg = PhysicsOptions(
            offseason_swd_cap=bool(opt_block.get("offseason_swd_cap", False)),
            cutoff_inclusive=bool(opt_block.get("cutoff_inclusive", True)),
            include_net_app_output=bool(opt_block.get("include_net_app_output", False)),
            default_pad_year_type=str(opt_block.get("default_pad_year_type", "AN")),
            emulate_legacy_bugs=bool(opt_block.get("emulate_legacy_bugs", False)),
        )

        for req_key in (
            "landuse",
            "et_correction",
            "critical",
            "noncritical",
            "detaw_output",
        ):
            if req_key not in params or params[req_key] is None:
                raise ValueError(f"Missing required DETAW config key: {req_key!r}")

        return cls(
            start_water_year=start_wy,
            end_water_year=end_wy,
            landuse_path=_abs(params["landuse"]),  # type: ignore[arg-type]
            et_correction_path=_abs(params["et_correction"]),  # type: ignore[arg-type]
            critical_crop_path=_abs(params["critical"]),  # type: ignore[arg-type]
            noncritical_crop_path=_abs(params["noncritical"]),  # type: ignore[arg-type]
            detaw_output_path=_abs(params["detaw_output"]),  # type: ignore[arg-type]
            precip_output_path=_abs(params.get("precip_output")),
            et_output_path=_abs(params.get("et_output")),
            crop_behavior_path=_abs(params.get("crop_behavior")),
            daily_output=bool(int(params.get("daily_output", 1))),
            daily_output_unit=int(params.get("daily_output_unit", 1)),
            target_model=params.get("target_model"),
            meteorology=met_cfg,
            et=et_cfg,
            options=opt_cfg,
            column_map=dict(params.get("column_map") or {}),
            base_dir=base,
        )
