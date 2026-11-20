"""Declarative crop behavior metadata, alias registries, and physical constants.

References
----------
Snyder, R. L., Orang, M. N., Geng, S., & Matyac, J. S. (2006). *DETAW: Delta
Evapotranspiration of Applied Water*, Version 1.0 Final Report. California
Department of Water Resources and University of California, Davis.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

# Unit conversion constants matching legacy DETAW / Snyder et al. (2006, p. 28).
# Note: Legacy detaw.py converts hectares -> acres via 2.471 and mm -> ft via
# 0.00328084, giving 2.471 * 0.00328084 = 0.00810695564 acre-ft / (ha * mm)
# in the daily output array, whereas Snyder et al. (2006, p. 28) states
# 10,000 m2/ha * 0.001 m/mm / 1233.482 m3/af = 0.00810713194.
HECTARES_TO_ACRES: float = 2.471
MM_TO_FEET: float = 0.00328084
ACRE_FT_PER_HA_MM_LEGACY: float = HECTARES_TO_ACRES * MM_TO_FEET
ACRE_FT_PER_HA_MM_EXACT: float = 10.0 / 1233.48183754752


def normalize_token(name: str) -> str:
    """Normalize a column, station, or crop token for case/space-insensitive matching."""
    return re.sub(r"[^a-z0-9]+", "", str(name).strip().lower())


@dataclass(frozen=True)
class CropBehavior:
    """Declarative hydrological behavior for a land-use / crop category.

    Replaces hardcoded ``crop_i == 6, 7, 12, 13, 14, 15, 16`` branches in the
    legacy C++-transcribed DETAW loop so new crops can be added purely via data.

    Attributes
    ----------
    code : str
        Canonical short code (e.g., ``"TO"``, ``"GR"``, ``"WS"``, ``"DP"``).
    label : str
        Canonical NetCDF output label expected by downstream DCD (e.g.,
        ``"Tomato"``, ``"Water surface"``, ``"Duck pond"``).
    aliases : tuple of str
        Alternative codes or labels encountered in historical CSVs.
    seepage_rate_per_month : float
        Monthly potential seepage per unit effective root depth ``erd`` in
        lowland subareas (``Region == 0``). Standard agricultural crops use
        ``0.025`` (0.3 in / ft / month; Snyder et al., 2006, p. 37), native &
        riparian vegetation use ``0.15``, duck pond uses ``0.06``, and open
        water uses ``0.0``.
    is_irrigated : bool
        Whether net irrigation applications (``NetApp``) are applied when
        ``SWD + Dsw > SWDx`` during the growing season.
    is_open_water : bool
        Open-water surface mode: ``Kc`` fixed to ``open_water_kc`` (``1.1``),
        ``Espg = 0``, ``ERn = PCP``, ``Dsw = ETc - PCP``, ``NetApp = 0``, and
        excluded from agricultural applied water in DCD.
    open_water_kc : float
        Fixed ``Kc`` used when ``is_open_water`` is ``True`` (``1.1``; Snyder
        et al., 2006, Table 7, p. 19).
    riparian_full_seepage : bool
        Riparian vegetation mode: in lowlands (``Region == 0``), ``Espg = ETc``,
        ``ERn = 0``, ``Dsw = 0``, ``NetApp = 0`` (Snyder et al., 2006, p. 31).
    continuous_flooding : bool
        Rice paddy mode: during the in-season period, ``Espg = 0``, ``ERn = 0``,
        and daily irrigation ``NetApp = SWD + Dsw`` up to
        ``pre_harvest_drain_days`` before harvest (Snyder et al., 2006, p. 30).
    pre_harvest_drain_days : int
        Days prior to ``end_doy`` when continuous flooding stops (20 days for
        rice).
    use_offseason_swdx_in_season : bool
        Use ``osSWDx`` (top 0.3 m allowable depletion) as the in-season
        scheduling threshold ``SWDx`` (used for Duck Pond / managed wetlands).
    year_type_unshifted_in_legacy : bool
        In legacy ``detaw.py``, fall-planted grain (``GR``, ``DG``/``DGR``)
        takes the water-year type of ``y`` rather than ``y - 1`` when indexing
        calendar years.
    grain_cutoff_rule : bool
        Apply the split fall/spring ``DOYGrainLIrrig`` cutoff rule in legacy
        smoothed ``et_aw`` post-processing (``GR``).
    """

    code: str
    label: str
    aliases: tuple[str, ...] = field(default_factory=tuple)
    seepage_rate_per_month: float = 0.025
    is_irrigated: bool = True
    is_open_water: bool = False
    open_water_kc: float = 1.1
    riparian_full_seepage: bool = False
    continuous_flooding: bool = False
    pre_harvest_drain_days: int = 0
    use_offseason_swdx_in_season: bool = False
    year_type_unshifted_in_legacy: bool = False
    grain_cutoff_rule: bool = False


DEFAULT_CROP_BEHAVIORS: tuple[CropBehavior, ...] = (
    CropBehavior(code="UR", label="Urban", aliases=("urban", "u")),
    CropBehavior(
        code="PA",
        label="Irrig pasture",
        aliases=("pasture", "irrigpasture", "irrigatedpasture"),
    ),
    CropBehavior(code="AL", label="Alfalfa", aliases=("alfalfa",)),
    CropBehavior(
        code="FI",
        label="All field",
        aliases=("fl", "f", "field", "allfield", "fieldcrops", "irrigfield"),
    ),
    CropBehavior(
        code="SB",
        label="Sugar beets",
        aliases=("sugarbeets", "sugarbeet"),
    ),
    CropBehavior(
        code="GR",
        label="Irrig grain",
        aliases=("grain", "grains", "irriggrain", "irriggrains"),
        year_type_unshifted_in_legacy=True,
        grain_cutoff_rule=True,
    ),
    CropBehavior(
        code="RI",
        label="Rice",
        aliases=("rice",),
        continuous_flooding=True,
        pre_harvest_drain_days=20,
        use_offseason_swdx_in_season=True,
    ),
    CropBehavior(
        code="TR",
        label="Truck crops",
        aliases=("truck", "truckcrops"),
    ),
    CropBehavior(
        code="TO",
        label="Tomato",
        aliases=("tomato", "tomatoes", "tomatos"),
    ),
    CropBehavior(
        code="OR",
        label="Orchard",
        aliases=("orchards", "orchard"),
    ),
    CropBehavior(
        code="VI",
        label="Vineyard",
        aliases=("vineyards", "vineyard"),
    ),
    CropBehavior(
        code="RV",
        label="Riparian vegetation",
        aliases=("riparian", "riparianvegetation", "nativeriparian"),
        seepage_rate_per_month=0.15,
        is_irrigated=True,  # Upland riparian (where channel Espg == 0) triggers irrigation when SWD+Dsw > osSWDx; lowland has Dsw == 0 so NetApp == 0
        riparian_full_seepage=True,
        use_offseason_swdx_in_season=True,
    ),
    CropBehavior(
        code="NV",
        label="Native vegetation",
        aliases=("native", "nativevegetation"),
        seepage_rate_per_month=0.15,
        is_irrigated=False,
    ),
    CropBehavior(
        code="DGR",
        label="Non-irrig grain",
        aliases=("dg", "drygrain", "nonirriggrain", "nonirrigatedgrain"),
        is_irrigated=False,
        year_type_unshifted_in_legacy=True,
    ),
    CropBehavior(
        code="WS",
        label="Water surface",
        aliases=("water", "watersurface", "watersurfaces", "openwater"),
        seepage_rate_per_month=0.0,
        is_irrigated=False,
        is_open_water=True,
        open_water_kc=1.1,
        use_offseason_swdx_in_season=True,
    ),
    CropBehavior(
        code="DP",
        label="Duck pond",
        aliases=("rvw", "rv_w", "duckpond", "duck"),
        seepage_rate_per_month=0.06,
        is_irrigated=True,
        use_offseason_swdx_in_season=True,
    ),
)


def build_crop_registry(
    extra_behaviors: list[CropBehavior] | tuple[CropBehavior, ...] = (),
) -> tuple[dict[str, CropBehavior], dict[str, CropBehavior]]:
    """Build lookup maps by canonical code and by normalized token/alias."""
    by_code: dict[str, CropBehavior] = {}
    by_token: dict[str, CropBehavior] = {}
    for cb in (*DEFAULT_CROP_BEHAVIORS, *extra_behaviors):
        by_code[cb.code] = cb
        by_token[normalize_token(cb.code)] = cb
        by_token[normalize_token(cb.label)] = cb
        for alias in cb.aliases:
            by_token[normalize_token(alias)] = cb
    return by_code, by_token


def resolve_crop_behavior(
    token: str,
    extra_behaviors: list[CropBehavior] | tuple[CropBehavior, ...] = (),
) -> CropBehavior:
    """Resolve a crop code or display name to its ``CropBehavior``."""
    _, by_token = build_crop_registry(extra_behaviors)
    key = normalize_token(token)
    if key not in by_token:
        known = sorted({cb.code for cb in (*DEFAULT_CROP_BEHAVIORS, *extra_behaviors)})
        raise ValueError(
            f"Unknown crop code or label {token!r}. Known canonical crop codes: {known}"
        )
    return by_token[key]


# Canonical station alias registry so station columns in precipitation CSVs and
# station weight columns in et_correction / Percentage CSVs match by name
# regardless of column ordering or historical naming variants.
STATION_ALIASES: dict[str, str] = {
    "brentwood": "brentwood",
    "davis": "davis",
    "galt": "galt",
    "lodi": "lodi",
    "stockton": "stockton",
    "riovista": "riovista",
    "tracy": "tracy_carbona",
    "tracycarbona": "tracy_carbona",
    "concord": "concord_riovista",
    "concordriovista": "concord_riovista",
}


def canonical_station_id(name: str, custom_aliases: dict[str, str] | None = None) -> str:
    """Normalize a meteorological station column header to a canonical station ID."""
    norm = normalize_token(name)
    if custom_aliases:
        normalized_custom = {
            normalize_token(k): normalize_token(v) for k, v in custom_aliases.items()
        }
        if norm in normalized_custom:
            norm = normalized_custom[norm]
    return STATION_ALIASES.get(norm, norm)
