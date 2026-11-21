"""Order-independent pandas/xarray readers and CF NetCDF writers for DETAW.

Every reader in this module resolves columns, rows, stations, subareas, and
crops by normalized header/label names and alias registries—never by fixed
column positions, row numbers, or hardcoded 15/16 crop counts.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import xarray as xr

import deltacd
from deltacd.detaw.params import (
    CropBehavior,
    canonical_station_id,
    normalize_token,
    resolve_crop_behavior,
)


def _find_column(
    columns: Iterable[str],
    candidates: tuple[str, ...],
    column_map: dict[str, str] | None = None,
    required: bool = True,
    context: str = "",
) -> str | None:
    """Find a DataFrame column matching any candidate token (order-independent)."""
    col_list = list(columns)
    if column_map:
        for canonical, user_col in column_map.items():
            if normalize_token(canonical) in {normalize_token(c) for c in candidates}:
                for col in col_list:
                    if col == user_col or normalize_token(col) == normalize_token(user_col):
                        return col

    norm_to_col: dict[str, str] = {}
    for col in col_list:
        norm = normalize_token(col)
        if norm not in norm_to_col:
            norm_to_col[norm] = col

    for cand in candidates:
        norm_cand = normalize_token(cand)
        if norm_cand in norm_to_col:
            return norm_to_col[norm_cand]

    if required:
        raise ValueError(
            f"Could not find required column matching {candidates} in {context or 'CSV'}. "
            f"Available columns: {col_list}"
        )
    return None


def _coerce_subarea_id(series: pd.Series) -> pd.Series:
    """Normalize subarea IDs to integers when all values are numeric, else stripped strings."""
    cleaned = series.astype(str).str.strip()
    try:
        return pd.to_numeric(cleaned, errors="raise").astype(np.int32)
    except (ValueError, TypeError):
        return cleaned


def read_crop_behavior_csv(filepath: str | Path) -> tuple[CropBehavior, ...]:
    """Read an optional custom crop-behavior CSV (order-independent columns)."""
    df = pd.read_csv(filepath)
    code_col = _find_column(df.columns, ("code", "crop_code", "crop"), context=str(filepath))
    label_col = _find_column(
        df.columns, ("label", "name", "crop_name"), required=False, context=str(filepath)
    )
    aliases_col = _find_column(
        df.columns, ("aliases", "alias"), required=False, context=str(filepath)
    )
    seep_col = _find_column(
        df.columns,
        ("seepage_rate_per_month", "seepage_fraction", "seepage_rate"),
        required=False,
        context=str(filepath),
    )
    irrig_col = _find_column(
        df.columns, ("is_irrigated", "irrigated"), required=False, context=str(filepath)
    )
    ow_col = _find_column(
        df.columns,
        ("is_open_water", "water_surface", "open_water"),
        required=False,
        context=str(filepath),
    )
    ow_kc_col = _find_column(df.columns, ("open_water_kc",), required=False, context=str(filepath))
    rv_col = _find_column(
        df.columns,
        ("riparian_full_seepage", "riparian"),
        required=False,
        context=str(filepath),
    )
    flood_col = _find_column(
        df.columns,
        ("continuous_flooding", "flooded", "rice_mode"),
        required=False,
        context=str(filepath),
    )
    drain_col = _find_column(
        df.columns,
        ("pre_harvest_drain_days", "drain_days"),
        required=False,
        context=str(filepath),
    )
    os_swdx_col = _find_column(
        df.columns,
        ("use_offseason_swdx_in_season", "use_offseason_swdx"),
        required=False,
        context=str(filepath),
    )

    behaviors: list[CropBehavior] = []
    for _, row in df.iterrows():
        code = str(row[code_col]).strip()  # type: ignore[index]
        label = str(row[label_col]).strip() if label_col and pd.notna(row[label_col]) else code
        aliases: tuple[str, ...] = ()
        if aliases_col and pd.notna(row[aliases_col]):
            aliases = tuple(
                a.strip() for a in str(row[aliases_col]).replace(";", "|").split("|") if a.strip()
            )
        behaviors.append(
            CropBehavior(
                code=code,
                label=label,
                aliases=aliases,
                seepage_rate_per_month=(
                    float(row[seep_col]) if seep_col and pd.notna(row[seep_col]) else 0.025
                ),
                is_irrigated=(
                    bool(row[irrig_col]) if irrig_col and pd.notna(row[irrig_col]) else True
                ),
                is_open_water=bool(row[ow_col]) if ow_col and pd.notna(row[ow_col]) else False,
                open_water_kc=(
                    float(row[ow_kc_col]) if ow_kc_col and pd.notna(row[ow_kc_col]) else 1.1
                ),
                riparian_full_seepage=(
                    bool(row[rv_col]) if rv_col and pd.notna(row[rv_col]) else False
                ),
                continuous_flooding=(
                    bool(row[flood_col]) if flood_col and pd.notna(row[flood_col]) else False
                ),
                pre_harvest_drain_days=(
                    int(row[drain_col]) if drain_col and pd.notna(row[drain_col]) else 0
                ),
                use_offseason_swdx_in_season=(
                    bool(row[os_swdx_col]) if os_swdx_col and pd.notna(row[os_swdx_col]) else False
                ),
            )
        )
    return tuple(behaviors)


def read_landuse(
    filepath: str | Path,
    water_years: Iterable[int] | None = None,
    extra_behaviors: tuple[CropBehavior, ...] = (),
    column_map: dict[str, str] | None = None,
) -> xr.Dataset:
    """Read subarea land-use acreages (in hectares) and water-year types.

    Resolves ``area_id``, ``DATE`` (water year), ``TYPE`` (water-year type), and
    all crop columns by header name regardless of column order or row ordering.
    Fixes Category A5 bug where legacy ``detaw.py`` used ``iloc[:, 1]`` and read
    ``DATE`` instead of ``TYPE`` in ``landuse_dsm2.csv``.

    Returns
    -------
    xarray.Dataset
        Dataset with coordinates ``(water_year, subarea, crop)`` and ``crop_code``:
        - ``area_ha(water_year, subarea, crop)`` : float64 (hectares)
        - ``year_type(water_year)`` : str (``"W"``, ``"AN"``, ``"BN"``, ``"D"``, ``"C"``)
    """
    df = pd.read_csv(filepath)
    area_col = _find_column(
        df.columns,
        ("area_id", "subarea", "subarea_id", "island"),
        column_map=column_map,
        context=f"landuse CSV ({filepath})",
    )
    wy_col = _find_column(
        df.columns,
        ("date", "water_year", "wy", "year"),
        column_map=column_map,
        context=f"landuse CSV ({filepath})",
    )
    type_col = _find_column(
        df.columns,
        ("type", "year_type", "wy_type", "wyt"),
        column_map=column_map,
        context=f"landuse CSV ({filepath})",
    )
    assert area_col is not None and wy_col is not None and type_col is not None

    meta_cols = {area_col, wy_col, type_col}
    # Ignore optional non-crop metadata columns like 'extension' if present
    ext_col = _find_column(df.columns, ("extension",), required=False)
    if ext_col is not None:
        meta_cols.add(ext_col)

    raw_crop_cols = [c for c in df.columns if c not in meta_cols]
    if not raw_crop_cols:
        raise ValueError(f"No crop columns found in landuse CSV: {filepath}")

    crop_behaviors: list[CropBehavior] = []
    seen_codes: set[str] = set()
    for col in raw_crop_cols:
        cb = resolve_crop_behavior(col, extra_behaviors=extra_behaviors)
        if cb.code in seen_codes:
            raise ValueError(
                f"Duplicate crop columns resolving to {cb.code!r} in {filepath}: {raw_crop_cols}"
            )
        seen_codes.add(cb.code)
        crop_behaviors.append(cb)

    work = df.copy()
    work[area_col] = _coerce_subarea_id(work[area_col])
    work[wy_col] = pd.to_numeric(work[wy_col], errors="raise").astype(int)
    work[type_col] = work[type_col].astype(str).str.strip().str.upper()

    if water_years is not None:
        requested_wy = [int(y) for y in water_years]
        available_wy = set(work[wy_col].unique())
        missing_wy = sorted(set(requested_wy) - available_wy)
        if missing_wy:
            raise ValueError(
                f"Landuse CSV {filepath} is missing requested water years: {missing_wy}"
            )
        work = work[work[wy_col].isin(requested_wy)].copy()
    else:
        requested_wy = sorted(work[wy_col].unique().tolist())

    if work.duplicated(subset=[area_col, wy_col]).any():
        dups = work.loc[work.duplicated(subset=[area_col, wy_col]), [area_col, wy_col]]
        raise ValueError(f"Duplicate (subarea, water_year) rows in {filepath}:\n{dups.head()}")

    subareas = sorted(work[area_col].unique().tolist(), key=lambda x: (isinstance(x, str), x))

    # Validate year_type consistency across subareas per water year
    wy_types = work.groupby(wy_col)[type_col].nunique()
    if (wy_types > 1).any():
        bad_years = wy_types[wy_types > 1].index.tolist()
        raise ValueError(
            f"Inconsistent water-year TYPE across subareas in {filepath} for years: {bad_years}"
        )
    year_type_series = work.groupby(wy_col)[type_col].first().reindex(requested_wy)

    # Pivot explicitly by (water_year, subarea) so row order never matters
    area_arr = np.empty((len(requested_wy), len(subareas), len(raw_crop_cols)), dtype=np.float64)
    for idx, col in enumerate(raw_crop_cols):
        pivoted = work.pivot(index=wy_col, columns=area_col, values=col).reindex(
            index=requested_wy, columns=subareas
        )
        if pivoted.isna().any().any():
            raise ValueError(
                f"Missing landuse entries for crop {col!r} in {filepath} after aligning "
                f"(water_year, subarea)."
            )
        area_arr[:, :, idx] = pivoted.to_numpy(dtype=np.float64)

    crop_labels = [cb.label for cb in crop_behaviors]
    crop_codes = [cb.code for cb in crop_behaviors]

    ds = xr.Dataset(
        data_vars={
            "area_ha": (
                ("water_year", "subarea", "crop"),
                area_arr,
                {"units": "hectare", "long_name": "planted land-use area"},
            ),
            "year_type": (
                ("water_year",),
                year_type_series.to_numpy(dtype="<U4"),
                {"long_name": "Sacramento Valley water year classification"},
            ),
        },
        coords={
            "water_year": np.asarray(requested_wy, dtype=np.int32),
            "subarea": np.asarray(subareas),
            "crop": np.asarray(crop_labels, dtype=object),
            "crop_code": ("crop", np.asarray(crop_codes, dtype=object)),
        },
    )
    return ds


# Canonical parameter names expected in crop parameter datasets
CROP_PARAM_SPECS: dict[str, tuple[str, ...]] = {
    "crop_type": ("type", "croptype", "crop_type"),
    "begin_doy": ("bd_doy", "begindoy", "begin_doy"),
    "end_doy": ("ed_doy", "enddoy", "end_doy"),
    "irrig_freq": ("f", "if", "irrig_freq", "irrigation_frequency"),
    "kc1": ("kc1",),
    "kc2": ("kc2",),
    "kc3": ("kc3",),
    "ab_pct": ("a-b", "ab_pct", "pct_ab"),
    "ac_pct": ("a-c", "ac_pct", "pct_ac"),
    "ad_pct": ("a-d", "ad_pct", "pct_ad"),
    "max_soil_depth_mm": ("sdx", "max_soil_depth", "max_soil_depth_mm", "soil_depth"),
    "root_depth_lowland_mm": ("rdx_lo", "rdxlo", "root_depth_lowland_mm", "root_depth_lo"),
    "root_depth_upland_mm": ("rdx_up", "rdxup", "root_depth_upland_mm", "root_depth_up"),
    "aw_lowland": ("aw_lo", "awlo", "aw_lowland"),
    "aw_upland": ("aw_up", "awup", "aw_upland"),
    "allowable_depletion_pct": ("ad %", "ad%", "allowable_depletion_pct", "ad_percent"),
}


def _parse_single_crop_param_file(
    filepath: str | Path,
    extra_behaviors: tuple[CropBehavior, ...] = (),
) -> pd.DataFrame:
    """Parse either a legacy transposed ``crop_info_*.csv`` or a tidy row-per-crop CSV.

    Returns a DataFrame indexed by canonical ``crop_code`` with standardized
    numeric parameter columns in ``CROP_PARAM_SPECS``.
    """
    raw_preview = pd.read_csv(filepath, header=None, dtype=str)

    # Check if this is a tidy row-per-crop table (first row has 'code'/'crop' + 'kc1')
    first_row_tokens = {normalize_token(v) for v in raw_preview.iloc[0].dropna().tolist()}
    if "kc1" in first_row_tokens and (
        "code" in first_row_tokens or "crop" in first_row_tokens or "cropcode" in first_row_tokens
    ):
        tidy_df = pd.read_csv(filepath)
        code_col = _find_column(
            tidy_df.columns, ("code", "crop_code", "crop"), context=str(filepath)
        )
        assert code_col is not None
        records: dict[str, dict[str, float]] = {}
        for _, row in tidy_df.iterrows():
            cb = resolve_crop_behavior(str(row[code_col]), extra_behaviors=extra_behaviors)
            rec: dict[str, float] = {}
            for param_name, aliases in CROP_PARAM_SPECS.items():
                col = _find_column(tidy_df.columns, aliases, context=str(filepath))
                assert col is not None
                rec[param_name] = float(row[col])
            records[cb.code] = rec
        return pd.DataFrame.from_dict(records, orient="index")

    # Otherwise parse the legacy transposed DETAW CSV by row labels (no magic row offsets!)
    # Locate the row containing crop short codes ("Code" in col 0, or row whose cells match known crops)
    code_row_idx: int | None = None
    for r_idx in range(len(raw_preview)):
        first_cell = normalize_token(raw_preview.iloc[r_idx, 0] or "")
        if first_cell == "code":
            code_row_idx = r_idx
            break

    if code_row_idx is None:
        raise ValueError(f"Could not locate 'Code' row in legacy crop info CSV: {filepath}")

    # Identify which column indices correspond to crop columns on the 'Code' row.
    # In legacy files, col 0 is 'Code', col 1 is 'Type' (a secondary label column),
    # and cols >= 2 (or any col with a valid crop token) hold crop codes.
    col_to_crop_code: dict[int, str] = {}
    for c_idx in range(1, raw_preview.shape[1]):
        cell = str(raw_preview.iloc[code_row_idx, c_idx] or "").strip()
        if not cell or normalize_token(cell) in {"type", "#"}:
            continue
        cb = resolve_crop_behavior(cell, extra_behaviors=extra_behaviors)
        if cb.code in col_to_crop_code.values():
            raise ValueError(f"Duplicate crop code {cb.code!r} in {filepath}")
        col_to_crop_code[c_idx] = cb.code

    if not col_to_crop_code:
        raise ValueError(f"No crop codes found on 'Code' row in {filepath}")

    crop_cols = list(col_to_crop_code.keys())

    # Scan all other rows and classify by parameter label in column 0.
    # Note: Legacy files have TWO 'BD' rows and TWO 'ED' rows:
    # one calendar string ('14-Apr') and one numeric DOY ('105').
    # Selecting the row whose crop values parse as numbers handles them cleanly in any order!
    param_rows: dict[str, dict[str, float]] = {}
    label_lookup: dict[str, str] = {}
    for canonical_param, aliases in CROP_PARAM_SPECS.items():
        if canonical_param not in ("begin_doy", "end_doy"):
            for alias in aliases:
                label_lookup[normalize_token(alias)] = canonical_param

    for r_idx in range(len(raw_preview)):
        if r_idx == code_row_idx:
            continue
        raw_label = str(raw_preview.iloc[r_idx, 0] or "").strip()
        if not raw_label:
            continue
        lower_compact = re.sub(r"\s+", "", raw_label.lower())
        norm_label = normalize_token(raw_label)

        # Extract values across crop columns and check if all are numeric
        vals: dict[str, float] = {}
        all_numeric = True
        for c_idx in crop_cols:
            val_str = str(raw_preview.iloc[r_idx, c_idx] or "").strip()
            try:
                vals[col_to_crop_code[c_idx]] = float(val_str)
            except ValueError:
                all_numeric = False
                break

        if not all_numeric:
            continue

        if lower_compact == "ad%" or norm_label in ("allowabledepletionpct", "adpercent"):
            param_rows["allowable_depletion_pct"] = vals
        elif lower_compact == "a-b" or norm_label in ("abpct", "pctab"):
            param_rows["ab_pct"] = vals
        elif lower_compact == "a-c" or norm_label in ("acpct", "pctac"):
            param_rows["ac_pct"] = vals
        elif lower_compact == "a-d" or norm_label in ("adpct", "pctad"):
            param_rows["ad_pct"] = vals
        elif norm_label in ("bd", "begindoy", "bddoy"):
            param_rows["begin_doy"] = vals
        elif norm_label in ("ed", "enddoy", "eddoy"):
            param_rows["end_doy"] = vals
        elif norm_label in label_lookup:
            target = label_lookup[norm_label]
            # Avoid matching header metadata row 'Sub-Area / TYPE' if any before Code
            if target == "crop_type" and any(v not in (1.0, 2.0, 3.0, 4.0) for v in vals.values()):
                continue
            param_rows[target] = vals

    missing_params = sorted(set(CROP_PARAM_SPECS.keys()) - set(param_rows.keys()))
    if missing_params:
        raise ValueError(
            f"Crop parameter CSV {filepath} is missing required numeric rows: {missing_params}"
        )

    df_out = pd.DataFrame(param_rows)
    df_out.index.name = "crop_code"
    return df_out


def read_crop_params(
    critical_filepath: str | Path,
    noncritical_filepath: str | Path,
    target_crops: Iterable[str] | None = None,
    extra_behaviors: tuple[CropBehavior, ...] = (),
) -> xr.Dataset:
    """Read critical and non-critical crop/soil parameters into a tidy ``xr.Dataset``.

    Aligns crops by canonical ``crop_code`` (matching ``target_crops`` in any
    order) and attaches declarative ``CropBehavior`` coordinates/variables.
    """
    df_crit = _parse_single_crop_param_file(critical_filepath, extra_behaviors=extra_behaviors)
    df_noncrit = _parse_single_crop_param_file(
        noncritical_filepath, extra_behaviors=extra_behaviors
    )

    if target_crops is not None:
        ordered_behaviors = [
            resolve_crop_behavior(c, extra_behaviors=extra_behaviors) for c in target_crops
        ]
    else:
        shared_codes = [c for c in df_noncrit.index if c in df_crit.index]
        ordered_behaviors = [
            resolve_crop_behavior(c, extra_behaviors=extra_behaviors) for c in shared_codes
        ]

    ordered_codes = [cb.code for cb in ordered_behaviors]
    ordered_labels = [cb.label for cb in ordered_behaviors]

    missing_crit = sorted(set(ordered_codes) - set(df_crit.index))
    missing_noncrit = sorted(set(ordered_codes) - set(df_noncrit.index))
    if missing_crit or missing_noncrit:
        raise ValueError(
            f"Missing crop parameter definitions for landuse crops. "
            f"Missing in critical ({critical_filepath}): {missing_crit}; "
            f"missing in noncritical ({noncritical_filepath}): {missing_noncrit}."
        )

    df_crit = df_crit.reindex(ordered_codes)
    df_noncrit = df_noncrit.reindex(ordered_codes)

    regimes = np.asarray(["noncritical", "critical"], dtype=object)
    data_vars: dict[str, tuple[tuple[str, ...], np.ndarray]] = {}

    int_params = {"crop_type", "begin_doy", "end_doy", "irrig_freq"}
    for param_name in CROP_PARAM_SPECS:
        stacked = np.stack(
            [
                df_noncrit[param_name].to_numpy(dtype=np.float64),
                df_crit[param_name].to_numpy(dtype=np.float64),
            ],
            axis=0,
        )
        if param_name in int_params:
            stacked = np.rint(stacked).astype(np.int32)
        data_vars[param_name] = (("regime", "crop"), stacked)

    # Attach declarative behavior variables along ("crop",)
    data_vars["seepage_rate_per_month"] = (
        ("crop",),
        np.asarray([cb.seepage_rate_per_month for cb in ordered_behaviors], dtype=np.float64),
    )
    data_vars["is_irrigated"] = (
        ("crop",),
        np.asarray([cb.is_irrigated for cb in ordered_behaviors], dtype=bool),
    )
    data_vars["is_open_water"] = (
        ("crop",),
        np.asarray([cb.is_open_water for cb in ordered_behaviors], dtype=bool),
    )
    data_vars["open_water_kc"] = (
        ("crop",),
        np.asarray([cb.open_water_kc for cb in ordered_behaviors], dtype=np.float64),
    )
    data_vars["riparian_full_seepage"] = (
        ("crop",),
        np.asarray([cb.riparian_full_seepage for cb in ordered_behaviors], dtype=bool),
    )
    data_vars["continuous_flooding"] = (
        ("crop",),
        np.asarray([cb.continuous_flooding for cb in ordered_behaviors], dtype=bool),
    )
    data_vars["pre_harvest_drain_days"] = (
        ("crop",),
        np.asarray([cb.pre_harvest_drain_days for cb in ordered_behaviors], dtype=np.int32),
    )
    data_vars["use_offseason_swdx_in_season"] = (
        ("crop",),
        np.asarray([cb.use_offseason_swdx_in_season for cb in ordered_behaviors], dtype=bool),
    )
    data_vars["year_type_unshifted_in_legacy"] = (
        ("crop",),
        np.asarray([cb.year_type_unshifted_in_legacy for cb in ordered_behaviors], dtype=bool),
    )
    data_vars["grain_cutoff_rule"] = (
        ("crop",),
        np.asarray([cb.grain_cutoff_rule for cb in ordered_behaviors], dtype=bool),
    )

    return xr.Dataset(
        data_vars=data_vars,
        coords={
            "regime": regimes,
            "crop": np.asarray(ordered_labels, dtype=object),
            "crop_code": ("crop", np.asarray(ordered_codes, dtype=object)),
        },
    )


def _parse_region_code(val: object) -> int:
    """Convert a region indicator (``0``/``1`` or ``'Lowlands'``/``'Uplands'``) to int ``0`` or ``1``."""
    s = str(val).strip().lower()
    if s in ("0", "0.0", "low", "lowland", "lowlands"):
        return 0
    if s in ("1", "1.0", "up", "upland", "uplands"):
        return 1
    raise ValueError(
        f"Unrecognized subarea REGION value: {val!r} (expected 0/Lowlands or 1/Uplands)"
    )


def read_subarea_params(
    filepath: str | Path,
    subareas: Iterable[int | str] | None = None,
    station_aliases: dict[str, str] | None = None,
    column_map: dict[str, str] | None = None,
) -> xr.Dataset:
    """Read subarea ETo spatial correction factors, lowland/upland region, and station weights.

    Aligns subareas by ID (never row order) and normalizes station weight
    columns via ``canonical_station_id`` (never alphabetical column position).
    Also handles historical header typos such as ``'ETo Corection Factor'`` and
    duplicate ``'REGION'`` / ``'REGION.1'`` columns.
    """
    df = pd.read_csv(filepath, index_col=False)
    area_col = _find_column(
        df.columns,
        ("area_id", "subarea", "subarea_id", "island"),
        column_map=column_map,
        context=f"et_correction CSV ({filepath})",
    )
    eto_col = _find_column(
        df.columns,
        ("eto_correction_factor", "etocorrectionfactor", "etocorectionfactor", "correction_factor"),
        column_map=column_map,
        context=f"et_correction CSV ({filepath})",
    )
    assert area_col is not None and eto_col is not None

    region_cols = [c for c in df.columns if normalize_token(c) in {"region", "region1", "uplow"}]
    if not region_cols:
        mapped_reg = _find_column(
            df.columns,
            ("region",),
            column_map=column_map,
            context=f"et_correction CSV ({filepath})",
        )
        assert mapped_reg is not None
        region_cols = [mapped_reg]

    ext_col = _find_column(df.columns, ("extension",), required=False)
    meta_cols = {area_col, eto_col, *region_cols}
    if ext_col is not None:
        meta_cols.add(ext_col)

    station_cols = [c for c in df.columns if c not in meta_cols]
    if not station_cols:
        raise ValueError(f"No precipitation station weight columns found in {filepath}")

    work = df.copy()
    work[area_col] = _coerce_subarea_id(work[area_col])
    if work[area_col].duplicated().any():
        dups = work.loc[work[area_col].duplicated(), area_col].tolist()
        raise ValueError(f"Duplicate subarea IDs in {filepath}: {dups}")

    work = work.set_index(area_col)
    if subareas is not None:
        target_ids = list(subareas)
        missing = [s for s in target_ids if s not in work.index]
        if missing:
            raise ValueError(f"Subareas {missing} missing from et_correction CSV {filepath}")
        work = work.reindex(target_ids)

    region_vals = np.asarray(
        [_parse_region_code(v) for v in work[region_cols[0]].tolist()], dtype=np.int32
    )
    eto_corr = work[eto_col].to_numpy(dtype=np.float64)

    canon_stations: list[str] = []
    weights_list: list[np.ndarray] = []
    seen_stations: dict[str, int] = {}
    for col in station_cols:
        st_id = canonical_station_id(col, custom_aliases=station_aliases)
        col_vals = work[col].to_numpy(dtype=np.float64)
        if st_id in seen_stations:
            weights_list[seen_stations[st_id]] = weights_list[seen_stations[st_id]] + col_vals
        else:
            seen_stations[st_id] = len(canon_stations)
            canon_stations.append(st_id)
            weights_list.append(col_vals)

    weights = np.stack(weights_list, axis=0)  # (station, subarea)

    return xr.Dataset(
        data_vars={
            "precip_weights": (
                ("station", "subarea"),
                weights,
                {"units": "1", "long_name": "Thiessen polygon station precipitation weight"},
            ),
            "eto_correction": (
                ("subarea",),
                eto_corr,
                {"units": "1", "long_name": "spatial ETo regression slope correction factor"},
            ),
            "region": (
                ("subarea",),
                region_vals,
                {"long_name": "Delta region code (0=Lowlands, 1=Uplands)"},
            ),
        },
        coords={
            "station": np.asarray(canon_stations, dtype=object),
            "subarea": np.asarray(work.index.tolist()),
        },
    )


def _extract_datetime_index(
    df: pd.DataFrame,
    filepath: str | Path,
    column_map: dict[str, str] | None = None,
) -> tuple[pd.DatetimeIndex, set[str]]:
    """Extract a clean ``DatetimeIndex`` from named date or (year, month, day) columns."""
    date_col = _find_column(
        df.columns,
        ("date", "datetime", "time"),
        column_map=column_map,
        required=False,
        context=str(filepath),
    )
    year_col = _find_column(df.columns, ("year",), column_map=column_map, required=False)
    month_col = _find_column(df.columns, ("month", "mon"), column_map=column_map, required=False)
    day_col = _find_column(df.columns, ("day",), column_map=column_map, required=False)
    doy_col = _find_column(df.columns, ("doy", "dayofyear"), column_map=column_map, required=False)

    used_cols: set[str] = {
        c for c in (date_col, year_col, month_col, day_col, doy_col) if c is not None
    }

    if date_col is not None:
        dt = pd.to_datetime(df[date_col], errors="raise")
    elif year_col is not None and month_col is not None and day_col is not None:
        dt = pd.to_datetime(
            {
                "year": pd.to_numeric(df[year_col], errors="raise"),
                "month": pd.to_numeric(df[month_col], errors="raise"),
                "day": pd.to_numeric(df[day_col], errors="raise"),
            }
        )
    elif year_col is not None and doy_col is not None:
        years = pd.to_numeric(df[year_col], errors="raise").astype(int)
        doys = pd.to_numeric(df[doy_col], errors="raise").astype(int)
        dt = pd.to_datetime(years.astype(str) + "-01-01") + pd.to_timedelta(doys - 1, unit="D")
    else:
        raise ValueError(
            f"Could not find a 'Date' column or ('year', 'month', 'day') columns in {filepath}. "
            f"Columns: {list(df.columns)}"
        )

    return pd.DatetimeIndex(dt, name="time"), used_cols


def read_station_precipitation(
    filepath: str | Path,
    start_date: str | pd.Timestamp | None = None,
    end_date: str | pd.Timestamp | None = None,
    station_aliases: dict[str, str] | None = None,
    column_map: dict[str, str] | None = None,
) -> xr.DataArray:
    """Read daily station precipitation CSV into ``xr.DataArray(time, station)`` in mm/d."""
    df = pd.read_csv(filepath, index_col=False)
    dt_index, time_cols = _extract_datetime_index(df, filepath, column_map=column_map)

    station_cols = [c for c in df.columns if c not in time_cols]
    if not station_cols:
        raise ValueError(f"No precipitation station columns found in {filepath}")

    work = df.loc[:, station_cols].copy()
    work.index = dt_index
    work = work.sort_index()
    if work.index.duplicated().any():
        raise ValueError(f"Duplicate timestamps in precipitation CSV: {filepath}")

    if start_date is not None or end_date is not None:
        s_ts = pd.Timestamp(start_date) if start_date is not None else work.index.min()
        e_ts = pd.Timestamp(end_date) if end_date is not None else work.index.max()
        work = work.loc[s_ts:e_ts]

    canon_names: list[str] = []
    arrays: list[np.ndarray] = []
    seen: set[str] = set()
    for col in station_cols:
        st_id = canonical_station_id(col, custom_aliases=station_aliases)
        if st_id in seen:
            raise ValueError(
                f"Multiple precipitation columns map to station {st_id!r} in {filepath}: {station_cols}"
            )
        seen.add(st_id)
        canon_names.append(st_id)
        arrays.append(pd.to_numeric(work[col], errors="raise").to_numpy(dtype=np.float64))

    data = np.stack(arrays, axis=1)  # (time, station)
    return xr.DataArray(
        data,
        dims=("time", "station"),
        coords={"time": work.index.values, "station": np.asarray(canon_names, dtype=object)},
        name="station_precip",
        attrs={"units": "mm day-1", "long_name": "station daily precipitation"},
    )


WEATHER_VAR_CANDIDATES: dict[str, tuple[str, ...]] = {
    "tmax": ("tx(oc)", "tx", "tmax", "tmax_c", "temp_max"),
    "tmin": ("tn(oc)", "tn", "tmin", "tmin_c", "temp_min"),
    "precip": ("pcp(mm)", "pcp", "precip", "ppt", "prcp"),
    "solar_rad": ("rs", "rs(mjm-2d-1)", "solar_rad", "srad", "sol_rad"),
    "wind_speed": ("u2", "wind", "wind_speed", "ws", "u2(ms-1)"),
    "dewpoint": ("td", "tdew", "dewpoint", "td(oc)"),
    "rh_max": ("rhx", "rhmax", "rh_max"),
    "rh_min": ("rhn", "rhmin", "rh_min"),
    "et0": ("et0", "eto", "et_0", "eto(mm)"),
}


def read_weather_station_csv(
    filepath: str | Path,
    start_date: str | pd.Timestamp | None = None,
    end_date: str | pd.Timestamp | None = None,
    column_map: dict[str, str] | None = None,
) -> xr.Dataset:
    """Read daily station meteorological variables (e.g. Lodi Tx/Tn, CIMIS PM inputs) by header name."""
    df = pd.read_csv(filepath, index_col=False)
    dt_index, _ = _extract_datetime_index(df, filepath, column_map=column_map)

    work = df.copy()
    work.index = dt_index
    work = work.sort_index()
    if work.index.duplicated().any():
        raise ValueError(f"Duplicate timestamps in weather CSV: {filepath}")

    if start_date is not None or end_date is not None:
        s_ts = pd.Timestamp(start_date) if start_date is not None else work.index.min()
        e_ts = pd.Timestamp(end_date) if end_date is not None else work.index.max()
        work = work.loc[s_ts:e_ts]

    data_vars: dict[str, tuple[tuple[str], np.ndarray]] = {}
    for canonical_var, candidates in WEATHER_VAR_CANDIDATES.items():
        col = _find_column(
            work.columns,
            candidates,
            column_map=column_map,
            required=False,
            context=str(filepath),
        )
        if col is not None:
            data_vars[canonical_var] = (
                ("time",),
                pd.to_numeric(work[col], errors="raise").to_numpy(dtype=np.float64),
            )

    if not data_vars:
        raise ValueError(
            f"No recognized meteorological variables found in {filepath}. Columns: {list(df.columns)}"
        )

    return xr.Dataset(data_vars=data_vars, coords={"time": work.index.values})


def write_weather_netcdf(
    ds_meteo: xr.Dataset,
    precip_output_path: str | Path | None,
    et_output_path: str | Path | None,
) -> None:
    """Write ``precip`` and ``ET0`` NetCDF files matching the DeltaCD schema."""
    if precip_output_path is not None:
        p_path = Path(precip_output_path)
        p_path.parent.mkdir(parents=True, exist_ok=True)
        da_pcp = ds_meteo["precip"].transpose("time", "subarea").copy()
        da_pcp.name = "precip"
        da_pcp.attrs = {"units": "mm"}
        da_pcp = da_pcp.assign_coords(subarea=da_pcp["subarea"].astype("i4"))
        da_pcp.to_netcdf(p_path)

    if et_output_path is not None:
        e_path = Path(et_output_path)
        e_path.parent.mkdir(parents=True, exist_ok=True)
        da_et0 = ds_meteo["et0"].transpose("time", "subarea").copy()
        da_et0.name = "ET0"
        da_et0.attrs = {"units": "mm"}
        da_et0 = da_et0.assign_coords(subarea=da_et0["subarea"].astype("i4"))
        da_et0.to_netcdf(e_path)


def write_detaw_netcdf(
    ds: xr.Dataset,
    output_path: str | Path,
    include_net_app: bool = False,
) -> xr.Dataset:
    """Write DETAW output Dataset to NetCDF matching downstream ``dcd.py`` contract."""
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    core_vars = ["et_c", "s_e", "precip", "et_aw", "d_sw", "e_r"]
    if include_net_app and "net_app" in ds.data_vars:
        core_vars.append("net_app")

    ds_out = ds[core_vars].transpose("time", "subarea", "crop").copy()
    # Drop extra helper coordinates (such as crop_code or water_year) from NetCDF if needed
    extra_coords = [c for c in ds_out.coords if c not in ("time", "subarea", "crop")]
    if extra_coords:
        ds_out = ds_out.drop_vars(extra_coords)

    try:
        ds_out = ds_out.assign_coords(subarea=ds_out["subarea"].astype("i4"))
    except (ValueError, TypeError):
        pass
    ds_out = ds_out.assign_coords(crop=ds_out["crop"].astype(str))

    ds_out.coords["subarea"].attrs = {"description": "subarea id"}
    ds_out.coords["crop"].attrs = {"description": "land-use category"}
    ds_out.attrs["title"] = (
        "Output from Delta Evapotranspiration of Applied Water (DETAW) module of DeltaCD"
    )
    created_at = pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")
    ds_out.attrs["history"] = (
        f"{created_at} - generated by DeltaCD version {deltacd.__version__}, DETAW module"
    )

    ds_out.to_netcdf(out_path)
    return ds_out
