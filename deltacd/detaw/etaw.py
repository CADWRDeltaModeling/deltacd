"""Smoothed daily Evapotranspiration of Applied Water ($ET_{aw}$) post-processing.

In DETAW v2.0 / DeltaCD, downstream ``dcd.py`` consumes a smoothed daily
``et_aw`` rather than spiky discrete field-scale irrigation events (``net_app``).
For each water year (Oct 1 - Sep 30), subarea, and crop:

1. Compute the annual change in root-zone soil water content ``SumDelSWC`` and
   distribute ``MonthlySWC = SumDelSWC / 12`` across each month's days:
   ``et_aw_raw[d] = d_sw[d] - MonthlySWC / days_in_month[d]``.
2. Within each month of the water year, iteratively zero negative daily values
   and redistribute their negative sum equally across positive days in the same
   month until no negative days remain (or the month becomes all zero).
3. Zero ``et_aw`` after the last irrigation day of the season (using
   ``cutoff_inclusive=True`` for legacy ``DOY >= DOYLIrrig`` baseline, or
   ``False`` for ``DOY > DOYLIrrig``), including the split fall/spring cutoff
   rule for fall-planted irrigated grain.
"""

from __future__ import annotations

import numba
import numpy as np
import pandas as pd
import xarray as xr


@numba.njit(cache=True)
def _smoothed_etaw_kernel(
    dsw_mm: np.ndarray,
    swc_mm: np.ndarray,
    fc_mm: np.ndarray,
    net_app_mm: np.ndarray,
    wy_starts: np.ndarray,
    wy_ends: np.ndarray,
    months: np.ndarray,
    doys: np.ndarray,
    cal_years: np.ndarray,
    days_in_month_clean: np.ndarray,
    grain_cutoff_flags: np.ndarray,
    cutoff_inclusive: bool,
    emulate_legacy_bugs: bool,
) -> np.ndarray:
    """Compute smoothed daily ``et_aw_mm(time, subarea, crop)`` via Numba."""
    n_days, n_subareas, n_crops = dsw_mm.shape
    n_wy = len(wy_starts)
    etaw_out = np.zeros((n_days, n_subareas, n_crops), dtype=np.float64)

    num_days_legacy = np.array([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=np.int32)

    for s in range(n_subareas):
        for c in range(n_crops):
            is_grain_cutoff = grain_cutoff_flags[c]

            # Track last irrigation DOY per IrrigYear across the simulation
            # Note: In legacy main_calc_loop, IrrigYear starts at 0 on WY1 Oct 1
            # and increments by 1 every Jan 1 (DOY == 1).
            # Thus WY index w (0-based, starting Oct 1 of cal_year Y0):
            # - Oct..Dec of WY w has IrrigYear = w
            # - Jan..Sep of WY w has IrrigYear = w + 1
            # Meanwhile in legacy calc_etaw_daily:
            # yDaily[dd] = year_i + model_start_year - 1:
            # - Oct..Dec of WY w is in year_i = w + 1 -> yDaily - model_start_year = w!
            # - Jan..Sep of WY w is in year_i = w + 2 -> yDaily - model_start_year = w + 1!
            # Both map to exact calendar-year index (cal_year - model_start_year)!
            base_cal_year = cal_years[0]
            n_cal_slots = cal_years[-1] - base_cal_year + 2
            doy_last_irrig = np.zeros(n_cal_slots, dtype=np.float64)
            doy_grain_last_irrig = np.zeros(n_cal_slots, dtype=np.float64)

            for t in range(n_days):
                if net_app_mm[t, s, c] > 1.0e-6:
                    cy_slot = cal_years[t] - base_cal_year
                    doy_val = float(doys[t])
                    doy_last_irrig[cy_slot] = doy_val
                    if is_grain_cutoff and doy_val < 152.0:
                        doy_grain_last_irrig[cy_slot] = doy_val

            prev_swc = 0.0
            for w in range(n_wy):
                i_start = wy_starts[w]
                i_end = wy_ends[w]  # exclusive

                # 1. Annual soil-water content change SumDelSWC across the water year
                sum_del_swc = 0.0
                for t in range(i_start, i_end):
                    if w == 0 and t == i_start:
                        if emulate_legacy_bugs:
                            prev_swc = int(fc_mm[t, s, c] * 10.0) / 10.0
                        else:
                            prev_swc = fc_mm[t, s, c]

                    if emulate_legacy_bugs:
                        swc_cmp = int(swc_mm[t, s, c] * 1000.0) / 1000.0
                        del_swc = prev_swc - swc_cmp
                        prev_swc = swc_mm[t, s, c]
                    else:
                        del_swc = prev_swc - swc_mm[t, s, c]
                        prev_swc = swc_mm[t, s, c]

                    sum_del_swc += del_swc

                monthly_swc = sum_del_swc / 12.0

                # 2. Subtract monthly_swc / days_in_month and redistribute negatives within each month
                t_cursor = i_start
                while t_cursor < i_end:
                    cur_m = months[t_cursor]
                    m_end = t_cursor + 1
                    while m_end < i_end and months[m_end] == cur_m:
                        m_end += 1

                    m_len = m_end - t_cursor
                    if emulate_legacy_bugs:
                        denom_m = float(num_days_legacy[cur_m])
                    else:
                        denom_m = float(days_in_month_clean[t_cursor])

                    for t in range(t_cursor, m_end):
                        dsw_val = (
                            int(dsw_mm[t, s, c] * 1000.0) / 1000.0
                            if emulate_legacy_bugs
                            else dsw_mm[t, s, c]
                        )
                        etaw_out[t, s, c] = dsw_val - (monthly_swc / denom_m)

                    # Iterative negative redistribution within month [t_cursor, m_end)
                    # Note: On leap-year Feb in legacy DETAW, NumDay[2] == 28, so
                    # Feb 29 (etDay=29) was excluded from calc_etawavg's 1..28 loop!
                    active_len = int(denom_m) if emulate_legacy_bugs else m_len
                    if active_len > m_len:
                        active_len = m_len

                    keep_going = True
                    avg_neg = 0.0
                    while keep_going:
                        sum_neg = 0.0
                        count_neg = 0
                        for k in range(active_len):
                            val = etaw_out[t_cursor + k, s, c]
                            if val < 0.0:
                                sum_neg += val
                                count_neg += 1
                                etaw_out[t_cursor + k, s, c] = 0.0

                        keep_going = abs(sum_neg) > 0.0
                        if active_len > count_neg:
                            avg_neg = sum_neg / float(active_len - count_neg)
                        else:
                            avg_neg = 0.0

                        for k in range(active_len):
                            if etaw_out[t_cursor + k, s, c] > 0.0:
                                etaw_out[t_cursor + k, s, c] += avg_neg

                    if not emulate_legacy_bugs:
                        for k in range(m_len):
                            if etaw_out[t_cursor + k, s, c] < 0.0:
                                etaw_out[t_cursor + k, s, c] = 0.0

                    t_cursor = m_end

                # 3. Zero values after (or on/after) last irrigation day
                for t in range(i_start, i_end):
                    cy_slot = cal_years[t] - base_cal_year
                    doy_val = float(doys[t])
                    last_irr = doy_last_irrig[cy_slot]

                    if not is_grain_cutoff:
                        should_zero = (
                            (doy_val >= last_irr) if cutoff_inclusive else (doy_val > last_irr)
                        )
                        if should_zero:
                            etaw_out[t, s, c] = 0.0
                    else:
                        grain_last = doy_grain_last_irrig[cy_slot]
                        if last_irr < 275.0:
                            cond = (
                                (doy_val >= grain_last)
                                if cutoff_inclusive
                                else (doy_val > grain_last)
                            )
                            if cond and doy_val < 275.0:
                                etaw_out[t, s, c] = 0.0
                        else:
                            cond = (
                                (doy_val >= last_irr) if cutoff_inclusive else (doy_val > last_irr)
                            )
                            if cond:
                                etaw_out[t, s, c] = 0.0

    return etaw_out


def compute_smoothed_etaw(
    ds_wb: xr.Dataset,
    ds_crop_params: xr.Dataset,
    cutoff_inclusive: bool = True,
    emulate_legacy_bugs: bool = False,
) -> xr.DataArray:
    """Compute smoothed daily ``et_aw_mm(time, subarea, crop)`` in ``mm day-1``."""
    time_index = pd.DatetimeIndex(ds_wb["time"].values)
    water_years = np.where(time_index.month >= 10, time_index.year + 1, time_index.year)

    unique_wys = np.unique(water_years)
    wy_starts = np.empty(len(unique_wys), dtype=np.int64)
    wy_ends = np.empty(len(unique_wys), dtype=np.int64)
    for idx, wy in enumerate(unique_wys):
        positions = np.flatnonzero(water_years == wy)
        wy_starts[idx] = int(positions[0])
        wy_ends[idx] = int(positions[-1]) + 1

    dsw_arr = np.ascontiguousarray(
        ds_wb["d_sw_mm"].transpose("time", "subarea", "crop").values, dtype=np.float64
    )
    swc_arr = np.ascontiguousarray(
        ds_wb["swc_mm"].transpose("time", "subarea", "crop").values, dtype=np.float64
    )
    fc_arr = np.ascontiguousarray(
        ds_wb["fc_mm"].transpose("time", "subarea", "crop").values, dtype=np.float64
    )
    na_arr = np.ascontiguousarray(
        ds_wb["net_app_mm"].transpose("time", "subarea", "crop").values, dtype=np.float64
    )

    etaw_mm = _smoothed_etaw_kernel(
        dsw_arr,
        swc_arr,
        fc_arr,
        na_arr,
        wy_starts,
        wy_ends,
        time_index.month.to_numpy(dtype=np.int32),
        time_index.dayofyear.to_numpy(dtype=np.int32),
        time_index.year.to_numpy(dtype=np.int32),
        time_index.days_in_month.to_numpy(dtype=np.int32),
        np.ascontiguousarray(ds_crop_params["grain_cutoff_rule"].values, dtype=np.bool_),
        bool(cutoff_inclusive),
        bool(emulate_legacy_bugs),
    )

    return xr.DataArray(
        etaw_mm,
        dims=("time", "subarea", "crop"),
        coords={
            "time": time_index,
            "subarea": ds_wb["subarea"],
            "crop": ds_wb["crop"],
        },
        name="et_aw_mm",
        attrs={
            "units": "mm day-1",
            "long_name": "smoothed daily evapotranspiration of applied water depth",
        },
    )
