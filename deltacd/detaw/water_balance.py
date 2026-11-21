"""Daily root-zone soil water balance and irrigation scheduling for DETAW.

Implements the daily water balance recurrence and management allowable
depletion ($SWD_x$) scheduling described in **Snyder et al. (2006), pp. 27-38**:

- Effective root depth ``erd = min(root_depth, max_soil_depth)``
- Plant available water ``PAW = erd * aw`` and yield threshold depletion
  ``YTD = PAW * (allowable_depletion_pct / 100)``
- Seasonal scheduling thresholds ``NA1`` (initial growth), ``NA2`` (mid/late
  growth), and ``NA3`` (end-of-season ``YTD``)
- Daily potential seepage ``Spg = seepage_rate_per_month / days_in_month * erd``
  in lowland subareas (``region == 0``), fixing both the Feb 29 leap-year bug
  (**Bug A1**) and the 1st-of-month ``Mon`` off-by-one bug (**Bug A2**).
- Compact ``@numba.njit(cache=True)`` sequential recurrence driven strictly by
  declarative crop-behavior flags (zero hardcoded crop index numbers).
"""

from __future__ import annotations

import numba
import numpy as np
import pandas as pd
import xarray as xr

from deltacd.detaw.et import _is_critical_year_type


def _compute_soil_and_schedule_arrays(
    time_index: pd.DatetimeIndex,
    da_etc_mm: xr.DataArray,
    da_ikc: xr.DataArray,
    ds_crop_params: xr.Dataset,
    ds_subarea: xr.Dataset,
    da_year_type: xr.DataArray,
    pad_year_type: str = "AN",
    season_mode: str = "legacy_wrapped",
    emulate_legacy_bugs: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Precompute daily arrays ``(spg_mm, swdx_mm, in_season_mask, rice_irrig_mask, fc_mm, os_swdx_2d)``.

    All daily arrays have shape ``(n_days, n_subareas, n_crops)`` in ``mm``;
    ``os_swdx_2d`` has shape ``(n_subareas, n_crops)``.
    """
    n_days = len(time_index)
    n_subareas = ds_subarea.sizes["subarea"]
    n_crops = ds_crop_params.sizes["crop"]

    regions = ds_subarea["region"].values.astype(np.int32)  # (n_subareas,) 0=Lowland, 1=Upland
    etc_vals = np.ascontiguousarray(
        da_etc_mm.transpose("time", "subarea", "crop").values, dtype=np.float64
    )
    if "subarea" in da_ikc.dims:
        ikc_3d = np.ascontiguousarray(
            da_ikc.transpose("time", "subarea", "crop").values, dtype=np.float64
        )
    else:
        ikc_2d = np.ascontiguousarray(da_ikc.transpose("time", "crop").values, dtype=np.float64)
        ikc_3d = np.broadcast_to(ikc_2d[:, np.newaxis, :], (n_days, n_subareas, n_crops))

    in_season_mask = ikc_3d != 0.0

    # Soil properties per (subarea, crop): in DETAW CSVs soil depth, root depth,
    # AW, and AD% are identical between critical and noncritical, though we read
    # from the noncritical table (matching legacy calc_kc_vals final year / regime).
    p_nc = ds_crop_params.sel(regime="noncritical")
    sdx = p_nc["max_soil_depth_mm"].values  # (n_crops,)
    rdx_lo = p_nc["root_depth_lowland_mm"].values
    rdx_up = p_nc["root_depth_upland_mm"].values
    aw_lo = p_nc["aw_lowland"].values
    aw_up = p_nc["aw_upland"].values
    adep = p_nc["allowable_depletion_pct"].values

    is_upland = regions[:, np.newaxis] == 1  # (n_subareas, 1)
    rdx = np.where(is_upland, rdx_up[np.newaxis, :], rdx_lo[np.newaxis, :])
    aw = np.where(is_upland, aw_up[np.newaxis, :], aw_lo[np.newaxis, :])
    erd = np.minimum(rdx, sdx[np.newaxis, :])  # (n_subareas, n_crops)
    paw = erd * aw
    ytd_2d = paw * (adep[np.newaxis, :] / 100.0)
    os_swdx_2d = aw * 300.0 * (adep[np.newaxis, :] / 100.0)
    fc_2d = (erd * aw * (adep[np.newaxis, :] / 100.0)) * 4.0  # FC = MaxSWD * 4 = 2 * PAW

    # Daily potential seepage Spg(time, subarea, crop)
    seep_rates = ds_crop_params["seepage_rate_per_month"].values  # (n_crops,)
    if emulate_legacy_bugs:
        # Reproduce Bug A1 (Feb=28 always) and Bug A2 (1st of Feb..Dec uses previous month!)
        # Note: On t_i == 0 (Oct 1, ii = 274 non-leap), legacy year_i=1 ran ii=1..273 first;
        # on ii=273 (Sep 30), ii > NI[8] (273 > 273) is False! So on Oct 1 (ii=274), Mon was
        # still 8 (September = 30 days) and only incremented to 9 AFTER computing Spg!
        num_days_legacy = np.array(
            [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=np.float64
        )
        ni_nonleap = np.array(
            [31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365], dtype=np.int32
        )
        ni_leap = np.array(
            [31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366], dtype=np.int32
        )
        denom_days = np.empty(n_days, dtype=np.float64)
        cur_mon = 0
        for t_i, ts in enumerate(time_index):
            doy = ts.dayofyear
            if doy == 1:
                cur_mon = 0
            elif t_i == 0:
                cur_mon = 8  # Sep (30 days) on Oct 1!
            if cur_mon < 12:
                denom_days[t_i] = num_days_legacy[cur_mon]
            else:
                denom_days[t_i] = 31.0
            thresh = ni_leap[cur_mon] if (ts.year % 4 == 0) else ni_nonleap[cur_mon]
            if doy > thresh:
                cur_mon += 1
    else:
        denom_days = time_index.days_in_month.to_numpy(dtype=np.float64)

    spg_mm = (
        seep_rates[np.newaxis, np.newaxis, :]
        / denom_days[:, np.newaxis, np.newaxis]
        * erd[np.newaxis, :, :]
    )
    # Uplands (region == 1) have zero channel seepage
    spg_mm[:, regions == 1, :] = 0.0

    # Precompute SWDx(time, subarea, crop) and rice_irrig_mask(time, crop)
    swdx_mm = np.broadcast_to(os_swdx_2d[np.newaxis, :, :], (n_days, n_subareas, n_crops)).copy()
    rice_irrig_mask = np.zeros((n_days, n_crops), dtype=bool)
    fc_mm = np.broadcast_to(fc_2d[np.newaxis, :, :], (n_days, n_subareas, n_crops)).copy()

    wy_lookup: dict[int, str] = {
        int(wy): str(yt) for wy, yt in zip(da_year_type["water_year"].values, da_year_type.values)
    }

    cal_years = sorted(time_index.year.unique().tolist())
    start_cal_year = cal_years[0] - 1 if season_mode.lower() == "continuous" else cal_years[0]
    all_cal_years = list(range(start_cal_year, cal_years[-1] + 1))

    # Map (year, doy) -> integer position in time_index for fast vectorized lookup
    year_arr = time_index.year.to_numpy(dtype=np.int32)
    doy_arr = time_index.dayofyear.to_numpy(dtype=np.int32)
    min_yr = start_cal_year
    max_yr = cal_years[-1] + 2
    yr_doy_to_pos = np.full((max_yr - min_yr + 1, 368), -1, dtype=np.int32)
    for idx in range(n_days):
        yr_doy_to_pos[year_arr[idx] - min_yr, doy_arr[idx]] = idx

    unshifted_flags = ds_crop_params["year_type_unshifted_in_legacy"].values
    use_os_swdx_flags = ds_crop_params["use_offseason_swdx_in_season"].values
    cont_flood_flags = ds_crop_params["continuous_flooding"].values
    drain_days_arr = ds_crop_params["pre_harvest_drain_days"].values

    # Pre-extract regime parameter arrays (2, n_crops) to avoid xarray .sel() inside loops
    bd_reg = ds_crop_params["begin_doy"].values
    ed_reg = ds_crop_params["end_doy"].values
    f_reg = ds_crop_params["irrig_freq"].values
    ab_reg = ds_crop_params["ab_pct"].values

    for c_idx in range(n_crops):
        unshifted = bool(unshifted_flags[c_idx])
        use_os_swdx = bool(use_os_swdx_flags[c_idx])
        is_flood = bool(cont_flood_flags[c_idx])
        drain_days = int(drain_days_arr[c_idx])
        ytd_crop = ytd_2d[:, c_idx]

        for cy in all_cal_years:
            dpy = 366 if (cy % 4 == 0 and (cy % 100 != 0 or cy % 400 == 0)) else 365
            prev_cy = cy - 1
            dpy_prev = (
                366 if (prev_cy % 4 == 0 and (prev_cy % 100 != 0 or prev_cy % 400 == 0)) else 365
            )
            if emulate_legacy_bugs and cy == cal_years[0]:
                dpy_prev = 365

            if season_mode.lower() == "continuous":
                nc_bd = int(bd_reg[0, c_idx])
                wy_for_season = (cy + 1) if nc_bd >= 274 else cy
                yt = wy_lookup.get(wy_for_season, pad_year_type)
            else:
                wy_for_season = (cy + 1) if unshifted else cy
                yt = wy_lookup.get(wy_for_season, pad_year_type)

            reg_i = 1 if _is_critical_year_type(yt) else 0
            bd = int(bd_reg[reg_i, c_idx])
            ed = int(ed_reg[reg_i, c_idx])
            irrig_f = int(f_reg[reg_i, c_idx]) or 30
            ab = float(ab_reg[reg_i, c_idx])

            end_adj = ed + 365 if ed <= bd else ed
            ee = end_adj - 365 if end_adj > 365 else end_adj
            b_day = int(0.01 * ab * (end_adj - bd) + bd)
            num_i1 = int((b_day - bd) / irrig_f)
            bi = bd + irrig_f * num_i1 if num_i1 != 0 else b_day

            season_days = np.arange(bd, end_adj + 1, dtype=np.int32)

            if is_flood:
                rice_end = end_adj - drain_days
                if season_mode.lower() == "continuous":
                    yrs_s = np.where(season_days > dpy, cy + 1, cy)
                    doys_s = np.where(season_days > dpy, season_days - dpy, season_days)
                    pos_s = yr_doy_to_pos[yrs_s - min_yr, doys_s]
                    valid = pos_s >= 0
                    rice_irrig_mask[pos_s[valid], c_idx] = season_days[valid] <= rice_end
                else:
                    doys_s = np.where(season_days > dpy, season_days - dpy, season_days)
                    in_bounds = (doys_s >= 1) & (doys_s <= dpy)
                    pos_s = yr_doy_to_pos[cy - min_yr, np.clip(doys_s, 0, 367)]
                    valid = in_bounds & (pos_s >= 0)
                    rice_irrig_mask[pos_s[valid], c_idx] = doys_s[valid] <= rice_end

            if use_os_swdx:
                continue

            # Vectorized index lookup for season_days
            if season_mode.lower() == "continuous":
                yrs_c = np.where(season_days > dpy, cy + 1, cy)
                doys_c = np.where(season_days > dpy, season_days - dpy, season_days)
                pos_c = yr_doy_to_pos[yrs_c - min_yr, np.clip(doys_c, 0, 367)]
                valid_c = (doys_c >= 1) & (pos_c >= 0)
                season_etc_cetc = np.zeros((len(season_days), n_subareas), dtype=np.float64)
                season_etc_cetc[valid_c, :] = etc_vals[pos_c[valid_c], :, c_idx]
                season_etc_et12 = season_etc_cetc
            else:
                if emulate_legacy_bugs:
                    jj_cetc = np.where(season_days > dpy_prev, season_days - dpy, season_days)
                    jj_et12 = np.where(season_days > dpy_prev, season_days - dpy_prev, season_days)
                else:
                    jj_cetc = np.where(season_days > dpy, season_days - dpy, season_days)
                    jj_et12 = jj_cetc

                pos_c = yr_doy_to_pos[cy - min_yr, np.clip(jj_cetc, 0, 367)]
                valid_c = (jj_cetc >= 1) & (jj_cetc <= dpy) & (pos_c >= 0)
                season_etc_cetc = np.zeros((len(season_days), n_subareas), dtype=np.float64)
                season_etc_cetc[valid_c, :] = etc_vals[pos_c[valid_c], :, c_idx]

                if emulate_legacy_bugs and dpy != dpy_prev:
                    pos_e = yr_doy_to_pos[cy - min_yr, np.clip(jj_et12, 0, 367)]
                    valid_e = (jj_et12 >= 1) & (jj_et12 <= dpy) & (pos_e >= 0)
                    season_etc_et12 = np.zeros((len(season_days), n_subareas), dtype=np.float64)
                    season_etc_et12[valid_e, :] = etc_vals[pos_e[valid_e], :, c_idx]
                else:
                    season_etc_et12 = season_etc_cetc

            cetc = season_etc_cetc.sum(axis=0)
            et1_steps = bi - bd + 1
            et1 = season_etc_et12[:et1_steps, :].sum(axis=0)
            na1 = (et1 / num_i1) if num_i1 != 0 else np.zeros(n_subareas, dtype=np.float64)

            na2 = np.zeros(n_subareas, dtype=np.float64)
            na3 = np.zeros(n_subareas, dtype=np.float64)
            li_arr = np.zeros(n_subareas, dtype=np.int32)

            if len(season_days) > et1_steps:
                cum_et2 = np.cumsum(season_etc_et12[et1_steps:, :], axis=0)
                rem_days = season_days[et1_steps:]
                cond_mat = (cum_et2 + et1[np.newaxis, :]) > (cetc - ytd_crop)[np.newaxis, :]
                any_exceed = cond_mat.any(axis=0)
                first_ks = np.argmax(cond_mat, axis=0)
                for s_idx in np.flatnonzero(any_exceed):
                    fk = int(first_ks[s_idx])
                    et2_val = float(cum_et2[fk, s_idx])
                    li_arr[s_idx] = int(rem_days[fk])
                    ni2 = int(et2_val / ytd_crop[s_idx]) + 1 if ytd_crop[s_idx] > 0 else 1
                    na2[s_idx] = et2_val / ni2
                    if num_i1 == 0:
                        na1[s_idx] = na2[s_idx]
                    na3[s_idx] = ytd_crop[s_idx]

            # Vectorized assignment of daily SWDx across the season / calendar year
            if season_mode.lower() == "continuous":
                yrs_s = np.where(season_days > dpy, cy + 1, cy)
                doys_s = np.where(season_days > dpy, season_days - dpy, season_days)
                pos_s = yr_doy_to_pos[yrs_s - min_yr, np.clip(doys_s, 0, 367)]
                valid_s = (doys_s >= 1) & (pos_s >= 0)
                sd_v = season_days[valid_s][:, np.newaxis]
                p_v = pos_s[valid_s]
                block = np.broadcast_to(na1[np.newaxis, :], (len(p_v), n_subareas)).copy()
                block = np.where(sd_v >= bi, na2[np.newaxis, :], block)
                block = np.where(
                    (li_arr[np.newaxis, :] > 0) & (sd_v >= li_arr[np.newaxis, :]),
                    na3[np.newaxis, :],
                    block,
                )
                block = np.where(sd_v == end_adj, na3[np.newaxis, :], block)
                swdx_mm[p_v, :, c_idx] = block
            else:
                year_doys = np.arange(1, dpy + 1, dtype=np.int32)
                pos_y = yr_doy_to_pos[cy - min_yr, year_doys]
                valid_y = pos_y >= 0
                yd_v = year_doys[valid_y][:, np.newaxis]
                p_y = pos_y[valid_y]
                block = swdx_mm[p_y, :, c_idx].copy()
                block = np.where(yd_v >= bd, na1[np.newaxis, :], block)
                block = np.where(yd_v >= bi, na2[np.newaxis, :], block)
                block = np.where(yd_v >= li_arr[np.newaxis, :], na3[np.newaxis, :], block)
                block = np.where(yd_v == ee, na3[np.newaxis, :], block)
                if ee < bd:
                    block = np.where(yd_v < bi, na2[np.newaxis, :], block)
                swdx_mm[p_y, :, c_idx] = block

    # Off-season days always set SWDx = os_swdx_2d
    swdx_mm = np.where(in_season_mask, swdx_mm, os_swdx_2d[np.newaxis, :, :])
    return spg_mm, swdx_mm, in_season_mask, rice_irrig_mask, fc_mm, os_swdx_2d


@numba.njit(cache=True)
def _water_balance_recurrence_kernel(
    etc_mm: np.ndarray,
    pcp_mm: np.ndarray,
    spg_mm: np.ndarray,
    swdx_mm: np.ndarray,
    in_season_mask: np.ndarray,
    rice_irrig_mask: np.ndarray,
    os_swdx_2d: np.ndarray,
    regions: np.ndarray,
    is_irrigated: np.ndarray,
    is_open_water: np.ndarray,
    riparian_full_seepage: np.ndarray,
    continuous_flooding: np.ndarray,
    offseason_swd_cap: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sequential daily root-zone soil water balance recurrence across ``(time, subarea, crop)``.

    Returns ``(espg_mm, ern_mm, dsw_mm, net_app_mm, swd_mm)`` in ``mm``.
    """
    n_days, n_subareas, n_crops = etc_mm.shape
    espg_out = np.empty((n_days, n_subareas, n_crops), dtype=np.float64)
    ern_out = np.empty((n_days, n_subareas, n_crops), dtype=np.float64)
    dsw_out = np.empty((n_days, n_subareas, n_crops), dtype=np.float64)
    na_out = np.empty((n_days, n_subareas, n_crops), dtype=np.float64)
    swd_out = np.empty((n_days, n_subareas, n_crops), dtype=np.float64)

    for s in range(n_subareas):
        is_upland = regions[s] == 1
        for c in range(n_crops):
            irrig = is_irrigated[c]
            open_w = is_open_water[c]
            riparian = riparian_full_seepage[c]
            flood = continuous_flooding[c]
            os_cap = os_swdx_2d[s, c]

            swd = 0.0
            for t in range(n_days):
                etc = etc_mm[t, s, c]
                pcp = pcp_mm[t, s]
                spg = spg_mm[t, s, c]
                in_seas = in_season_mask[t, s, c]

                # 1. Effective seepage (Espg)
                if (swd + etc - spg) >= 0.0:
                    espg = spg
                else:
                    espg = 0.0 if is_upland else (swd + etc)

                if flood and in_seas:
                    espg = 0.0
                if riparian:
                    espg = 0.0 if is_upland else etc
                if open_w:
                    espg = 0.0

                # 2. Effective rainfall (ERn)
                if (swd + etc - espg - pcp) >= 0.0:
                    ern = pcp
                else:
                    ern = swd + etc - espg

                if flood and in_seas:
                    ern = 0.0
                if riparian:
                    ern = 0.0
                if open_w:
                    ern = pcp

                # 3. Daily change in soil water depletion (Dsw)
                if not open_w:
                    dsw = etc - ern - espg
                else:
                    dsw = etc - ern + espg

                # 4. Net irrigation application (NetApp)
                net_app = 0.0
                if in_seas and irrig:
                    if flood:
                        if rice_irrig_mask[t, c]:
                            net_app = swd + dsw
                        else:
                            net_app = 0.0
                    else:
                        if (swd + dsw) > swdx_mm[t, s, c]:
                            net_app = swd + dsw

                # 5. Update soil water depletion (SWD)
                new_swd = swd + dsw - net_app
                if offseason_swd_cap and (not in_seas) and dsw > 0.0:
                    # Snyder et al. (2006, pp. 30, 37): in off-season, SWD cannot
                    # increase beyond os_cap, though rainfall/seepage (dsw < 0)
                    # can decrease it.
                    cap_limit = swd if swd > os_cap else os_cap
                    if new_swd > cap_limit:
                        new_swd = cap_limit

                if new_swd < 0.0:
                    new_swd = 0.0
                swd = new_swd

                espg_out[t, s, c] = espg
                ern_out[t, s, c] = ern
                dsw_out[t, s, c] = dsw
                na_out[t, s, c] = net_app
                swd_out[t, s, c] = swd

    return espg_out, ern_out, dsw_out, na_out, swd_out


def run_soil_water_balance(
    ds_et: xr.Dataset,
    ds_crop_params: xr.Dataset,
    ds_subarea: xr.Dataset,
    da_year_type: xr.DataArray,
    pad_year_type: str = "AN",
    season_mode: str = "legacy_wrapped",
    offseason_swd_cap: bool = False,
    emulate_legacy_bugs: bool = False,
) -> xr.Dataset:
    """Execute the daily root-zone soil water balance over ``ds_et["time"]``.

    Note: In legacy DETAW, ``year_i = 1`` starts the water balance loop on
    January 1 of ``model_start_year`` (initializing ``SWD = 0`` on Jan 1 and
    resetting ``SWD = 0`` again on Oct 1, ``ii = 274``!), so starting the daily
    ``SWD`` recurrence on October 1 of ``model_start_year`` with ``SWD = 0`` is
    mathematically identical to legacy DETAW while avoiding 273 wasted days!
    """
    time_index = pd.DatetimeIndex(ds_et["time"].values)

    spg_mm, swdx_mm, in_season_mask, rice_irrig_mask, fc_mm, os_swdx_2d = (
        _compute_soil_and_schedule_arrays(
            time_index=time_index,
            da_etc_mm=ds_et["et_c_mm"],
            da_ikc=ds_et["ikc"],
            ds_crop_params=ds_crop_params,
            ds_subarea=ds_subarea,
            da_year_type=da_year_type,
            pad_year_type=pad_year_type,
            season_mode=season_mode,
            emulate_legacy_bugs=emulate_legacy_bugs,
        )
    )

    etc_arr = np.ascontiguousarray(
        ds_et["et_c_mm"].transpose("time", "subarea", "crop").values, dtype=np.float64
    )
    pcp_arr = np.ascontiguousarray(
        ds_et["precip"].transpose("time", "subarea").values, dtype=np.float64
    )
    regions = np.ascontiguousarray(ds_subarea["region"].values, dtype=np.int32)

    espg_mm, ern_mm, dsw_mm, net_app_mm, swd_mm = _water_balance_recurrence_kernel(
        etc_arr,
        pcp_arr,
        spg_mm,
        swdx_mm,
        in_season_mask,
        rice_irrig_mask,
        os_swdx_2d,
        regions,
        np.ascontiguousarray(ds_crop_params["is_irrigated"].values, dtype=np.bool_),
        np.ascontiguousarray(ds_crop_params["is_open_water"].values, dtype=np.bool_),
        np.ascontiguousarray(ds_crop_params["riparian_full_seepage"].values, dtype=np.bool_),
        np.ascontiguousarray(ds_crop_params["continuous_flooding"].values, dtype=np.bool_),
        bool(offseason_swd_cap),
    )

    swc_mm = fc_mm - swd_mm
    coords = {
        "time": time_index,
        "subarea": ds_et["subarea"],
        "crop": ds_crop_params["crop"],
        "crop_code": ds_crop_params["crop_code"],
    }
    dims = ("time", "subarea", "crop")

    return xr.Dataset(
        data_vars={
            "et_c_mm": (dims, etc_arr, {"units": "mm day-1"}),
            "precip_mm": (("time", "subarea"), pcp_arr, {"units": "mm day-1"}),
            "spg_mm": (dims, spg_mm, {"units": "mm day-1"}),
            "s_e_mm": (dims, espg_mm, {"units": "mm day-1"}),
            "e_r_mm": (dims, ern_mm, {"units": "mm day-1"}),
            "d_sw_mm": (dims, dsw_mm, {"units": "mm day-1"}),
            "net_app_mm": (dims, net_app_mm, {"units": "mm day-1"}),
            "swd_mm": (dims, swd_mm, {"units": "mm"}),
            "swc_mm": (dims, swc_mm, {"units": "mm"}),
            "fc_mm": (dims, fc_mm, {"units": "mm"}),
            "swdx_mm": (dims, swdx_mm, {"units": "mm"}),
        },
        coords=coords,
    )
