===============================================================
DETAW (Delta Evapotranspiration of Applied Water) Overview
===============================================================

Delta Evapotranspiration of Applied Water (DETAW) simulates daily reference evapotranspiration (:math:`ET_0`), crop evapotranspiration (:math:`ET_c`), root-zone soil water balance, and evapotranspiration of applied water (:math:`ET_{aw}`) across subareas (islands) of the Sacramento–San Joaquin River Delta and Suisun Marsh.

Theoretical foundations are documented in **Snyder et al. (2006)**, *DETAW: Delta Evapotranspiration of Applied Water, Version 1.0 Final Report* (`DETAW v1.0 report <https://data.cnra.ca.gov/dataset/detaw>`_), the `DETAW v2.1 user manual <https://github.com/CADWRDeltaModeling/DETAW-DCD/blob/master/DETAW/Documents/DETAW%20v2.1_user's%20manual.pdf>`_, and DWR Bay-Delta Office Annual Report chapters.

Modular Architecture (`deltacd.detaw`)
======================================

In DeltaCD, ``deltacd.detaw`` is organized into modular ``xarray`` and ``pandas`` components with Numba-accelerated sequential daily kernels:

* ``deltacd.detaw.config``: Typed YAML configuration (backward-compatible with legacy ``detaw:`` configs plus extensible ``meteorology:``, ``et:``, ``options:``, and ``column_map:`` sections).
* ``deltacd.detaw.params``: Declarative ``CropBehavior`` metadata and station/crop alias registries (supporting Delta 15-crop, Suisun 16-crop with Duck Pond, and custom land-use categories without hardcoded indices).
* ``deltacd.detaw.io``: Order-independent CSV/NetCDF readers and CF-styled NetCDF writers (columns, rows, stations, subareas, and crops are matched by normalized header/label names rather than fixed positions).
* ``deltacd.detaw.meteorology``: Pluggable ``MeteorologyProvider`` implementations and vectorized :math:`ET_0` equations:

  - **Hargreaves–Samani** (``hargreaves_samani``, default; Snyder et al., 2006, Appendix B, Eq. 23)
  - **ASCE Standardized Penman–Monteith** (``asce_penman_monteith``; Snyder et al., 2006, Appendix B, Eq. 5–26)
  - **Precomputed / Gridded Datasets** (``DatasetMeteorologyProvider`` for NetCDF or CSV ``et0`` and ``precip``)

* ``deltacd.detaw.et``: Pluggable ``ETcProvider`` implementations:

  - ``KcETcProvider`` (default): Two-stage bare-soil evaporation :math:`OK_c` (Ventura et al., 2006) combined with critical/non-critical seasonal :math:`IK_c` curves (:math:`ET_c = K_c \times ET_0`).
  - ``ExternalETProvider`` (**OpenET integration**): Ingests daily or monthly external :math:`ET_c` or :math:`ET_{aw}` from NetCDF or CSV (with uniform or :math:`ET_0`-scaled monthly disaggregation).

* ``deltacd.detaw.water_balance`` & ``deltacd.detaw.etaw``: Daily root-zone soil moisture depletion (:math:`SWD`), effective seepage (:math:`S_e`), effective rainfall (:math:`E_r`), net irrigation applications (:math:`NA`), and smoothed water-year :math:`ET_{aw}` for downstream DCD routing.
* ``notebooks/``: Educational Jupyter notebooks walking through DETAW theory (``01_detaw_theory_and_walkthrough.ipynb``) and modular :math:`ET_0` / OpenET extensions (``02_detaw_modular_et_and_openet.ipynb``).

Main Input File
===============

The main input file of DETAW uses YAML format. All relative paths inside the YAML file are resolved relative to the directory containing the YAML file:

..  code-block:: yaml
    :caption: detaw_dsm2.yaml

    detaw:
      start_water_year: 1922
      end_water_year: 2024
      input_pcp: 'inputs/mm_pcp.csv'
      input_temperature: 'inputs/LODI_PT.csv'
      landuse: 'inputs/landuse_dsm2.csv'
      et_correction: 'inputs/et_correction_dsm2.csv'
      critical: 'inputs/crop_info_critical.csv'
      noncritical: 'inputs/crop_info_noncritical.csv'

      # Optional modular meteorology configuration
      meteorology:
        provider: 'station'              # 'station' (default) or 'dataset'
        eto_method: 'hargreaves_samani'  # 'hargreaves_samani', 'asce_penman_monteith', or 'direct'
        latitude_deg: 38.5
        temp_qc_mode: 'clamp'            # 'clamp' (legacy) or 'swap'

      # Optional modular ET / OpenET configuration
      et:
        provider: 'kc_eto'               # 'kc_eto' (default) or 'openet'
        # path: 'inputs/openet_etc.nc'
        # variable: 'et_c'
        # units: 'mm/d'                  # 'mm/d', 'mm/month', 'in/d', 'in/month'
        # temporal_resolution: 'daily'   # 'daily' or 'monthly'
        # disaggregation: 'et0_scaled'   # 'et0_scaled' or 'uniform'

      # Optional physics switches
      options:
        offseason_swd_cap: false
        include_net_app_output: false

      # Model output control
      detaw_output: 'outputs/detawoutput_dsm2.nc'
      precip_output: 'outputs/precip_dsm2.nc'
      et_output: 'outputs/ET0.nc'
      daily_output: 1
      daily_output_unit: 1

Core Configuration Keys
-----------------------

* ``start_water_year`` / ``end_water_year``: Simulation period in California water years (Oct 1 of ``start_water_year - 1`` through Sep 30 of ``end_water_year``).
* ``input_pcp``: Daily station precipitation CSV (mm/d). Station columns may appear in any order and are matched to ``et_correction`` weights by normalized station name/alias.
* ``input_temperature``: Daily station weather CSV (``Tx(oC)``, ``Tn(oC)`` for Hargreaves–Samani, plus ``solar_rad``, ``wind_speed``, and ``dewpoint`` or ``rh_max``/``rh_min`` when using ``asce_penman_monteith``). If ``(start_water_year - 1)-09-30`` is present, it is used automatically as a 1-day bare-soil evaporation spin-up; if data begin on ``10-01``, DETAW runs cleanly without requiring a spin-up row.
* ``landuse``: Annual land-use CSV in hectares containing ``area_id`` (or ``SubArea``), ``DATE`` (water year), ``TYPE`` (water-year type: ``C``/``D`` critical vs. ``W``/``AN``/``BN`` non-critical), and crop columns in any order.
* ``et_correction``: Subarea ETo spatial correction factors, ``REGION`` (``0``/``Lowlands`` vs. ``1``/``Uplands``), and precipitation station weights.
* ``critical`` / ``noncritical``: Crop and soil parameter CSVs (supports both legacy transposed ``crop_info_*.csv`` and tidy row-per-crop CSV tables).
* ``crop_behavior`` *(optional)*: Custom crop hydrological behavior CSV for introducing new land-use categories.
* ``detaw_output``: Output NetCDF path containing ``(time, subarea, crop)`` daily variables ``et_c``, ``s_e``, ``precip``, ``et_aw``, ``d_sw``, ``e_r`` (and ``net_app`` when ``include_net_app_output: true``) in ``acre-foot day-1`` (``daily_output_unit: 1``) or ``mm day-1`` (``daily_output_unit: 0``).
* ``precip_output`` / ``et_output`` *(optional)*: Output NetCDF paths for daily subarea precipitation and :math:`ET_0` (mm/d) over the simulation period.

Input files
-----------
Below are sample input file formats.

- *input_pcp*: Input precipitation sample format

  CSV file containing station precipitation in mm.

..  code-block:: text
    :caption: input_pcp

      year,month,day,DOY,Brentwood,Davis,Galt,Lodi,RioVista,Stockton,Tracy
      1921,9,30,273,0.0,0.0,0.0,0.0,0.0,0.0,0.0
      1921,10,1,274,0.0,0.0,0.0,0.0,0.0,0.0,0.0
      ....
      ....
      2024,9,29,273,0.0,0.5,0.0,0.0,0.0,0.0,0.0
      2024,9,30,274,0.0,0.0,0.0,0.0,0.0,0.0,0.0

- *input_temperature*: Input temperature sample format

  CSV file containing station temperature in deg C.

..  code-block:: text
    :caption: input_temperature

      Date,Year,Month,DOY,Pcp(mm),Tx(oC),Tn(oC)
      9/30/1921,1921,9,273,0,27.8,9.4
      10/1/1921,1921,10,274,0,27.8,9.4
      ....
      ....
      9/29/2024,2024,9,273,0,20.9,6.5
      9/30/2024,2024,9,274,0,23.6,4

- *landuse*: Landuse sample format

  CSV file containing landuse for each area_id for different years.

..  code-block:: text
    :caption: landuse

      DATE,TYPE,UR,PA,AL,FI,SB,GR,RI,TR,TO,OR,VI,RV,NV,DGR,WS,area_id
      1922,AN,13.0,75.0,387.0,1070.0,147,879.0,0,1642.0,0.0,227.0,0.0,36,178,0,141,1
      1923,BN,13.0,90.0,385.0,1082.0,160,825.0,0,1459.0,0.0,236.0,0.0,49,302,0,195,1
      ....
      ....
      2023,BN,0.0,113.032,28.91,194.648,0,7.11,0,63.0,470.718,389.27,271.25,0,0,0,0,174
      2024,W,0.0,113.032,28.91,194.648,0,7.11,0,63.0,470.718,389.27,271.25,0,0,0,0,174


Output files
------------

- *detaw_output*: DETAW output format

  After successful detaw model run a netCDF file containing output is created whereever the detaw_output points. Below is a sample of the what the output netCDF header might look like.

..  code-block:: text
    :caption: detaw_output

    detawoutput_dsm2 {
      dimensions:
        subarea = 168 ;
        crop = 15 ;
        time = 37621 ;
      variables:
        int subarea(subarea) ;
                subarea:description = "subarea id" ;
        string crop(crop) ;
                crop:description = "land-use category" ;
        int64 time(time) ;
                time:units = "days since 1921-10-01 00:00:00" ;
                time:calendar = "proleptic_gregorian" ;
        double et_c(time, subarea, crop) ;
                et_c:_FillValue = NaN ;
                et_c:long_name = "crop evapotranspiration" ;
                et_c:units = "acre-foot day-1" ;
        double s_e(time, subarea, crop) ;
                s_e:_FillValue = NaN ;
                s_e:long_name = "Effective seepage" ;
                s_e:units = "acre-foot day-1" ;
        double precip(time, subarea, crop) ;
                precip:_FillValue = NaN ;
                precip:long_name = "precipitation" ;
                precip:units = "acre-foot day-1" ;
        double et_aw(time, subarea, crop) ;
                et_aw:_FillValue = NaN ;
                et_aw:long_name = "evapotranspiration of applied water" ;
                et_aw:units = "acre-foot day-1" ;
        double d_sw(time, subarea, crop) ;
                d_sw:_FillValue = NaN ;
                d_sw:long_name = "change in soil water content" ;
                d_sw:units = "acre-foot day-1" ;
        double e_r(time, subarea, crop) ;
                e_r:_FillValue = NaN ;
                e_r:long_name = "effective rainfall" ;
                e_r:units = "acre-foot day-1" ;
      }

- *precip_output*: Precipitation output format

  After successful detaw model run a netCDF file containing output is created whereever the precip_output points. Below is a sample of the what the output netCDF header might look like.

..  code-block:: text
    :caption: precip_output

    netcdf precip_dsm2 {
      dimensions:
        area_id = 168 ;
        time = 37621 ;
      variables:
        int64 time(time) ;
                time:units = "days since 1921-10-01 00:00:00" ;
                time:calendar = "proleptic_gregorian" ;
        int subarea(subarea) ;
        double precip(time, subarea) ;
                precip:_FillValue = NaN ;
                precip:units = "mm" ;
      }

- *et_output*: ET output format

  After successful detaw model run a netCDF file containing output is created whereever the et_output points. Below is a sample of the what the output netCDF header might look like.

..  code-block:: text
    :caption: et_output

    netcdf ET0 {
      dimensions:
        time = 37621 ;
        subarea = 168 ;
      variables:
        int64 time(time) ;
                time:units = "days since 1921-10-01 00:00:00" ;
                time:calendar = "proleptic_gregorian" ;
        int subarea(subarea) ;
        double ET0(subarea, time) ;
                ET0:_FillValue = NaN ;
                ET0:units = "mm" ;
      }
