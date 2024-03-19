=============================
DCD (Delta Channel Depletion)
=============================

Delta Channel Depletion (DCD) is an extension of DETAW to estimate daily channel depletions, including diversions, drainages, and seepages of the Delta islands. DCD also distributes these estimates for other model inputs such as DSM2, SCHISM, and CalSIM.

This version of DCD is based on the previous version of DCD from DETAW-DCD, but it is totally rewritten in Python. The algorithms in and the calibration of the previous version DCD are described in an annual report chapter, `Calibrating and Validating Delta Channel Depletion Estimates <https://data.cnra.ca.gov/dataset/dcd/resource/24890484-11a6-4ada-a61a-3d1f7fdd948e>`_, and they are the same in the current version. The changes in the current DCD in DeltaCD are described in Chapter 3, Updating DETAW-DCD to DeltaCD, of `our 2023 annual report <https://data.cnra.ca.gov/dataset/methodology-for-flow-and-salinity-estimates-in-the-sacramento-san-joaquin-delta-and-suisun-marsh/resource/dcabdb20-e638-4cf5-b199-78e78f0d482f>`_.

Main input file
~~~~~~~~~~~~~~~

The main input file of DCD follows YAML format, just like DETAW. The following is an example of the DCD main input file. Currently, the parameters in the main input file keep the names from the previous DETAW-DCD. Note that some of the options are not well maintained, and they are subject to change.

..  code-block:: yaml
    :caption: dcd_dsm2.yaml

    dcd:
      path_subarea_info: "inputs/subarea_info_dsm2.csv"
      path_irrigation_efficiency: "inputs/irrigation_efficiencies_dsm2.csv"
      path_leach_applied: "inputs/leach_applied_dsm2.csv"
      path_leach_drained: "inputs/leach_drained_dsm2.csv"
      path_detaw_output: "output/detawoutput_dsm2.nc"
      path_groundwater_rates: "inputs/gwrates_dsm2.csv"
      path_dsm2_diversion_rates: "inputs/diversion_factors_dsm2.csv"
      path_dsm2_drainage_rates: "inputs/drainage_factors_dsm2.csv"

      # Output files
      path_dcd_output: "output/dcd_areas_dsm2.nc"
      path_dcd_node_output: "output/dcd_dsm2.nc"

      # DCD parameters
      start_water_year: 1922
      end_water_year: 2024
      leach_scale: 5.0
      runoff_rate: 0.75
      deep_percolation_rate: 0.25
      is_adding_waterbody_evaporation: True

The main input file starts with a `dcd` key. The `dcd` key contains the following. File paths in the input file are based on the current working directory where a user runs the model:

* path_subarea_info: subarea (or island) information file.
* path_irrigation_efficiency: irrigation efficiency file.
* path_leach_applied: applied leach water file.
* path_leach_drained: drained leach water file.
* path_groundwater_rates: groundwater rates file.
* path_detaw_output: DETAW output file as an input for DCD.
* path_dsm2_diversion_rates: DSM2 diversion rates file for DSM2 nodes.
* path_dsm2_drainage_rates: DSM2 drainage rates file for DSM2 nodes.
* start_water_year: start water year.
* end_water_year: end water year.
* leach_scale: leach scale.
* runoff_rate: runoff rate.
* deep_percolation_rate: deep percolation late.
* is_adding_waterbody_evaporation: a flag to add evaporation from open water bodies.
* path_dcd_output: DCD output file in the NetCDF format.
* path_dcd_node_output: DCD nodes output file  in the NetCDF format.

Input files
-----------

Below are sample input file formats. The input files can are in CSV or NetCDF format. Sample files can be found in the 'examples/inputs' directory.

- *path_subarea_info*: Input subarea information sample format

  CSV file containing subarea information.

..  code-block:: csv
    :caption: path_subarea_info

      area_id,name,uplow,acreage,docregion
      1,UNION ISLAND  (EAST),low,11853.0,Lower
      2,UNION ISLAND  (WEST),low,13711.0,Lower
      ...
      ...
      173,BBID in UNDESIGNATED AREA,up,13901.142857142857,Lower
      174,BBID in UNDESIGNATED AREA,up,34512.15126249121,Lower

- *path_irrigation_efficiency*: Input irrigation efficiency sample format

  CSV file containing irrigation efficiency information.

..  code-block:: csv
    :caption: path_irrigation_efficiency

      area_id,irrigation_efficiency
      1,0.85
      2,0.85
      ...
      ...
      173,0.85
      174,0.85

- *path_leach_applied*: Input leach applied sample format

  CSV file containing leach applied information.

..  code-block:: csv
    :caption: path_leach_applied

      area_id,10,11,12,1,2,3,4,5,6,7,8,9
      1,0.0,464.0,464.0,464.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
      2,0.0,529.0,529.0,529.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
      ...
      ...
      173,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
      174,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0

- *path_leach_drained*: Input leach drained sample format

  CSV file containing leach drained information.

..  code-block:: csv
    :caption: path_leach_drained

      area_id,10,11,12,1,2,3,4,5,6,7,8,9
      1,0.0,0.0,0.0,0.0,779.0,597.0,16.0,0.0,0.0,0.0,0.0,0.0
      2,0.0,0.0,0.0,0.0,889.0,681.0,18.0,0.0,0.0,0.0,0.0,0.0
      ...
      ...
      173,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
      174,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0

- *path_groundwater_rates*: Input groundwater rates sample format

  CSV file containing groundwater rates information.

..  code-block:: csv
    :caption: path_groundwater_rates

      year,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174
      1921,0.35,0.35,0.3,0.0,0.35,0.35,0.3,0.25,0.3,0.35,0.3,0.3,0.3,0.3,0.3,0.25,0.0,0.3,0.25,0.3,0.3,0.3,0.3,0.35,0.3,0.35,0.3,0.3,0.25,0.35,0.35,0.25,0.25,0.25,0.3,0.0,0.25,0.35,0.35,0.3,0.0,0.25,0.25,0.0,0.0,0.3,0.3,0.25,0.25,0.25,0.25,0.3,0.3,0.25,0.3,0.3,0.25,0.25,0.3,0.3,0.3,0.3,0.25,0.25,0.35,0.0,0.25,0.25,0.35,0.0,0.25,0.25,0.0,0.3,0.25,0.0,0.0,0.25,0.0,0.3,0.0,0.35,0.3,0.0,0.0,0.0,0.3,0.0,0.0,0.35,0.0,0.35,0.0,0.3,0.0,0.3,0.3,0.25,0.3,0.3,0.0,0.25,0.0,0.25,0.25,0.0,0.3,0.25,0.3,0.3,0.3,0.25,0.3,0.25,0.25,0.25,0.25,0.3,0.25,0.25,0.25,0.0,0.25,0.3,0.35,0.0,0.3,0.0,0.25,0.0,0.25,0.35,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.0,0.0,0.0,0.35,0.0,0.0,0.25,0.25,0.25,0.25,0.0,0.0,0.0,0.0,0.0,0.25,0.0,0.0,0.0,0.0,0.0,0.0,0.25,0.25,0.25,0.0,0.25,0.25,0.0,0.0,0.0,0.0
      1922,0.35,0.35,0.3,0.0,0.35,0.35,0.3,0.25,0.3,0.35,0.3,0.3,0.3,0.3,0.3,0.25,0.0,0.3,0.25,0.3,0.3,0.3,0.3,0.35,0.3,0.35,0.3,0.3,0.25,0.35,0.35,0.25,0.25,0.25,0.3,0.0,0.25,0.35,0.35,0.3,0.0,0.25,0.25,0.0,0.0,0.3,0.3,0.25,0.25,0.25,0.25,0.3,0.3,0.25,0.3,0.3,0.25,0.25,0.3,0.3,0.3,0.3,0.25,0.25,0.35,0.0,0.25,0.25,0.35,0.0,0.25,0.25,0.0,0.3,0.25,0.0,0.0,0.25,0.0,0.3,0.0,0.35,0.3,0.0,0.0,0.0,0.3,0.0,0.0,0.35,0.0,0.35,0.0,0.3,0.0,0.3,0.3,0.25,0.3,0.3,0.0,0.25,0.0,0.25,0.25,0.0,0.3,0.25,0.3,0.3,0.3,0.25,0.3,0.25,0.25,0.25,0.25,0.3,0.25,0.25,0.25,0.0,0.25,0.3,0.35,0.0,0.3,0.0,0.25,0.0,0.25,0.35,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.0,0.0,0.0,0.35,0.0,0.0,0.25,0.25,0.25,0.25,0.0,0.0,0.0,0.0,0.0,0.25,0.0,0.0,0.0,0.0,0.0,0.0,0.25,0.25,0.25,0.0,0.25,0.25,0.0,0.0,0.0,0.0
      ...
      ...
      2023,0.35,0.35,0.3,0.4,0.35,0.35,0.3,0.25,0.3,0.35,0.3,0.3,0.3,0.3,0.3,0.25,0.4,0.3,0.25,0.3,0.3,0.3,0.3,0.35,0.3,0.35,0.3,0.3,0.25,0.35,0.35,0.25,0.25,0.25,0.3,0.4,0.25,0.35,0.35,0.3,0.4,0.25,0.25,0.4,0.4,0.3,0.3,0.25,0.25,0.25,0.25,0.3,0.3,0.25,0.3,0.3,0.25,0.25,0.3,0.3,0.3,0.3,0.25,0.25,0.35,0.4,0.25,0.25,0.35,0.4,0.25,0.25,0.4,0.3,0.25,0.4,0.4,0.25,0.4,0.3,0.4,0.35,0.3,0.4,0.4,0.4,0.3,0.4,0.4,0.35,0.4,0.35,0.4,0.3,0.4,0.3,0.3,0.25,0.3,0.3,0.4,0.25,0.4,0.25,0.25,0.4,0.3,0.25,0.3,0.3,0.3,0.25,0.3,0.25,0.25,0.25,0.25,0.3,0.25,0.25,0.25,0.4,0.25,0.3,0.35,0.4,0.3,0.4,0.25,0.4,0.25,0.35,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.4,0.4,0.4,0.35,0.4,0.4,0.25,0.25,0.25,0.25,0.4,0.4,0.4,0.4,0.4,0.25,0.4,0.4,0.4,0.4,0.4,0.4,0.25,0.25,0.25,0.4,0.25,0.25,0.4,0.4,0.4,0.4
      2024,0.35,0.35,0.3,0.4,0.35,0.35,0.3,0.25,0.3,0.35,0.3,0.3,0.3,0.3,0.3,0.25,0.4,0.3,0.25,0.3,0.3,0.3,0.3,0.35,0.3,0.35,0.3,0.3,0.25,0.35,0.35,0.25,0.25,0.25,0.3,0.4,0.25,0.35,0.35,0.3,0.4,0.25,0.25,0.4,0.4,0.3,0.3,0.25,0.25,0.25,0.25,0.3,0.3,0.25,0.3,0.3,0.25,0.25,0.3,0.3,0.3,0.3,0.25,0.25,0.35,0.4,0.25,0.25,0.35,0.4,0.25,0.25,0.4,0.3,0.25,0.4,0.4,0.25,0.4,0.3,0.4,0.35,0.3,0.4,0.4,0.4,0.3,0.4,0.4,0.35,0.4,0.35,0.4,0.3,0.4,0.3,0.3,0.25,0.3,0.3,0.4,0.25,0.4,0.25,0.25,0.4,0.3,0.25,0.3,0.3,0.3,0.25,0.3,0.25,0.25,0.25,0.25,0.3,0.25,0.25,0.25,0.4,0.25,0.3,0.35,0.4,0.3,0.4,0.25,0.4,0.25,0.35,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.4,0.4,0.4,0.35,0.4,0.4,0.25,0.25,0.25,0.25,0.4,0.4,0.4,0.4,0.4,0.25,0.4,0.4,0.4,0.4,0.4,0.4,0.25,0.25,0.25,0.4,0.25,0.25,0.4,0.4,0.4,0.4



Output files
------------

Coming soon.