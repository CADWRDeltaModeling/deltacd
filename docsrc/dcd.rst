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

Output files
------------

Coming soon.