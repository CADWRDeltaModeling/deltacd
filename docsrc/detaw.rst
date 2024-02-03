=====
DETAW
=====

Main input file
---------------

The main input file of DETAW follows YAML format. YAML is a human-readable data language format, common for configuration files (`Wikipedia YAML page <https://en.wikipedia.org/wiki/YAML>`_). The following is an example of the DETAW main input file. Currently, the parameters in the main input file keep the names from the previous DETAW-DCD. Note that some of the options are not well maintained, and they are subject to change.

..  code-block:: yaml
    :caption: detaw_dsm2.yaml

    detaw:
      # Model inputs
      start_water_year: 1922
      end_water_year: 2024
      input_pcp: 'inputs/mm_pcp.csv'
      input_temperature: 'inputs/LODI_PT.csv'
      landuse : 'inputs/landuse_dsm2.csv'
      et_correction: 'inputs/Percentage.csv'
      critical: 'inputs/critical.csv'
      noncritical: 'inputs/noncritical.csv'
      # Model output control
      detaw_output: 'output/detawoutput_dsm2.nc'
      precip_output: 'output/precip_dsm2.nc'
      et_output: 'output/ET0.nc'
      daily_output: 1
      monthly_output: 1
      yearly_output: 0
      delta_output: 0
      daily_output_unit: 1
      for_dsm2_only: 1

The main input file starts with a `detaw` key. The `detaw` key contains the following:

* `start_water_year`: The model start year in water year. The water year starts from October 1st and ends on September 30th of the following year.
* `end_water_year`: The model end year in water year.
* `input_pcp`: The precipitation CSV file.
* `input_temperature`: The temperature CSV file.
* `landuse`: The land use CSV file.
* `et_correction`: The ET correction CSV file.
* `critical`: The crop coefficient CSV file for critical (dry, below average) water years.
* `noncritical`: The crop coefficient CSV file for non-critical water years.
* `detaw_output`: The DETAW netCDF output filename.
* `precip_output``: The precipitation netCDF output filename.
* `et_output`: The ET netCDF output filename.
* `daily_output`: The flag to output daily data. If 1, the model outputs daily data. If 0, the model does not output daily data.
* `monthly_output`: The flag to output monthly data. If 1, the model outputs monthly data. If 0, the model does not output monthly data.
* `yearly_output`: The flag to output yearly data. If 1, the model outputs yearly data. If 0, the model does not output yearly data.
* `delta_output`: The flag to output delta data. If 1, the model outputs total of the values. If 0, the model does not output it.
* `daily_output_unit`: The unit of daily output. If 1, the unit is in ft. If 0, the unit is in mm. Not actively maintained.
* `for_dsm2_only`: The flag to output data for DSM2. If 1, the model outputs data only relevant to DSM2. If 0, the outputs are written as they are.
