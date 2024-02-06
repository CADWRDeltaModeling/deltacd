=================================
DETAW (Delta Evapotranspiration of Applied Water)
=================================

Delta Evapotranspiration of Applied Water (DETAW) is designed to simulate the daily evapotranspiration(ET) and root zone water balance in subareas of the Sacramento-San Joaquin River Delta. The current version of DETAW in DeltaCD, is roughly the same as DETAW v2.1 from the previous DETAW-DCD. Hence, though the formats of many input files change slightly, `DETAW v2.1 user manual <https://github.com/CADWRDeltaModeling/DETAW-DCD/blob/master/DETAW/Documents/DETAW%20v2.1_user's%20manual.pdf>`_ is a good reference to use. The algorithms and the theories behind of DETAW are described in `DETAW v1.0 report <https://og-production-open-data-cnra-892364687672.s3.amazonaws.com/resources/6539f894-325d-4092-b34a-3139fd35c5b1/08detaw.pdf?Content-Type=application%2Fpdf&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAJJIENTAPKHZMIPXQ%2F20240206%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20240206T180739Z&X-Amz-Expires=3600&X-Amz-SignedHeaders=host&X-Amz-Signature=38bcbf3f4cf3de7b97e5b9c8ee0e1b5fb9443a4d911e4d7a041802191e977678>`_ and annual report chapters, `Implementing DETAW in Modeling Hydrodynamics and Water Quality in the Sacramento-San Joaquin Delta <https://og-production-open-data-cnra-892364687672.s3.amazonaws.com/resources/61bf1927-14a9-44ae-96ef-00fc8af9b725/method_2017_chapter3.pdf?Content-Type=application%2Fpdf&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAJJIENTAPKHZMIPXQ%2F20240206%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20240206T180930Z&X-Amz-Expires=3600&X-Amz-SignedHeaders=host&X-Amz-Signature=d2139c9152f4b8b7817ff7185513d142631deaa3efe7e3b21b52d6eb4e26f00e>`_ and `Estimates for Consumptive Water Demands in the Delta using DETAW <https://og-production-open-data-cnra-892364687672.s3.amazonaws.com/resources/91c46590-7383-41fd-8446-db60dae5e874/2006ch7.pdf?Content-Type=application%2Fpdf&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAJJIENTAPKHZMIPXQ%2F20240206%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20240206T180924Z&X-Amz-Expires=3600&X-Amz-SignedHeaders=host&X-Amz-Signature=3bb64bc65a2b348dd2b442b60c6f03730154ae7b36f258e9dc2bc84cc3f1f39a>`_.

Input files
-----------

Main input file
^^^^^^^^^^^^^^^

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

The main input file starts with a `detaw` key. The `detaw` key contains the following. File paths in the input file are based on the current working directory where a user runs the model:

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

Output files
------------

Coming soon.
