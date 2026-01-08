===============================================================
DETAW (Delta Evapotranspiration of Applied Water) Overview
===============================================================

Delta Evapotranspiration of Applied Water (DETAW) is designed to simulate the daily evapotranspiration(ET) and root zone water balance in subareas of the Sacramento-San Joaquin River Delta. The current version of DETAW in DeltaCD, is roughly the same as DETAW v2.1 from the previous DETAW-DCD. Hence, though the formats of many input files change slightly, `DETAW v2.1 user manual <https://github.com/CADWRDeltaModeling/DETAW-DCD/blob/master/DETAW/Documents/DETAW%20v2.1_user's%20manual.pdf>`_ is a good reference to use. The algorithms and the theories behind of DETAW are described in `DETAW v1.0 report <https://og-production-open-data-cnra-892364687672.s3.amazonaws.com/resources/6539f894-325d-4092-b34a-3139fd35c5b1/08detaw.pdf?Content-Type=application%2Fpdf&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAJJIENTAPKHZMIPXQ%2F20240206%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20240206T180739Z&X-Amz-Expires=3600&X-Amz-SignedHeaders=host&X-Amz-Signature=38bcbf3f4cf3de7b97e5b9c8ee0e1b5fb9443a4d911e4d7a041802191e977678>`_ and annual report chapters, `Implementing DETAW in Modeling Hydrodynamics and Water Quality in the Sacramento-San Joaquin Delta <https://og-production-open-data-cnra-892364687672.s3.amazonaws.com/resources/61bf1927-14a9-44ae-96ef-00fc8af9b725/method_2017_chapter3.pdf?Content-Type=application%2Fpdf&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAJJIENTAPKHZMIPXQ%2F20240206%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20240206T180930Z&X-Amz-Expires=3600&X-Amz-SignedHeaders=host&X-Amz-Signature=d2139c9152f4b8b7817ff7185513d142631deaa3efe7e3b21b52d6eb4e26f00e>`_ and `Estimates for Consumptive Water Demands in the Delta using DETAW <https://og-production-open-data-cnra-892364687672.s3.amazonaws.com/resources/91c46590-7383-41fd-8446-db60dae5e874/2006ch7.pdf?Content-Type=application%2Fpdf&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAJJIENTAPKHZMIPXQ%2F20240206%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20240206T180924Z&X-Amz-Expires=3600&X-Amz-SignedHeaders=host&X-Amz-Signature=3bb64bc65a2b348dd2b442b60c6f03730154ae7b36f258e9dc2bc84cc3f1f39a>`_.

Main input file
===============

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

The main input file starts with a `detaw` key. The `detaw` key contains the following item-value pairs. The item and values they represent are described below. File paths in the input file are based on the current working directory where a user runs the model:

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

Input files
-----------
Below are sample input file formats.

- *input_pcp*: Input precipitation sample format

  CSV file containing station precipitation in mm.

..  code-block:: csv
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

..  code-block:: csv
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

..  code-block:: csv
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

..  code-block:: txt
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

..  code-block:: txt
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

..  code-block:: txt
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
