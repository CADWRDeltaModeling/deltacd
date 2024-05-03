====================
Post-processing
====================

After the simulation is complete, the data can be post-processed using the `deltacd2dsm2` method. This method takes the following arguments:

.. code:: bash

    $    deltacd2dsm2 --input INPUT --output_dss OUTPUT

(`$` denotes the command prompt. Do not type it in.)

where `INPUT` is the path to the NetCDF file and `OUTPUT` is the path to the dss output file. Below is an example of how to use the `deltacd2dsm2` method:

.. code:: bash

    $   deltacd2dsm2 --input .\dcd_dsm2.nc --output_dss dcd_dsm2.dss