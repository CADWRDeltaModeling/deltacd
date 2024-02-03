=====
Usage
=====

DeltaCD consists of two parts: `DETAW` and `DCD`. `DETAW` calculates the amount of water demand from estimating evapotranspiration. `DCD` calculates the amount of water diversion based on the water demand from `DETAW`. Hence, `DETAW` is run first and `DCD` follows with the output from `DETAW` typically. This section describes how to use `DETAW` and `DCD`.

DETAW
-----

To run DETAW, navigate into a directory with input files, and launch the model in a command line interface as the following.

..  code-block:: bash

    detaw detaw_dsm2.yaml

