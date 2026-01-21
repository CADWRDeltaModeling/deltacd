=====
Usage
=====

DeltaCD consists of two parts: `DETAW` and `DCD`. `DETAW` calculates the amount of water demand from estimating evapotranspiration. `DCD` calculates the amount of water diversion based on the water demand from `DETAW`. Hence, `DETAW` is run first and `DCD` follows with the output from `DETAW` typically. This section describes how to use `DETAW` and `DCD`.

DETAW
-----

To run DETAW, navigate into a directory with input yaml files, and launch the model in a command line interface similar to the following with your input file. This assumes that DeltaCD is installed to your Python environment and readily available in your command line interface.

..  code-block:: bash

    detaw detaw_dsm2.yaml

DCD
---

To run DCD, navigate into a directory with input yaml files, and launch the model in a command line interface similar to the following with your input file. This assumes that DeltaCD is installed to your Python environment and readily available in your command line interface.

..  code-block:: bash

    dcd dcd_dsm2.yaml
