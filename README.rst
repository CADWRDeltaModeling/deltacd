===============================
DeltaCD
===============================

Code repo for DETAW and Delta Channel Depletions Model.

Installation
===============================

Developer Installation
----------------------

First clone the deltacd package.

It is a good practice to create a new environment to setup the package. Create an environment with a suitable *my_env_name* e.g. deltacd.

``conda create -n my_env_name python=3.8``

Once the *my_env_name* environment is created, activate it using

``conda activate my_env_name``

Installing the package can be accomplished by

``pip install -e .``

These steps create an environment and install the deltacd package which can then be accessed from the command line.

Usage
===============================

The deltacd package has two models, detaw and dcd. They can be accessed from the command line once the environment is activated.

The detaw model is run first

``detaw .\detaw_dsm2.yaml``

Followed by dcd

``dcd .\deltacd_dsm2.yaml``

Example yaml and input files are provided in the *examples* directory in the deltacd package. While not required, it would be a good practice to make a copy of the *examples* directory to a project or working directory and rename it suitably.
