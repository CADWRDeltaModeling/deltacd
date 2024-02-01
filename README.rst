===============================
DeltaCD
===============================

DeltaCD is a numerical model to estimate the Delta channel depletion. DeltaCD is a successor of DETAW-DCD, re-written in pure Python. More information about DeltaCD can be found at `Chapter 3 of one of our Annual Reports <https://og-production-open-data-cnra-892364687672.s3.amazonaws.com/resources/dcabdb20-e638-4cf5-b199-78e78f0d482f/2023-bay-delta-annual-report.pdf?Content-Type=application%2Fpdf&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAJJIENTAPKHZMIPXQ%2F20240131%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20240131T230242Z&X-Amz-Expires=3600&X-Amz-SignedHeaders=host&X-Amz-Signature=e806d752d23a4c74d87efedb159ded16cb7ceb3bede78628ed759986fd4179cd>`_.

An online documentation is available at https://cadwrdeltamodeling.github.io/deltacd/.

Installation
===============================

Developer Installation
----------------------

For developer to use the source code directly, follow this section.

First, clone the DeltaCD repository from GitHub, navigate into it.

Create an environment with a suitable name, for example, *deltacd*:

``conda create -n deltacd python=3.11``

Once a new environment is created, activate it:

``conda activate deltacd``

Install DeltaCD with ``-e``, editable option:

``pip install -e .``

Usage
===============================

DeltaCD package has two models, DETAW and DCD. They can be launched from the command line (after a Python environment with DeltaCD is activated).

DETAW needs to run first, for example:

``detaw detaw_dsm2.yaml``

Followed by DCD:

``dcd dcd_dsm2.yaml``

Example yaml and input files are provided in the *examples* directory in DeltaCD package.

While not required, it would be a good practice to make a copy of the *examples* directory to a project or working directory and rename it suitably.


Getting Help
===============================

Any questions and issues? We would like to help and improve the package.

For any questions, please contact DeltaCD developers. Please report bugs and issues at `GitHub issues <https://github.com/CADWRDeltaModeling/deltacd/issues>`_.
