================
Installation
================

Release available on PyPI
--------------------------
Temporarily unavailable. Use the install from source code instructions below.
To install DeltaCD from PyPI, run the following command in a terminal:

.. code:: bash

    $ pip install deltacd

(`$` denotes the command prompt. Do not type it in.)

Using the source code
---------------------

If you want to install DeltaCD from source for development, you can do so by running the following commands in a terminal:

.. code:: bash

    $ git clone https://github.com/CADWRDeltaModeling/DeltaCD.git
    $ cd DeltaCD
    $ conda create -n deltacd python=3.11
    $ conda activate deltacd
    $ pip install versioneer
    $ pip install -e .
