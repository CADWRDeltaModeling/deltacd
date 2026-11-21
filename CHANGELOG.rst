=========
Changelog
=========

v0.10.0 - Unreleased
====================
Modular rewrite of ``deltacd.detaw`` using ``xarray``, ``pandas``, and compact ``@numba.njit`` kernels.

Added
-----
* Modular ``deltacd.detaw`` subpackage (``config``, ``params``, ``io``, ``meteorology``, ``et``, ``water_balance``, ``etaw``, ``model``).
* Pluggable :math:`ET_0` equations (Hargreaves–Samani and ASCE Standardized Penman–Monteith) and precomputed/gridded ``DatasetMeteorologyProvider``.
* ``ExternalETProvider`` supporting external / **OpenET** daily or monthly :math:`ET_c` and :math:`ET_{aw}` NetCDF and CSV datasets.
* Declarative ``CropBehavior`` system supporting custom land-use categories without hardcoded crop indices.
* Educational Jupyter notebooks (``notebooks/01_detaw_theory_and_walkthrough.ipynb`` and ``notebooks/02_detaw_modular_et_and_openet.ipynb``) and AI developer guides (``AGENTS.md``, ``INSTRUCTIONS.md``, ``PLANS.md``, ``SKILLS.md``).

Fixed
-----
* Fixed positional ``iloc[:, 1]`` bug in ``read_landuse`` that read ``DATE`` instead of ``TYPE`` in Delta land-use CSVs, restoring critical/dry water-year crop parameters.
* Fixed leap-year February 29 day-count, ``IKc`` curve length, and ``calc_cetc`` index-wrapping bugs.
* Fixed 1st-of-month ``Mon`` off-by-one bug in daily channel seepage calculations.
* Fixed +1 day timestamp shift in ``et_output`` and ``precip_output`` NetCDF files caused by the prior-September-30 spin-up row.
* Fixed unassigned ``numpy.clip`` no-op on station-weighted precipitation and removed C++ float-truncation drift.
* Eliminated fixed column-order, row-order, and crop-order assumptions across all DETAW readers and DCD crop selection.

v0.9.0 - 2024-XX-XX
===================
This release is a rather major release that breaks input file compatibility.

Changed
-------
* DeltaCD splits subareas in subarea depletion calculation codes, similarly to the original DETAW-DCD, not before the evapotranspiration calculation (DETAW) in the input files. This makes the input file preparation easy and similar to what the original DETAW-DCD required. This change, however, breaks the input file compatibility with the previous versions of DeltaCD. A few new input files are introduced to accommodate this change.
* When file names in the main YAML file are provided relatively, the paths are relative to the main YAML file, not the current working directory anymore.

Added
-----
* A simple example input is added.

Removed
-------
* The full model inputs are moved to a separate repository, 'deltacd_inputs'.

Fixed
-----
* A bug in codes to detect leap years in DETAW is fixed.

v0.8.1 - 2023-02-27
===================
A bug fix version of v0.8

NOTE: There are known differences between results from the DETAW-DCD and DeltaCD. The codes and the results are being reviewed though most of the differences are from improvements.

Added
-----
* The input precipitation and temperature files are updated till 01/31/2023, and the future is padded with 2020 water year climate.

Fixed
-----
* The incorrect lowlands/uplands values in the codes are fixed.

v0.8.0 - 2023-01-17
===================
This release is the first pre-release of DeltaCD, renamed and reimplemented DETAW-DCD model. We now use a new name, DeltaCD and a new versioning to make it user-friendly by using more user inputs and by following the Pythonic ways. We adopted common file formats such as CSV, YAML, and NetCDF, and the DSS files are not used anymore. DeltaCD aims to clean and vectorized implementation so that we can accommodate future requirements easily while improving the run speed. This release improves the run speed about ten times.

NOTE: There are known differences between results from the DETAW-DCD and DeltaCD. The codes and the results are being reviewed though most of the differences are from improvements.

Added
------
* The simulation period can be chosen in the user inputs, which is now in YAML.

Changed
-------
* The model is packaged in a proper Python package.
* The model control parameters are in YAML formats.
* The data files format for the previously DCD are all in comma-separated values (CSV).
* The outputs are in the NetCDF format.
* The data splitting for certain irrigation districts and Calsim

Removed
-------
* The option to choose a target model is dropped. Now the users need to provide model-corresponding input data sets.
* The support for and the use of DSS files are removed.
* The Fortran module called externally in the older DCD is all reimplemented in Python, and it is removed.

Others
------
* The adjustment for 1977 is not reimplemented yet.
* A tool to convert NetCDF outputs to DSS files is not implemented yet.
