=========
Changelog
=========

0.8.1 - 2023-02-27
=============================
A bug fix version of v0.8

NOTE: There are known differences between results from the DETAW-DCD and DeltaCD. The codes and the results are being reviewed though most of the differences are from improvements.

Added
-----
*  The input precip and temperature files are updated till 01/31/2023, and the future is padded with 2020 water year climate.

Fixed
-----
* The incorrect lowlands/uplands values in the codes are fixed.

0.8.0 - 2023-01-17
=============================
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
