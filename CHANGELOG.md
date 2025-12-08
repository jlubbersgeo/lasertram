# Release Notes

## Upgrading

You can always update to the latest version of `lasertram` using pip:

```
pip install --upgrade lasertram
```

or install a specific version:

```
pip install lasertram==version-number
```

To check which version of `lasertram` you have installed:

```
pip show lasertram
```

## Maintenance team

Current and past maintainers of `lasertram`:

- [@jlubbersgeo](https://github.com/jlubbersgeo)

## 1.0.3 (DATE HERE)

This release mostly revolves around improvements to the parts of the library that deal with type checking, formatting of input data, and other things that better help the user avoid common pitfalls related to improper data formatting. There has been some reorganizing as well:

- moving the test data to its own directory: `test_data`. This contains subdirectories `basic_SRM_example` and `computers_and_geosciences_examples` which contain data used for the docs and associated manuscript [https://doi.org/10.1016/j.acags.2025.100225](https://doi.org/10.1016/j.acags.2025.100225) , respectively.
- 

### :sparkles: Features
- added module for formatting checks:
  - added in functionality to check that input data to `LaserTRAM` is properly formatted. This includes checking for:
  - added functionality to check that the chosen calibration standard has published values for all analytes in the input dataset. If not throw error.
  - added functionality to check that the SRM dataset is uploaded in the correct format
  - added functionality to check that data loaded into LaserCalc are properly formatted.
- added ability to calculate concentrations using internal standard values that are in either ppm element, wt% element, or wt% oxide.
- added support for more internal standard analytes. While any element can be used as an internal standard in ppm units, supported oxides can be found by checking the `lasertram.helpers.conversions.supported_internal_standard_oxides` list. `mendeleev` was removed as a dependency because of this.
- `LaserCalc` now checks for and fixes duplicate spot names. Because a lot of the code uses pandas indexing and the spot names are the index, duplicates will break the `LaserCalc` process. `LaserCalc` now will append `_A`, `_B`, etc. to duplicate spot names to make them unique to prevent this break. 
- `make_LT_ready_folder()` function now deals with empty `.csv` files better by skipping them and throwing a message to the user rather than breaking. It also now displays more verbose messages to the user about what files are being processed so if something does break they have better context. Because of this `tqdm` and `tabulate` were added as dependencies.


### :test_tube: Tests

- `pytest` coverage has been updated to better reflect the way LaserTRAM and LaserCalc output data. This also increases the coverage to > 90%.
- switch to using computers and geosciences dataset for testing rather than basic SRM dataset. Better reflects actual use cases for users rather than "idealized" data.
- added data to be used for testing how the `preprocessing` module deals with raw data from ThermoFisher quadrupole ICP-MS. This is in the `tests/raw` directory.

### :bug: Bug Fixes

- fixed bad path handling with functions related to loading test data
- migrated how `b.d.l.` is handled internally to avoid deprecation warnings from `pandas`
- fixed bug with `make_LT_ready_folder()` where empty `csv` files would cause the function to break.

## 1.0.1 and 1.0.2

_1.0.2 is the same as 1.0.1 but I screwed up the PyPI distribution so I'm redoing it as 1.0.2_

### :sparkles: First release after manuscript acceptance :sparkles:

- [_`lasertram`: A Python library for the time resolved analysis of laser ablation inductively coupled plasma mass spectrometry data_](https://doi.org/10.1016/j.acags.2025.100225)

### :sparkles: Features

- added a module for basic plotting of time-series and uncertainties in `LaserTRAM`
- added a module for preprocessing which has functions related to loading in test data ready for `LaserTRAM` and for taking raw .csv files from Thermo or Agilent quadrupoles and converting them directly to `pandas.DataFrames` ready for `LaserTRAM`. 

### :books: Documentation

- added more example jupyter notebooks on signal de-spiking and region omission
- Added much better docstring coverage to _everything_
- Added much better line-by-line coverage to _everything_

### :bug: Bug Fixes

- fixed pytest bugs with move to new version. 

## 1.0.0 (07/11/2024)

### :sparkles: USGS Software Release approval :sparkles:

- Addresses reviewer comments from Issue !1 related to USGS internal code review

### :sparkles: Features

- updated syntax to reflect `pylint` recommendations (e.g., snake case variables, correct negation/equivalence referencing)

### :books: Documentation

- transitioned to solely using a `pyproject.toml` rather than combination of `setup.cfg` + `pyproject.toml`

## v0.1.2 (5/1/2024)

### :bug: Bug Fixes

- fixed indexing issue in `despike_data`

### :books: Documentation

- added more coverage in tests to now include `get_secondary_standard_accuracies` and `despike_data`

## v0.1.1 (4/29/2024)

### :books: Documentation

- added `lasertram.__version__`

### :tada:

- changed naming of all functions that refer to standards for consistency
  - Any attribute or method that references an internal standard uses `int_std`
  - Any attribute or method that references the calibration standard uses `calibration_std`

## v0.1.0 (4/25/2024)

### :art: :tada: :mega: Code structure change

A re-organization of the package structure was completed to set up for further growth. This has created the following modules for `lasertram`: 

- `tram`: holds the class 'LaserTRAM`
- `calc`: holds the class `LaserCalc`
- `helpers`: holds the submodules `batch` and `conversions`, for batch processing and unit conversions, respectively. 
  
An example folder structure for this is as follows: 


```
└── 📁lasertram
    └── 📁calc
        └── calc.py
        └── __init__.py
    └── 📁helpers
        └── batch.py
        └── conversions.py
        └── __init__.py
    └── 📁tram
        └── tram.py
        └── __init__.py
    └── __init__.py
```
The only real effect on the user is the change in import statements to now be the following:
```python
from lasertram import LaserTRAM, LaserCalc, batch
```
or to import an entire module:
```python
from from lasertram.helpers import batch, conversions

```

## v0.0.13 (4/08/2024)

### :bug: Bug Fixes

- Fixed undeclared list instantiation in `LaserTRAM.despike()`.

## v0.0.12 (1/12/2024)

### :bug: Bug Fixes

- Fixed issue with `calculate_uncertainties()` internal error calculation dropping the wrong terms

## v0.0.11 (1/12/2024)

### :bug: Bug Fixes

- Fixed issue in `drift_check()` where data do not have a timestamp and the index of analysis was being treated as a pandas Series. It is now a numpy array for consistency and better compatibility with the statsmodels package.
- Updated analyte indexing for regex deprecation warning related to string literals.

### :books: Documentation

- Updated Basic Usage page to reflect new uncertainty output implemented in `v0.0.10`

## v0.0.10 (1/11/2024)

### :sparkles: Features

- added support for both internal and external precision on output.
  - Now on output uncertainties are displayed with column header suffixes as either `_exterr` or `interr` to reflect overall uncertainties that incorporate the uncertainty in the calibration standard accepted value or those that don't, respectively. The vast majority of use cases should utilize `_exterr` values unless one is comparing datasets processed with the same calibration standard.

## v0.0.9 (10/26/2023)

### :bug: Bug Fixes

- fixed type issue with internal standard compositions being defaulted to `int`
- fixed pandas indexing in `drift_check()` to support future move to `.iloc`

### :sparkles: Features

- added support for both GUI (LaserTRAM-DB) and API (`lasertram.LaserTRAM`) outputs to be used as inputs into `lasertram.LaserCalc()`
- changed all naming instances of `calibration_standard` to `calibration_std`. Previously some there were a mix of attributes that had one or the other. Now anything that deals with the calibration standard (e.g., `self.calibration_std_stats`) will always have `calibration_std` in the name for consistency.

## v0.0.8 (10/25/2023)

### :bug: Bug Fixes

- Fixed bug that overwrote and ignored omission interval in `LaserTRAM.normalize_interval()`
- Future proofed output report for deprecated pandas indexing. Now properly uses `.iloc`
- Fixed interval `dtype` inconsistencies so they are all type `float64` rather than `object`

### :books: Documentation

- Added general motivation to documentation landing page

### :rocket: New Features

- Added required dependencies to `setup.cfg`

## v0.0.7 (10/24/2023)

### :bug: Bug Fixes

- Fixed `omitted_region` being counted as experiment analyte

### :books: Documentation

- fixed filepath issues and typos in Basic Useage tutorial
- moved all test data to `test_data` folder
- created `tests` folder for future tests to be held
- created `CHANGELOG.md`

## :tada: v0.0.5 (10/20/2023)

- First release of `lasertram` for public use.
