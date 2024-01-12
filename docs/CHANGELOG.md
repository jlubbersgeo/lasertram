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
