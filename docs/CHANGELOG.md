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
