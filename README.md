# lasertram

Welcome to the repository for `lasertram`, a package for the time resolved analysis for processing laser ablation inductively coupled plasma mass spectrometry (LA-ICP-MS) data.

## Structure

`lasertram` is comprised of two classes:

1. `LaserTRAM`: taking raw counts per second data and normalizing it to an internal standard as well as calculating all the required attributes needed for determining concentrations (i.e., detection limits, uncertainties,). To help with relatability, we can think of the object created by the `LaserTRAM` class as a spot object in which all of the data and metadata for that analysis is stored:

```python
    output_report = True
    spot = LaserTRAM(name = 'BCR-2G')

    # assign data to the spot
    spot.get_data(raw_data)

    # despike the data if desired
    if despike is True:
        spot.despike_data(analyte_list="all")

    # assign the internal standard analyte
    internal_std = "29Si"
    # instantiate our 'spot' object
    spot.assign_int_std(internal_std)

    # assign intervals for background and ablation signal
    bkgd = (5,10) #background interval in seconds
    keep = (20,40) #ablation interval to calculate concentrations from

    spot.assign_intervals(bkgd=bkgd, interval=keep)

    # assign and save the median background values
    spot.get_bkgd_data()

    # remove the median background values from the ablation interval
    spot.subtract_bkgd()

    # calculate detection limits in cps based off background values
    spot.get_detection_limits()

    # normalize the ablation interval to the internal standard analyte,
    # get the median values, and the standard error
    spot.normalize_interval()

    # create an output report for the spot analysis that contains
    # all the metadata and data from the above methods

    if output_report is True:
        spot.make_output_report()

    # To accomplish many of the above functions in one line of code
    # we have created the `process_data()` function that calls on
    # many of the methods within the `LaserTRAM` class. Ex:

    # create the spot object
    spot = LaserTRAM(name = 'BCR-2G')
    # assign data, establish intervals, normalize, and create output
    # all in one line
    process_data(
        spot,
        raw_data = df,
        bkgd = (5,10),
        keep = (20,50),
        internal_std = "29Si",
        despike = False,
        output_report = True
        )

```

Attributes follow naming conventions similar to the methods that create them. For example:

```python
spot.get_detection_limits()
```

creates the attribute `detection_limits` which is simply a 1D array of values corresponding to the detection limit of each analyte in counts per second.

2. `LaserCalc`: taking the output from the `LaserTRAM` class and calculating concentrations according to [Longerich et al., (1996)](https://pubs.rsc.org/en/content/articlepdf/1996/ja/ja9961100899?casa_token=KagVZMK9AgAAAAAA:pPybAcUcksXzD8UYmpjc_MI4uWa2tLELI_jC9Gtc1ycNTyH_tPK2meMuJ1SXVZICLwky-NggmJQXRA).

## Installation

## Contributing
