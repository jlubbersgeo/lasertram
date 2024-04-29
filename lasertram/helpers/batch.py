"""
batch module:
For batch processing operations in laserTRAM
"""


def process_spot(
    spot,
    raw_data,
    bkgd,
    keep,
    int_std,
    omit=None,
    despike=False,
    output_report=True,
):
    """a function to incorporate all the methods of the `LaserTRAM` class
    so a spot can be processed in an efficient and compact way.

    Args:
        spot (LaserTRAM spot object): an empty `LaserTRAM` spot object to be processed
        raw_data (pandas DataFrame): the raw counts per second dataframe to be assigned to the spot. Shape is (m x n) where m is the number of cycles through the mass range
        bkgd (tuple): (start, stop) pair of values corresponding to the analysis time where the background signal starts and stops
        keep (tuple): (start, stop) pair of values correpsonding to the analysis time where the interval signal for concentrations starts and stops
        internal_std (str): column name for the internal standard analyte (e.g., 29Si)
        omit (tuple): (start, stop) pair of values corresponding to the analysis time to be omitted from the `keep` interval. Defaults to None.
        despike (bool, optional): Whether or not to despike all analyte signals using the standard deviation filter from `LaserTRAM.despike_data()`. Defaults to False
        output_report (bool, optional): Whether or not to create a 1-row pandas DataFrame output report in the following format. Defaults to True.


    """
    # assign data to the spot
    spot.get_data(raw_data)
    # despike the data if desired
    if despike is True:
        spot.despike_data(analyte_list="all")
    # assign the internal standard analyte
    spot.assign_int_std(int_std)
    # assign intervals for background and ablation signal
    spot.assign_intervals(bkgd=bkgd, keep=keep, omit=omit)
    # assign and save the median background values
    spot.get_bkgd_data()
    # remove the median background values from the ablation interval
    spot.subtract_bkgd()
    # calculate detection limits based off background values
    spot.get_detection_limits()
    # normalize the ablation interval to the internal standard analyte,
    # get the median values, and the standard error
    spot.normalize_interval()

    if output_report is True:
        spot.make_output_report()
