"""

TRAM module: (T)ime (R)esolved (A)nalysis (M)odule.

For taking raw counts per second data from a Laser Ablation Inductively
Coupled Plasma Mass Spectrometry (LA-ICP-MS) experiment, choosing an
interval to be turned into a concentration, normalizing
that interval to an internal standard and outputting that value + other
metadata


"""

import numpy as np
import pandas as pd


def _z_filter(signal, std_devs, window=50):
    """compute a z filter for a given signal

    Parameters
    ----------
    signal : pandas.Series
        the timeseries signal to be filtered
    std_devs : float
        number of standard deviations away from the mean to
        be considered an outlier
    window : int, optional
        number of points in the moving average window, by default 50
    """

    # rolling calculation setup
    roll = signal.rolling(window=window, min_periods=1, center=True)
    # rolling average
    avg = roll.mean()
    # rolling standard deviation
    std = roll.std(ddof=0)
    # z score
    z = (signal - avg) / std
    # boolean if a given z score is within the bounds set up by
    # the standard deviation limits
    m = z.between(-std_devs, std_devs)

    # original signal where z score is between bounds
    # mean of rolling window otherwise
    despiked = signal.where(m, avg)

    return despiked


class LaserTRAM:
    """
    # LaserTRAM
    The class `LaserTRAM` which is devoted to the "time resolved analysis"
    operations during the laser data reduction process. To be used in
    conjunction with the `LaserCalc` class. The general idea is that
    this creates an object that contains all the information related
    to one individual spot analysis.

    """

    def __init__(self, name):
        """

        Args:
            name (str): your sample name i.e. the value in the `SampleLabel` column of the LT_ready file
        """
        # all attributes in relative chronological order that they are created in
        # if everything is done correctly. These all will get rewritten throughout the
        # data processing pipeline but this allows us to see what all the potential attributes
        # are going to be from the beginning (PEP convention)

        # for the math involved please see:

        # name of the lasertram spot object
        self.name = name

        # boolean flag for whether or not the data have been
        # despiked
        self.despiked = False

        # list of elements that have been despiked. Also may be 'all'
        self.despiked_elements = None

        # data from a single spot to be processed. 2D pandas dataframe
        self.data = None

        # self.data but as a 2D numpy matrix. Equivalent to self.data.values
        self.data_matrix = None

        # list of analytes in the analysis
        self.analytes = None

        # datetime corresponding to the analysis
        self.timestamp = None

        # string representation internal standard analyte for the processing.
        # this is just the column header of the analyte chosen as the internal
        # standard e.g., "29Si"
        self.int_std = None

        # column number in self.data_matrix that denotes the internal standard analyte
        # data. Remember python starts counting at 0!
        self.int_std_loc = None

        # background interval start time
        self.bkgd_start = None

        # background interval stop time
        self.bkgd_stop = None

        # desired ablation interval start time
        self.int_start = None

        # desired ablation interval stop time
        self.int_stop = None

        # row in self.data corresponding to self.bkgd_start
        self.bkgd_start_idx = None

        # row in self.data corresponding to self.bkgd_stop
        self.bkgd_stop_idx = None

        # row in self.data corresponding to self.int_start
        self.int_start_idx = None

        # row in self.data corresponding to self.int_stop
        self.int_stop_idx = None

        # desired omitted region start time
        self.omit_start = None

        # desired omitted region stop time
        self.omit_stop = None

        # row in self.data corresponding to self.omit_start
        self.omit_start_idx = None

        # row in self.data corresponding to self.omit_stop
        self.omit_stop_idx = None

        #
        self.omitted_region = None

        # 1D array of median background values [self.bkgd_start - self.bkgd_stop)
        # that is len(analytes) in shape
        self.bkgd_data_median = None

        # 1D array of detection limits in counts per second
        # that is len(analytes) in shape
        self.detection_limits = None

        # 2D array of background corrected data over the self.int_start - self.int_stop
        # region
        self.bkgd_subtract_data = None

        # 2D array of background corrected data over the self.int_start - self.int_stop
        # region that is normalized to the internal standard
        self.bkgd_subtract_normal_data = None

        # 1D array of median background corrected normalized values over the self.int_start - self.int_stop
        # retion that is len(analytes) in shape
        self.bkgd_subtract_med = None

        # 1D array of 1 standard error of the mean values for each analyte over the
        # self.int_start - self.int_stop region
        self.bkgd_subtract_std_err = None

        #
        self.bkgd_subtract_std_err_rel = None

        # 1D pandas dataframe that contains many of the attributes created during the
        # LaserTRAM process:
        # |timestamp|Spot|despiked|omitted_region|bkgd_start|bkgd_stop|int_start|int_stop|norm|norm_cps|analyte vals and uncertainties -->|
        # |---------|----|--------|--------------|----------|---------|---------|--------|----|--------|----------------------------------|
        self.output_report = None

    def get_data(self, df, time_units="ms"):
        """assigns raw counts/sec data to the object

        Args:
            df (pandas DataFrame): raw data corresponding to the spot being processed i.e., `all_data.loc[spot,:]` if `all_data` is the LT_ready file
            time_units (str): string denoting the units for the `Time` column. Used to convert input time values to seconds. Defaults to 'ms'.
        """
        # get data and set index to "SampleLabel" column
        self.data = df.reset_index()
        self.data = self.data.set_index("SampleLabel")

        # convert time units from ms --> s if applicable
        if time_units == "ms":
            self.data["Time"] = self.data["Time"] / 1000
        elif time_units == "s":
            pass

        # just numpy matrix for data
        self.data_matrix = self.data.iloc[:, 1:].to_numpy()

        # list of analytes in experiment
        self.analytes = self.data.loc[:, "Time":].columns.tolist()[1:]

        # need to add check for if this exists otherwise there is no timestamp attribute
        self.timestamp = str(self.data.loc[:, "timestamp"].unique()[0])

    def assign_int_std(self, int_std):
        """assigns the spot an internal standard
        analyte

        Args:
            int_std (str): the name of the column for the internal standard analyte e.g., "29Si"
        """

        # set the internal standard analyte
        self.int_std = int_std

        # get the internal standard array index
        self.int_std_loc = np.where(np.array(self.analytes) == self.int_std)[0][0]

    def assign_intervals(self, bkgd, keep, omit=None):
        """assigns the intervals to be used as background
        as well as the portion of the ablation interval to
        be used in calculating concentrations

        Args:
            bkgd (tuple): (start, stop) pair of values corresponding to the analysis time where the background signal starts and stops
            keep (tuple): (start, stop) pair of values correpsonding to the analysis time where the interval signal for concentrations starts and stops
            omit (tuple): (start, stop) pair of values corresponding to the analysis time to be omitted from the `keep` interval. Defaults to None.
        """

        # set background and interval times in s
        self.bkgd_start = bkgd[0]
        self.bkgd_stop = bkgd[1]
        self.int_start = keep[0]
        self.int_stop = keep[1]

        # equivalent background and interval times but as indices
        # in their respective arrays
        self.bkgd_start_idx = np.where(self.data["Time"] > self.bkgd_start)[0][0]
        self.bkgd_stop_idx = np.where(self.data["Time"] > self.bkgd_stop)[0][0]
        self.int_start_idx = np.where(self.data["Time"] > self.int_start)[0][0]
        self.int_stop_idx = np.where((self.data["Time"] > self.int_stop))[0][0]

        # boolean whether or not there is an omitted region
        self.omitted_region = False
        # if omission is true, set those start and stop times like above
        if omit:
            self.omit_start = omit[0]
            self.omit_stop = omit[1]
            self.omit_start_idx = (
                np.where(self.data["Time"] > self.omit_start)[0][0] - self.int_start_idx
            )
            self.omit_stop_idx = (
                np.where(self.data["Time"] > self.omit_stop)[0][0] - self.int_start_idx
            )

            self.omitted_region = True

    def get_bkgd_data(self):
        """
        uses the intervals assigned in `assign_intervals` to take the median
        value of all analytes within that range and use them as the
        background signal that gets subtracted from the ablation signal
        """
        # median background data values
        self.bkgd_data_median = np.median(
            self.data_matrix[self.bkgd_start_idx : self.bkgd_stop_idx, 1:], axis=0
        )

    def get_detection_limits(self):
        """
        Calculates detection limits in counts per second for each analyte. This
        is defined as the value that is three standard deviations away from the
        background.
        """

        self.detection_limits = np.std(
            self.data_matrix[self.bkgd_start_idx : self.bkgd_stop_idx, 1:], axis=0
        ) * 3 + np.median(
            self.data_matrix[self.bkgd_start_idx : self.bkgd_stop_idx, 1:], axis=0
        )

    def subtract_bkgd(self):
        """
        subtract the median background values calculated in `get_bkgd_data`
        from the signal in the "keep" interval established in `assign_intervals`

        """
        self.bkgd_subtract_data = (
            self.data_matrix[self.int_start_idx : self.int_stop_idx, 1:]
            - self.bkgd_data_median
        )

    def normalize_interval(self):
        """
        normalize the analytes from the "keep" portion of the signal
        the internal standard analyte. This is done by simply
        dividing the analytes by the internal standard analyte.

        This also calculates the median normalized value, its
        standard error of the mean, and relative standard error
        of the mean.
        """

        # set the detection limit thresholds to be checked against
        # with the interval data. This basically takes the detection limits
        threshold = self.detection_limits - self.bkgd_data_median

        # if there's an omitted region, remove it from the data to be further processed
        # for the chosen interval
        if self.omitted_region is True:
            self.bkgd_subtract_normal_data = np.delete(
                self.bkgd_subtract_data,
                np.arange(self.omit_start_idx, self.omit_stop_idx),
                axis=0,
            ) / np.delete(
                self.bkgd_subtract_data[:, self.int_std_loc][:, None],
                np.arange(self.omit_start_idx, self.omit_stop_idx),
                axis=0,
            )

        else:
            self.bkgd_subtract_normal_data = (
                self.bkgd_subtract_data
                / self.bkgd_subtract_data[:, self.int_std_loc][:, None]
            )

        # get background corrected and normalized median values for an interval
        self.bkgd_subtract_med = np.median(self.bkgd_subtract_normal_data, axis=0)
        self.bkgd_subtract_med[
            np.median(self.bkgd_subtract_data, axis=0) <= threshold
        ] = -9999
        self.bkgd_subtract_med[np.median(self.bkgd_subtract_data, axis=0) == 0] = -9999

        # standard error of the mean for the interval region
        self.bkgd_subtract_std_err = self.bkgd_subtract_normal_data.std(
            axis=0
        ) / np.sqrt(abs(self.int_stop_idx - self.int_start_idx))

        self.bkgd_subtract_std_err_rel = 100 * (
            self.bkgd_subtract_std_err / self.bkgd_subtract_med
        )

    def make_output_report(self):
        """
        create an output report for the spot processing. This is a
        pandas DataFrame that has the following format:

        |timestamp|Spot|despiked|omitted_region|bkgd_start|bkgd_stop|int_start|int_stop|norm|norm_cps|analyte vals and uncertainties -->|
        |---------|----|--------|--------------|----------|---------|---------|--------|----|--------|----------------------------------|
        """
        if self.despiked is True:
            despike_col = self.despiked_elements
        else:
            despike_col = "None"

        if self.omitted_region is True:
            omitted_col = (
                self.data["Time"].iloc[self.omit_start_idx + self.int_start_idx],
                self.data["Time"].iloc[self.omit_stop_idx + self.int_start_idx],
            )
        else:
            omitted_col = "None"

        spot_data = pd.DataFrame(
            [
                self.timestamp,
                self.name,
                despike_col,
                omitted_col,
                self.data["Time"].iloc[self.bkgd_start_idx],
                self.data["Time"].iloc[self.bkgd_stop_idx],
                self.data["Time"].iloc[self.int_start_idx],
                self.data["Time"].iloc[self.int_stop_idx],
                self.int_std,
                np.median(self.bkgd_subtract_data[:, self.int_std_loc]),
            ]
        ).T
        spot_data.columns = [
            "timestamp",
            "Spot",
            "despiked",
            "omitted_region",
            "bkgd_start",
            "bkgd_stop",
            "int_start",
            "int_stop",
            "norm",
            "norm_cps",
        ]
        spot_data = pd.concat(
            [
                spot_data,
                pd.DataFrame(
                    self.bkgd_subtract_med[np.newaxis, :], columns=self.analytes
                ),
                pd.DataFrame(
                    self.bkgd_subtract_std_err_rel[np.newaxis, :],
                    columns=[f"{analyte}_se" for analyte in self.analytes],
                ),
            ],
            axis="columns",
        )

        for col in ["bkgd_start", "bkgd_stop", "int_start", "int_stop", "norm_cps"]:
            spot_data[col] = spot_data[col].astype(np.float64)

        self.output_report = spot_data

    def despike_data(self, analyte_list="all", std_devs=4, window=25):
        """
        despike counts per second normalized to an internal standard using a z score filter

        Parameters
        ----------
        analyte_list : str, optional
            list of analytes to despike. Accepts singular analytes e.g., "29Si"
            or numerous e.g., ["7Li", "29Si"]. by default "all"
        std_devs : int, optional
            number of standard deviations from the mean to be considered an outlier, by default 3
        window : int, optional
            size of the window to be used in the moving average, by default 50
        """

        assert (
            self.bkgd_subtract_normal_data is not None
        ), "please normalize your data prior to despiking"

        self.despiked = True

        if analyte_list == "all":
            filter_list = self.analytes
        else:
            if isinstance(analyte_list, list):
                pass
            else:
                analyte_list = [analyte_list]

            filter_list = analyte_list

        self.despiked_elements = filter_list

        df = pd.DataFrame(self.bkgd_subtract_normal_data, columns=self.analytes)

        for analyte in filter_list:

            filtered = _z_filter(df[analyte], window=window, std_devs=std_devs)

            # replaces data with despiked data
            df[analyte] = filtered

        self.bkgd_subtract_normal_data = df.loc[:, self.analytes].values

        # now recalculate uncertainties after despiking
        # standard error of the mean for the interval region
        self.bkgd_subtract_std_err = self.bkgd_subtract_normal_data.std(
            axis=0
        ) / np.sqrt(abs(self.int_stop_idx - self.int_start_idx))

        self.bkgd_subtract_std_err_rel = 100 * (
            self.bkgd_subtract_std_err / self.bkgd_subtract_med
        )
