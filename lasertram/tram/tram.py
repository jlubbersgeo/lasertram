"""

TRAM module: (T)ime (R)esolved (A)nalysis (M)odule.

For taking raw counts per second data from a Laser Ablation Inductively
Coupled Plasma Mass Spectrometry (LA-ICP-MS) experiment, choosing an
interval to be turned into a concentration, normalizing
that interval to an internal standard and outputting that value + other
metadata


"""

import warnings

import numpy as np
import pandas as pd


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
        self.name = name
        self.despiked = False
        self.despiked_elements = None
        self.data = None
        self.data_matrix = None
        self.analytes = None
        self.timestamp = None
        self.int_std = None
        self.bkgd_start = None
        self.bkgd_stop = None
        self.int_start = None
        self.int_stop = None
        self.bkgd_start_idx = None
        self.bkgd_stop_idx = None
        self.int_start_idx = None
        self.int_stop_idx = None
        self.omit_start = None
        self.omit_stop = None
        self.omit_start_idx = None
        self.omit_stop_idx = None
        self.omitted_region = None
        self.bkgd_data = None
        self.detection_limits = None
        self.bkgd_correct_data = None
        self.int_std_loc = None
        self.bkgd_subtract_normal_data = None
        self.bkgd_correct_med = None
        self.bkgd_correct_std_err = None
        self.bkgd_correct_std_err_rel = None
        self.output_report = None

    def get_data(self, df, time_units="ms"):
        """assigns raw counts/sec data to the object

        Args:
            df (pandas DataFrame): raw data corresponding to the spot being processed i.e., `all_data.loc[spot,:]` if `all_data` is the LT_ready file
            time_units (str): string denoting the units for the `Time` column. Used to convert input time values to seconds. Defaults to 'ms'.
        """
        self.data = df.reset_index()
        self.data = self.data.set_index("SampleLabel")
        if time_units == "ms":
            self.data["Time"] = self.data["Time"] / 1000
        elif time_units == "s":
            pass
        self.data_matrix = self.data.iloc[:, 1:].to_numpy()
        self.analytes = self.data.loc[:, "Time":].columns.tolist()[1:]

        # need to add check for if this exists otherwise there is no timestamp attribute
        self.timestamp = str(self.data.loc[:, "timestamp"].unique()[0])

    def assign_int_std(self, int_std):
        """assigns the spot an internal standard
        analyte

        Args:
            int_std (str): the name of the column for the internal standard analyte e.g., "29Si"
        """

        self.int_std = int_std

    def assign_intervals(self, bkgd, keep, omit=None):
        """assigns the intervals to be used as background
        as well as the portion of the ablation interval to
        be used in calculating concentrations

        Args:
            bkgd (tuple): (start, stop) pair of values corresponding to the analysis time where the background signal starts and stops
            keep (tuple): (start, stop) pair of values correpsonding to the analysis time where the interval signal for concentrations starts and stops
            omit (tuple): (start, stop) pair of values corresponding to the analysis time to be omitted from the `keep` interval. Defaults to None.
        """

        self.bkgd_start = bkgd[0]
        self.bkgd_stop = bkgd[1]
        self.int_start = keep[0]
        self.int_stop = keep[1]

        self.bkgd_start_idx = np.where(self.data["Time"] > self.bkgd_start)[0][0]
        self.bkgd_stop_idx = np.where(self.data["Time"] > self.bkgd_stop)[0][0]
        self.int_start_idx = np.where(self.data["Time"] > self.int_start)[0][0]
        self.int_stop_idx = np.where((self.data["Time"] > self.int_stop))[0][0]

        self.omitted_region = False

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

        self.bkgd_data = np.median(
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
        self.bkgd_correct_data = (
            self.data_matrix[self.int_start_idx : self.int_stop_idx, 1:]
            - self.bkgd_data
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
        self.int_std_loc = np.where(np.array(self.analytes) == self.int_std)[0][0]

        threshold = self.detection_limits - np.median(
            self.data_matrix[self.bkgd_start_idx : self.bkgd_stop_idx, 1:], axis=0
        )

        if self.omitted_region is True:
            self.bkgd_subtract_normal_data = np.delete(
                self.bkgd_correct_data,
                np.arange(self.omit_start_idx, self.omit_stop_idx),
                axis=0,
            ) / np.delete(
                self.bkgd_correct_data[:, self.int_std_loc][:, None],
                np.arange(self.omit_start_idx, self.omit_stop_idx),
                axis=0,
            )

        else:
            self.bkgd_subtract_normal_data = (
                self.bkgd_correct_data
                / self.bkgd_correct_data[:, self.int_std_loc][:, None]
            )

        self.bkgd_correct_med = np.median(self.bkgd_subtract_normal_data, axis=0)
        self.bkgd_correct_med[
            np.median(self.bkgd_correct_data, axis=0) <= threshold
        ] = -9999
        self.bkgd_correct_med[np.median(self.bkgd_correct_data, axis=0) == 0] = -9999

        self.bkgd_correct_std_err = self.bkgd_subtract_normal_data.std(
            axis=0
        ) / np.sqrt(abs(self.int_stop_idx - self.int_start_idx))
        self.bkgd_correct_std_err_rel = 100 * (
            self.bkgd_correct_std_err / self.bkgd_correct_med
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
                np.median(self.bkgd_correct_data[:, self.int_std_loc]),
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
                    self.bkgd_correct_med[np.newaxis, :], columns=self.analytes
                ),
                pd.DataFrame(
                    self.bkgd_correct_std_err_rel[np.newaxis, :],
                    columns=[f"{analyte}_se" for analyte in self.analytes],
                ),
            ],
            axis="columns",
        )

        for col in ["bkgd_start", "bkgd_stop", "int_start", "int_stop", "norm_cps"]:
            spot_data[col] = spot_data[col].astype(np.float64)

        self.output_report = spot_data

    def despike_data(self, analyte_list="all"):
        """
        apply a standard deviation filter to all specified
        analytes.


        Args:
            analyte_list (str or list, optional): analyte to despike (e.g., '7Li'). Or list of analytes to despike (e.g., ['7Li','88Sr']). If 'all', despikes all analytes in the experiment. Defaults to "all".
        """

        def despike_signal(data, analyte, passes=2):
            """
            apply a standard deviation filter to analyte signal

            Args:
                data (pandas DataFrame): dataframe representing the spot raw counts per second data.
                analyte (string): analyte to despike
                passes (int, optional): the number of iterations for the filter to complete. Defaults to 2.

            Returns:
                signal (ndarray): the filtered signal
            """
            window = 3
            sigma = 25
            kernel = np.ones(window) / window

            signal_raw = data[analyte].to_numpy()
            signal = signal_raw.copy()

            for i in range(passes):
                signal_mean = np.convolve(signal, kernel, "valid")
                signal_mean = np.insert(
                    signal_mean,
                    0,
                    signal_mean[0],
                )
                signal_mean = np.append(signal_mean, signal_mean[-1])
                signal_std = np.sqrt(signal_mean)

                spikes = signal > signal_mean + signal_std * sigma
                despiked_signal = signal.copy()
                despiked_signal[spikes] = signal_mean[spikes]
                signal = despiked_signal

            return signal

        self.despiked = True

        # analyte_list = "all"  # currently only supports despiking every element

        if analyte_list == "all":
            filter_list = self.analytes
        else:
            warnings.warn(
                "single element despiking currently not supported, falling back to 'all'"
            )
            # # this
            # if analyte_list is not type(list):
            #     filter_list = [analyte_list]
            # else:
            #     filter_list = analyte_list
            filter_list = self.analytes

        self.despiked_elements = filter_list
        despiked = []
        for analyte in filter_list:
            despiked.append(despike_signal(self.data, analyte))

        despiked = pd.DataFrame(
            np.array(despiked).T, columns=self.analytes, index=self.data.index
        )
        despiked.index.name = self.data.index.name
        despiked.insert(0, "Time", self.data["Time"])

        self.data = despiked
        self.data_matrix = despiked.to_numpy()
