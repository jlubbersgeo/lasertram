"""
    The API for LaserTRAM-DB
    This will largely be comprised of two classes:
    - LaserTRAM
    - LaserCalc

    Created and maintained by:
    Jordan Lubbers
    jlubbers@usgs.gov

"""
import os
import re

import mendeleev
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.tools.eval_measures import rmse


class LaserTRAM:
    """
    The class `LaserTRAM` which is devoted to the "time resolved analysis"
    operations during thelaser data reduction process. To be used in
    conjunction with the `LaserCalc` class. The general idea is that
    this creates an object that contains all the information related
    to one individual spot analysis. For relatability we can call this
    object `spot`:

    ```python
    spot = LaserTRAM(name = 'BCR-2G')

    # assign data to the spot
    spot.get_data(raw_data)

    # despike the data if desired
    if despike is True:
        spot.despike_data(analyte_list="all")

    # assign the internal standard analyte
    spot.assign_int_std(internal_std)

    # assign intervals for background and ablation signal
    spot.assign_intervals(bkgd=bkgd, interval=keep)

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

    # or utilizing the process_data() function that calls on many of the
    # methods within the LaserTRAM class. Ex:

    #create the spot object
    spot = LaserTRAM(name = 'BCR-2G')
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

    """

    def __init__(self, name):
        """

        Args:
            name (str): your sample name
            i.e., the value in the "SampleLabel" column
            of the LT_ready file
        """
        self.name = name
        self.despiked = False
        self.despiked_elements = None

    def get_data(self, df):
        """assigns raw counts/sec data to the object

        Args:
            df (pandas DataFrame): raw data corresponding
            to the spot being processed i.e.,
            `all_data.loc[spot,:]` if `all_data` is the
            LT_ready file
        """
        self.data = df.reset_index()
        self.data = self.data.set_index("SampleLabel")
        self.data["Time"] = self.data["Time"] / 1000
        self.data_matrix = self.data.iloc[:, 1:].to_numpy()
        self.analytes = self.data.loc[:, "Time":].columns.tolist()[1:]
        self.timestamp = str(self.data.loc[:, "timestamp"].unique()[0])

    def assign_int_std(self, int_std):
        """assigns the spot an internal standard
        analyte

        Args:
            int_std (str): the name of the column for the
            internal standard analyte e.g., "29Si"
        """
        self.int_std = int_std

    def assign_intervals(self, bkgd, keep):
        """assigns the intervals to be used as background
        as well as the portion of the ablation interval to
        be used in calculating concentrations

        Args:
            bkgd (tuple): (start, stop) pair of values corresponding
            to the analysis time where the background signal starts
            and stops
            interval (tuple): (start, stop) pair of values correpsonding
            to the analysis time where the interval signal for concentrations
            starts and stops
        """
        self.bkgd_start = bkgd[0]
        self.bkgd_stop = bkgd[1]
        self.int_start = keep[0]
        self.int_stop = keep[1]

        self.bkgd_start_idx = np.where(self.data["Time"] > self.bkgd_start)[0][0]
        self.bkgd_stop_idx = np.where(self.data["Time"] > self.bkgd_stop)[0][0]
        self.int_start_idx = np.where(self.data["Time"] > self.int_start)[0][0]
        self.int_stop_idx = np.where((self.data["Time"] > self.int_stop))[0][0]

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
        for e, i in zip(self.analytes, range(len(self.analytes))):
            if e == self.int_std:
                break
            self.int_std_loc = i

        threshold = self.detection_limits - np.median(
            self.data_matrix[self.bkgd_start_idx : self.bkgd_stop_idx, 1:], axis=0
        )

        self.bkgd_subtract_normal_data = self.bkgd_correct_data / self.bkgd_correct_data[:, self.int_std_loc][:, None]
        
        self.bkgd_correct_med = np.median(self.bkgd_subtract_normal_data, axis=0)
        self.bkgd_correct_med[
            np.median(self.bkgd_correct_data, axis=0) <= threshold
        ] = -9999
        self.bkgd_correct_med[np.median(self.bkgd_correct_data, axis=0) == 0] = -9999

        self.bkgd_correct_std_err = self.bkgd_subtract_normal_data.std(axis=0) / np.sqrt(
            abs(self.int_stop_idx - self.int_start_idx)
        )
        self.bkgd_correct_std_err_rel = 100 * (
            self.bkgd_correct_std_err / self.bkgd_correct_med
        )

    def make_output_report(self):
        """
        create an output report for the spot processing. This is a
        pandas dataframe that has the following format:

        |timestamp|Spot|bkgd_start|bkgd_stop|int_start|int_stop|norm|norm_cps|analyte vals and uncertainties -->|
        |---------|----|----------|---------|---------|--------|----|--------|----------------------------------|
        """
        if self.despiked is True:
            despike_col = self.despiked_elements
        else:
            despike_col = "None"
        spot_data = pd.DataFrame(
            [
                self.timestamp,
                self.name,
                despike_col,
                self.data["Time"][self.bkgd_start_idx],
                self.data["Time"][self.bkgd_stop_idx],
                self.data["Time"][self.int_start_idx],
                self.data["Time"][self.int_stop_idx],
                self.int_std,
                np.median(self.bkgd_correct_data[:, self.int_std_loc]),
            ]
        ).T
        spot_data.columns = [
            "timestamp",
            "Spot",
            "despiked",
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

        self.output_report = spot_data

    def despike_data(self, analyte_list="all"):
        """
        apply a standard deviation filter to all specified
        analytes.
        """

        def despike_signal(data, analyte, passes=2):
            """
            apply a standard deviation filter to analyte signal

            Args:
                data (pandas DataFrame): dataframe representing the
                spot raw counts per second data.
                analyte (string): analyte to despike
                passes (int, optional): the number of iterations
                for the filter to complete. Defaults to 2.

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

        if analyte_list == "all":
            filter_list = self.analytes
        else:
            if analyte_list is not type(list):
                filter_list = [analyte_list]
            else:
                filter_list = analyte_list

        self.despiked_elements = filter_list
        for analyte in filter_list:
            despiked.append(despike_signal(self.data, analyte))

        despiked = pd.DataFrame(np.array(despiked).T, columns=self.analytes)
        despiked.insert(0, "Time", self.data["Time"])

        self.data = despiked
        self.data_matrix = despiked.to_numpy()


def process_spot(
    spot, raw_data, bkgd, keep, internal_std, despike=False, output_report=True
):
    """a function to incorporate all the methods of the `LaserTRAM` class
    so a spot can be processed in an efficient and compact way.

    Args:
        spot (_type_): _description_
        raw_data (_type_): _description_
        bkgd (_type_): _description_
        keep (_type_): _description_
        internal_std (_type_): _description_
        output_report (bool, optional): _description_. Defaults to True.

    """
    # assign data to the spot
    spot.get_data(raw_data)
    # despike the data if desired
    if despike is True:
        spot.despike_data(analyte_list="all")
    # assign the internal standard analyte
    spot.assign_int_std(internal_std)
    # assign intervals for background and ablation signal
    spot.assign_intervals(bkgd=bkgd, keep=keep)
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


def oxide_to_ppm(wt_percent, int_std):
    """
    convert concentration internal standard analyte oxide in weight percent to
    concentration ppm for a 1D series of data

    Parameters
    ----------
    wt_percent : array-like
        the oxide values to be converted to ppm
    int_std : string
        the internal standard used in the experiment (e.g., '29Si', '43Ca',
                                                      '47Ti')
        currently supported elements:

            - 'SiO2'
            - 'TiO2'
            - 'Al2O3'
            - 'Cr2O3'
            - 'MnO'
            - 'FeO'
            - 'K2O'
            - 'CaO'
            - 'Na2O'
            - 'NiO'
            - 'MgO'

    Returns
    -------
    ppm : array-like
        concentrations in ppm the same shape as the wt_percent input

    """

    el = [i for i in int_std if not i.isdigit()]

    if len(el) == 2:
        element = el[0] + el[1]

    else:
        element = el[0]

    oxides = [
        "SiO2",
        "TiO2",
        "Al2O3",
        "Cr2O3",
        "MnO",
        "FeO",
        "K2O",
        "CaO",
        "Na2O",
        "NiO",
        "MgO",
    ]

    for o in oxides:
        if element in o:
            oxide = o

    s = oxide.split("O")
    cat_subscript = s[0]
    an_subscript = s[1]

    cat_subscript = [i for i in cat_subscript if i.isdigit()]
    if cat_subscript:
        cat_subscript = int(cat_subscript[0])
    else:
        cat_subscript = 1

    an_subscript = [i for i in an_subscript if i.isdigit()]
    if an_subscript:
        an_subscript = int(an_subscript[0])
    else:
        an_subscript = 1

    ppm = 1e4 * (
        (wt_percent * mendeleev.element(element).atomic_weight * cat_subscript)
        / (
            mendeleev.element(element).atomic_weight
            + mendeleev.element("O").atomic_weight * an_subscript
        )
    )
    return ppm


class LaserCalc:
    """
    The class `LaserCalc` which is devoted to calculating
    concentrations for laser ablation ICP-MS spot or
    line of spots data following the methodology of
    Longerich et al., (1996). It should be used in conjunction
    with the output from `LaserTRAM` class. The basic steps
    are as follows:
    1. upload SRM data
    2. upload `LaserTRAM` output
    3. set the calibration standard
    4. set the internal standard concentrations for the unknowns
    5. calculate the concentrations of the unknowns
    """

    def __init__(self, name):
        """


        Args:
            name (str): The name of the experiment
            to be processed
        """
        self.name = name

    def get_SRM_comps(self, df):
        """_summary_

        Args:
            df (_type_): _description_
        """

        self.standards_data = df.set_index("Standard")
        self.database_standards = self.standards_data.index.unique()
        # Get a list of all of the elements supported in the published standard datasheet
        # Get a second list for the same elements but their corresponding uncertainty columns
        self.standard_elements = [
            analyte
            for analyte in self.standards_data.columns.tolist()
            if not ("_std" in analyte)
        ]
        self.standard_element_uncertainties = [
            analyte + "_std" for analyte in self.standard_elements
        ]

    def get_data(self, df):
        """_summary_

        Args:
            df (_type_): _description_
        """

        data = df.set_index("Spot")
        data.insert(loc=1, column="index", value=np.arange(1, len(data) + 1))

        self.spots = data.index.unique().tolist()

        # Check for potential calibration standards. This will let us know what our options
        # are for choosing calibration standards by looking for spots that have the same string
        # as the standard spreadsheet

        stds_column = [
            [std for std in self.database_standards if std in spot]
            for spot in self.spots
        ]

        stds_column = [["unknown"] if not l else l for l in stds_column]

        stds_column = [std for sublist in stds_column for std in sublist]

        # standards that can be used as calibrations standards (must have more than 1 analysis)
        # potential_standards = list(np.unique(stds_column))
        potential_standards = [
            std for std in np.unique(stds_column) if stds_column.count(std) > 1
        ]
        potential_standards.remove("unknown")

        # all of the samples in your input sheet that are NOT potential standards
        all_standards = list(np.unique(stds_column))
        all_standards.remove("unknown")

        data["sample"] = stds_column

        data.reset_index(inplace=True)
        data.set_index("sample", inplace=True)

        self.data = data
        self.potential_calibration_standards = potential_standards
        self.samples_nostandards = list(np.setdiff1d(stds_column, all_standards))

        self.analytes = [
            analyte
            for analyte in data.columns.tolist()
            if not (
                "_se" in analyte
                or "norm" in analyte
                or "index" in analyte
                or "Spot" in analyte
                or "wt%" in analyte
                or "1stdev%" in analyte
                or "start" in analyte
                or "stop" in analyte
                or "long" in analyte
                or "timestamp" in analyte
                or "despiked" in analyte
            )
        ]

    def set_calibration_standard(self, std):
        """_summary_

        Args:
            std (_type_): _description_
        """
        self.calibration_standard = std

        self.calibration_std_data = self.data.loc[std, :]
        # Calibration standard information
        # mean
        self.calibration_std_means = self.calibration_std_data.loc[
            :, self.analytes + [analyte + "_se" for analyte in self.analytes]
        ].mean()
        # std deviation
        self.calibration_std_stdevs = self.calibration_std_data.loc[
            :, self.analytes + [analyte + "_se" for analyte in self.analytes]
        ].std()
        # relative standard error
        self.calibration_std_ses = 100 * (
            (self.calib_std_stdevs / self.calibration_std_means)
            / np.sqrt(self.calibration_std_data.shape[0])
        )

    def drift_check(self, pval=0.01):
        """_summary_"""
        calib_std_rmses = []
        calib_std_slopes = []
        calib_std_intercepts = []
        drift_check = []

        # For our calibration standard, calculate the concentration ratio of each analyte to the element used as the internal standard
        std_conc_ratios = []
        myanalytes_nomass = []
        for j in range(len(self.analytes)):
            # Getting regression statistics on analyte normalized ratios through time
            # for the calibration standard. This is what we use to check to see if it needs
            # to be drift corrected
            if "timestamp" in self.calibration_std_data.columns.tolist():
                # get an array in time units based on timestamp column. This is
                # is in seconds
                x = np.array(
                    [
                        np.datetime64(d, "m")
                        for d in self.calibration_std_data["timestamp"]
                    ]
                ).astype(np.float64)
                # x = np.cumsum(np.diff(x))
                # x = np.insert(x, 0, 0).astype(np.float64)

            else:
                x = self.calibration_std_data["index"]

            # x = self.self.calibration_std_data["index"]
            y = self.calibration_std_data[self.analytes[j]]

            X = sm.add_constant(x)
            # Note the difference in argument order
            model = sm.OLS(y, X).fit()
            # now generate predictions
            ypred = model.predict(X)

            # calc rmse
            RMSE = rmse(y, ypred)

            calib_std_rmses.append(RMSE)

            if model.params.shape[0] < 2:
                calib_std_slopes.append(model.params[0])
                calib_std_intercepts.append(0)

            else:
                calib_std_slopes.append(model.params[1])
                calib_std_intercepts.append(model.params[0])

            # new stuff
            # confidence limit 99%

            # f value stuff
            fvalue = model.fvalue
            f_pvalue = model.f_pvalue
            fcrit = stats.f.ppf(q=1 - pval, dfn=len(x) - 1, dfd=len(y) - 1)
            if (f_pvalue < pval) and (fvalue > fcrit):
                drift = "True"
                drift_check.append(drift)
            else:
                drift = "False"
                drift_check.append(drift)

        self.calibration_standard_stats = pd.DataFrame(
            {
                "drift_correct": drift_check,
                "rmse": calib_std_rmses,
                "slope": calib_std_slopes,
                "intercept": calib_std_intercepts,
            }
        )

    def get_calibration_std_ratios(self):
        # For our calibration standard, calculate the concentration ratio
        # of each analyte to the element used as the internal standard
        std_conc_ratios = []
        myanalytes_nomass = []

        for i in range(len(self.analytes)):
            # strip the atomic number from our analyte data
            nomass = re.split("(\d+)", self.analytes[i])[2]
            # make it a list
            myanalytes_nomass.append(nomass)

            # if our element is in the list of standard elements take the ratio
            if nomass in self.standard_elements:
                std_conc_ratios.append(
                    self.standards_data.loc[self.calibration_standard, nomass]
                    / self.standards_data.loc[
                        self.calibration_standard,
                        re.split(
                            "(\d+)", self.calibration_std_data["norm"].unique()[0]
                        )[2],
                    ]
                )

        # make our list an array for easier math going forward
        # std_conc_ratios = pd.DataFrame(np.array(std_conc_ratios)[np.newaxis,:],columns = myanalytes)
        self.calibration_standard_conc_ratios = np.array(std_conc_ratios)

    def set_internal_standard_concentrations(
        self, spots, concentrations, uncertainties
    ):
        """_summary_

        Args:
            spots (_type_): _description_
            concentrations (_type_): _description_
            uncertainties (_type_): _description_
        """
        self.data["internal_std_comp"] = 10
        self.data["internal_std_rel_unc"] = 1
        for spot, concentration, uncertainty in zip(
            spots, concentrations, uncertainties
        ):
            self.data.loc[spot, "internal_std_comp"] = concentration
            self.data.loc[spot, "internal_std_rel_unc"] = uncertainty

    def calculate_concentrations(self, SRM=False):
        # use this one function to calculate concentrations based on whether
        # it's an SRM or an unknown...functionally the same!

        secondary_standards = self.potential_calibration_standards.copy()
        secondary_standards.remove(self.calibration_standard)
        concentrations_list = []

        myuncertainties = [analyte + "_se" for analyte in self.analytes]

        if SRM is True:
            sample_list = secondary_standards
        else:
            sample_list = self.samples_nostandards

            for sample in sample_list:
                drift_concentrations_list = []

                for j, analyte, slope, intercept, drift in zip(
                    range(len(self.analytes)),
                    self.analytes,
                    self.calibration_standard_stats["slope"],
                    self.calibration_standard_stats["intercept"],
                    self.calibration_standard_stats["drift_check"],
                ):
                    if "True" in drift:
                        if "timestamp" in self.data.columns.tolist():
                            frac = (
                                slope
                                * np.array(
                                    [
                                        np.datetime64(d, "m")
                                        for d in self.data.loc[sample, "timestamp"]
                                    ]
                                ).astype(np.float64)
                                + intercept
                            )
                        else:
                            frac = slope * self.data.loc[sample, "index"] + intercept

                        if SRM is True:
                            drift_concentrations = (
                                (
                                    self.standards_data.loc[
                                        sample,
                                        re.split(
                                            "(\d+)",
                                            self.calibration_std_data["norm"].unique()[
                                                0
                                            ],
                                        )[2],
                                    ]
                                )
                                * (self.calibration_standard_conc_ratios[j] / frac)
                                * self.data.loc[sample, analyte]
                            )
                        else:
                            drift_concentrations = (
                                self.data.loc[sample, analyte]
                                * (self.calibration_standard_conc_ratios[j] / frac)
                                * (self.set_internal_standard_concentrations)
                            )

                        if type(drift_concentrations) == np.float64:
                            df = pd.DataFrame(
                                np.array([drift_concentrations]), columns=[analyte]
                            )

                        else:
                            df = pd.DataFrame(drift_concentrations, columns=[analyte])

                        drift_concentrations_list.append(df)

                if len(drift_concentrations_list) > 0:
                    drift_df = pd.concat(drift_concentrations_list, axis="columns")

                    if drift_df.shape[0] == 1:
                        drift_df["sample"] = sample
                        drift_df.set_index("sample", inplace=True)

                    concentrations = (
                        (
                            self.standards_data.loc[
                                sample,
                                re.split(
                                    "(\d+)",
                                    self.calibration_std_data["norm"].unique()[0],
                                )[2],
                            ]
                        )
                        * (
                            self.calibration_standard_conc_ratios
                            / self.calibration_standard_stats["mean"][self.analytes]
                        )
                        * self.data.loc[sample, self.analytes]
                    )

                    for column in drift_df.columns.tolist():
                        if type(concentrations) == pd.Series:
                            concentrations.loc[column] = drift_df[column].to_numpy()[0]

                        else:
                            concentrations[column] = drift_df[column]

                    if type(concentrations) == pd.Series:
                        concentrations = pd.DataFrame(concentrations).T
                        concentrations["sample"] = sample
                        concentrations.set_index("sample", inplace=True)

                    concentrations_list.append(concentrations)
                else:
                    concentrations = (
                        (
                            self.standards_data.loc[
                                sample,
                                re.split(
                                    "(\d+)",
                                    self.calibration_std_data["norm"].unique()[0],
                                )[2],
                            ]
                        )
                        * (
                            self.calibration_standard_conc_ratios
                            / self.calibration_standard_stats["mean"][self.analytes]
                        )
                        * self.data.loc[sample, self.analytes]
                    )
                    concentrations_list.append(concentrations)
            if SRM is True:
                self.SRM_concentrations = concentrations_list

            else:
                self.unknown_concentrations = concentrations_list

        # incorporate uncertainty in calibration standard
        calib_uncertainty = True

        stds_list = []
        unknowns_list = []
        # relative uncertainty in % of the concentration of the internal standard
        # in the unknown
        # unknown_int_std_unc = 1

        # use RMSE of regression for elements where drift correction is applied rather than the standard error
        # of the mean of all the calibration standard normalized ratios
        for j in range(len(self.calibration_standard_stats["drift_check"])):
            if "True" in self.calibration_standard_stats["drift_check"][j]:
                self.calibration_standard_stats["mean"][j] = (
                    100
                    * self.calibration_standard_stats["rmse"][j]
                    / self.calibration_standard_stats["mean"][j]
                )
        if SRM is True:
            # creates a list of dataframes that hold the uncertainty information for each secondary standard.
            for standard, concentration in zip(
                secondary_standards, concentrations_list
            ):
                # concentration of internal standard in unknown uncertainties
                t1 = (
                    self.standards_data.loc[
                        standard,
                        "{}_std".format(
                            re.split(
                                "(\d+)", self.calibration_std_data["norm"].unique()[0]
                            )[2]
                        ),
                    ]
                    / self.standards_data.loc[
                        standard,
                        "{}".format(
                            re.split(
                                "(\d+)", self.calibration_std_data["norm"].unique()[0]
                            )[2]
                        ),
                    ]
                ) ** 2

                # concentration of internal standard in calibration standard uncertainties
                t2 = (
                    self.standards_data.loc[
                        self.calibration_standard,
                        "{}_std".format(
                            re.split(
                                "(\d+)", self.calibration_std_data["norm"].unique()[0]
                            )[2]
                        ),
                    ]
                    / self.standards_data.loc[
                        self.calibration_standard,
                        "{}".format(
                            re.split(
                                "(\d+)",
                                self.self.calibration_std_data["norm"].unique()[0],
                            )[2]
                        ),
                    ]
                ) ** 2

                # concentration of each analyte in calibration standard uncertainties
                std_conc_stds = []
                for i in range(len(self.analytes)):
                    # strip the atomic number from our analyte data
                    nomass = re.split("(\d+)", self.analytes[i])[2]

                    # if our element is in the list of standard elements take the ratio
                    if nomass in self.standard_elements:
                        std_conc_stds.append(
                            (
                                self.standards_data.loc[
                                    self.calibration_standard, "{}_std".format(nomass)
                                ]
                                / self.standards_data.loc[
                                    self.calibration_standard, nomass
                                ]
                            )
                            ** 2
                        )

                std_conc_stds = np.array(std_conc_stds)

                # Overall uncertainties

                if calib_uncertainty == True:
                    stds_values = concentration * np.sqrt(
                        np.array(
                            t1
                            + t2
                            + std_conc_stds
                            + (
                                self.calibration_std_ses[self.analytes].to_numpy()[
                                    np.newaxis, :
                                ]
                                / 100
                            )
                            ** 2
                            + (
                                self.data.loc[standard, myuncertainties].to_numpy()
                                / 100
                            )
                            ** 2
                        ).astype(np.float64)
                    )

                    stds_values.columns = myuncertainties

                    stds_list.append(stds_values)
                else:
                    stds_values = concentration * np.sqrt(
                        t1
                        + t2
                        + std_conc_stds
                        + (
                            self.calibration_std_ses[self.analytes].to_numpy()[
                                np.newaxis, :
                            ]
                            / 100
                        )
                        ** 2
                        + (self.loc[standard, myuncertainties].to_numpy() / 100) ** 2
                    )
                    stds_values.columns = myuncertainties
                    stds_list.append(stds_values)

            self.SRM_uncertainties = stds_list

        else:
            # creates a list of dataframes that hold the uncertainty information for unknown spot.
            for sample, concentration in zip(
                self.samples_nostandards, concentrations_list
            ):
                # concentration of internal standard in unknown uncertainties
                t1 = (self.data["internal_std_rel_unc"] / 100) ** 2
                t1 = t1[:, np.newaxis]

                # concentration of internal standard in calibration standard uncertainties
                t2 = (
                    self.standards_data.loc[
                        self.calibration_standard,
                        "{}_std".format(
                            re.split(
                                "(\d+)",
                                self.self.calibration_std_data["norm"].unique()[0],
                            )[2]
                        ),
                    ]
                    / self.standards_data.loc[
                        self.calibration_standard,
                        "{}".format(
                            re.split(
                                "(\d+)",
                                self.self.calibration_std_data["norm"].unique()[0],
                            )[2]
                        ),
                    ]
                ) ** 2

                # concentration of each analyte in calibration standard uncertainties
                std_conc_stds = []
                for i in range(len(self.analytes)):
                    # strip the atomic number from our analyte data
                    nomass = re.split("(\d+)", self.analytes[i])[2]

                    # if our element is in the list of standard elements take the ratio
                    if nomass in self.standard_elements:
                        std_conc_stds.append(
                            (
                                self.standards_data.loc[
                                    self.calibration_standard, "{}_std".format(nomass)
                                ]
                                / self.standards_data.loc[
                                    self.calibration_standard, nomass
                                ]
                            )
                            ** 2
                        )

                std_conc_stds = np.array(std_conc_stds)

                if calib_uncertainty == True:
                    unknown_stds_values = concentration * np.sqrt(
                        t1
                        + t2
                        + std_conc_stds
                        + (
                            self.calibration_std_ses[self.analytes].to_numpy()[
                                np.newaxis, :
                            ]
                            / 100
                        )
                        ** 2
                        + (self.data.loc[sample, myuncertainties].to_numpy() / 100) ** 2
                    )
                    unknown_stds_values.columns = myuncertainties
                    unknowns_list.append(unknown_stds_values)
                else:
                    unknown_stds_values = concentration * np.sqrt(
                        t2
                        + std_conc_stds
                        + (
                            self.calibration_std_ses[self.analytes].to_numpy()[
                                np.newaxis, :
                            ]
                            / 100
                        )
                        ** 2
                        + (self.data.loc[sample, myuncertainties].to_numpy() / 100) ** 2
                    )
                    unknown_stds_values.columns = myuncertainties
                    unknowns_list.append(unknown_stds_values)

            self.unknown_uncertainties = unknowns_list
