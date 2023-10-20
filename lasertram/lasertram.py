"""
    The API for the package `lasertram`
    This will largely be comprised of two classes:

    - `LaserTRAM`: for taking raw counts per second data from a
    Laser Ablation Inductively Coupled Plasma Mass Spectrometry (LA-ICP-MS)
    experiment, choosing an interval to be turned into a concentration, normalizing
    that interval to an internal standard and outputting that value + other metadata

    - `LaserCalc`: for taking the output from `LaserTRAM` along with user input
    to calculate concentrations for a series of `LaserTRAM` spot objects along
    with statistics on calibration standards

    Created and maintained by:
    Jordan Lubbers
    jlubbers@usgs.gov

"""
import re

import mendeleev
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.tools.eval_measures import rmse


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
        self.name = name
        self.despiked = False
        self.despiked_elements = None

    def get_data(self, df):
        """assigns raw counts/sec data to the object

        Args:
            df (pandas DataFrame): raw data corresponding to the spot being processed i.e., `all_data.loc[spot,:]` if `all_data` is the LT_ready file
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
                self.data["Time"][self.omit_start_idx + self.int_start_idx],
                self.data["Time"][self.omit_stop_idx + self.int_start_idx],
            )
        else:
            omitted_col = "None"

        spot_data = pd.DataFrame(
            [
                self.timestamp,
                self.name,
                despike_col,
                omitted_col,
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

        self.output_report = spot_data

    def despike_data(self, analyte_list="all"):
        """
        apply a standard deviation filter to all specified
        analytes.


        Args:
            analyte_list (str or list, optional):analyte to despike (e.g., '7Li'). Or list of analytes to despike (e.g., ['7Li','88Sr']). If 'all', despikes all analytes in the experiment. Defaults to "all".
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
    spot,
    raw_data,
    bkgd,
    keep,
    internal_std,
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
    spot.assign_int_std(internal_std)
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


def oxide_to_ppm(wt_percent, int_std):
    """
    convert concentration internal standard analyte oxide in weight percent to
    concentration ppm for a 1D series of data

    Args:
    wt_percent (array-like): the oxide values to be converted to ppm
    int_std (str): the internal standard used in the experiment (e.g., '29Si', '43Ca','47Ti')

    Returns:
    ppm (array-like): concentrations in ppm the same shape as the wt_percent input

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
    # LaserCalc

    The class `LaserCalc` which is devoted to calculating
    concentrations for laser ablation ICP-MS spot or
    line of spots data following the methodology of
    Longerich et al., (1996) and Kent and Ungerer (2006). It should be used in conjunction
    with the output from `LaserTRAM` class. The basic steps are as follows:

    1. upload SRM data
    2. upload `LaserTRAM` output
    3. set the calibration standard
    4. set the internal standard concentrations for the unknowns
    5. calculate the concentrations and uncertainties of all analyses

    References


    - Longerich, H. P., Jackson, S. E., & Günther, D. (1996). Inter-laboratory note.
            Laser ablation inductively coupled plasma mass spectrometric transient signal
            data acquisition and analyte concentration calculation. Journal of analytical
            atomic spectrometry, 11(9), 899-904.
    - Kent, A. J., & Ungerer, C. A. (2006). Analysis of light lithophile elements
            (Li, Be, B) by laser ablation ICP-MS: comparison between magnetic sector and
            quadrupole ICP-MS. American Mineralogist, 91(8-9), 1401-1411.


    """

    def __init__(self, name):
        """


        Args:
            name (str): The name of the experiment to be processed
        """
        self.name = name

    def get_SRM_comps(self, df):
        """load in a database of standard reference material compositions

        Args:
            df (pandas DataFrame): pandas DataFrame of standard reference materials
        where each row represents data for a standard reference material.
        The first column should be named "Standard". All other columns are
        for different elemental concentrations.Standard names must be exact
        names found in GEOREM: http://georem.mpch-mainz.gwdg.de/sample_query_pref.asp
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
        """load in output from `LaserTRAM` for calculation of concentrations

        Args:
            df (pandas DataFrame): a 2D pandas DataFrame representing numerous concatenated calls to `LaserTRAM.make_output_report()`

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
        analytes_nomass = []
        for i in range(len(self.analytes)):
            # strip the atomic number from our analyte data
            nomass = re.split("(\d+)", self.analytes[i])[2]
            analytes_nomass.append(nomass)
        self.elements = analytes_nomass

    def set_calibration_standard(self, std):
        """Assign which standard reference material will be the calibration
        standard for calculating concentrations.

        Args:
            std (str): name of standard reference material (e.g., `NIST-612`,`BCR-2G`)
        """
        self.calibration_standard = std

        self.calibration_std_data = self.data.loc[std, :]
        # Calibration standard information
        # mean
        self.calibration_std_means = self.calibration_std_data.loc[
            :, self.analytes
        ].mean()
        # std deviation
        self.calibration_std_stdevs = self.calibration_std_data.loc[
            :, self.analytes
        ].std()
        # relative standard error
        self.calibration_std_ses = 100 * (
            (self.calibration_std_stdevs / self.calibration_std_means)
            / np.sqrt(self.calibration_std_data.shape[0])
        )

    def drift_check(self, pval=0.01):
        """For each analyte in the experiment, perform a linear regression to
        assess whether or not drift in the mass spectrometer is happening at a
        significant level. Significance is determined by setting the `pval` threshold.
        If the regression is statistically significant, it gets flagged for later
        correct treatment in `calculate_concentrations`

        """
        calib_std_rmses = []
        calib_std_slopes = []
        calib_std_intercepts = []
        drift_check = []

        # For our calibration standard, calculate the concentration ratio of
        # each analyte to the element used as the internal standard
        std_conc_ratios = []
        myanalytes_nomass = []

        f_pvals = []
        f_vals = []
        f_crits = []
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

            y = self.calibration_std_data[self.analytes[j]].astype("float64")

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
            f_vals.append(fvalue)
            f_pvalue = model.f_pvalue
            f_pvals.append(f_pvalue)
            fcrit = stats.f.ppf(q=1 - pval, dfn=len(x) - 1, dfd=len(y) - 1)
            f_crits.append(fcrit)
            if (f_pvalue < pval) and (fvalue > fcrit):
                drift = "True"
                drift_check.append(drift)
            else:
                drift = "False"
                drift_check.append(drift)

        self.calibration_standard_stats = pd.DataFrame(
            {
                "drift_correct": drift_check,
                "f_pval": f_pvals,
                "f_value": f_vals,
                "f_crit_value": f_crits,
                "rmse": calib_std_rmses,
                "slope": calib_std_slopes,
                "intercept": calib_std_intercepts,
                "mean": self.calibration_std_means[self.analytes].to_numpy(),
                "std_dev": self.calibration_std_stdevs[self.analytes].to_numpy(),
                "percent_std_err": self.calibration_std_ses[self.analytes].to_numpy(),
            },
            index=self.analytes,
        )

    def get_calibration_std_ratios(self):
        """
        For the calibration standard, calculate the concentration ratio between every analyte and the internal standard.
        """

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
        self,
        spots=None,
        concentrations=None,
        uncertainties=None,
    ):
        """Assign the concentration and uncertainty of the internal standard analyte to
        a series of spots.

        Briefly...a linear change in the concentration value reflects a linear change
        in the calculated concentration.

        Args:
            spots (pandas Series): pandas series containing the names of the spots tohave their internal standard concentration-uncertainty assigned. This is the `Spot` column from the output of `LaserTRAM`.

            concentrations (array-like): values representing the internal standard concentration. Must be the same shape as `spots`.
            uncertainties (array-like): values representing the internal standard relative uncertainty in percent. Must be the same shape as `spots`.
        """
        if spots is None:
            spots = (self.data["Spot"],)
            concentrations = (np.full(self.data["Spot"].shape[0], 10),)
            uncertainties = (np.full(self.data["Spot"].shape[0], 1),)

        self.data["internal_std_comp"] = 10
        self.data["internal_std_rel_unc"] = 1
        df = self.data.reset_index().set_index("Spot")

        for spot, concentration, uncertainty in zip(
            spots, concentrations, uncertainties
        ):
            df.loc[spot, "internal_std_comp"] = concentration
            df.loc[spot, "internal_std_rel_unc"] = uncertainty

        self.data["internal_std_comp"] = df["internal_std_comp"].to_numpy()
        self.data["internal_std_rel_unc"] = df["internal_std_rel_unc"].to_numpy()

    def calculate_concentrations(self):
        """
        Calculates the concentration and uncertainty of all spots in the experiment
        using the user specified calibration standard and internal standard
        concentrations/uncertainties.

        """

        secondary_standards = self.potential_calibration_standards.copy()
        secondary_standards.remove(self.calibration_standard)
        self.secondary_standards = secondary_standards
        secondary_standards_concentrations_list = []
        unknown_concentrations_list = []

        for sample in secondary_standards:
            Cn_u = self.standards_data.loc[
                sample,
                re.split(
                    "(\d+)",
                    self.calibration_std_data["norm"].unique()[0],
                )[2],
            ]
            Cin_std = self.calibration_standard_conc_ratios
            Ni_std = self.calibration_standard_stats["mean"][self.analytes]
            Ni_u = self.data.loc[sample, self.analytes]

            concentrations = Cn_u * (Cin_std / Ni_std) * Ni_u

            drift_concentrations_list = []

            for j, analyte, slope, intercept, drift in zip(
                range(len(self.analytes)),
                self.analytes,
                self.calibration_standard_stats["slope"],
                self.calibration_standard_stats["intercept"],
                self.calibration_standard_stats["drift_correct"],
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

                    Ni_std = frac

                    drift_concentrations = Cn_u * (Cin_std[j] / Ni_std) * Ni_u[analyte]

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
            else:
                drift_df = pd.DataFrame()

            for column in drift_df.columns.tolist():
                if type(concentrations) == pd.Series:
                    concentrations.loc[column] = drift_df[column].to_numpy()[0]

                else:
                    concentrations[column] = drift_df[column].to_numpy()

            if type(concentrations) == pd.Series:
                concentrations = pd.DataFrame(concentrations).T
                concentrations["sample"] = sample
                concentrations.set_index("sample", inplace=True)

            secondary_standards_concentrations_list.append(concentrations)

        ###############################
        for sample in self.samples_nostandards:
            Cn_u = oxide_to_ppm(
                self.data.loc[sample, "internal_std_comp"],
                self.data.loc[sample, "norm"].unique()[0],
            ).to_numpy()
            Cin_std = self.calibration_standard_conc_ratios
            Ni_std = self.calibration_standard_stats["mean"][self.analytes].to_numpy()
            Ni_u = self.data.loc[sample, self.analytes].to_numpy()

            concentrations = pd.DataFrame(
                Cn_u[:, np.newaxis] * (Cin_std / Ni_std) * Ni_u, columns=self.analytes
            )

            drift_concentrations_list = []

            for j, analyte, slope, intercept, drift in zip(
                range(len(self.analytes)),
                self.analytes,
                self.calibration_standard_stats["slope"],
                self.calibration_standard_stats["intercept"],
                self.calibration_standard_stats["drift_correct"],
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
                    frac = np.array(frac)
                    drift_concentrations = (
                        Cn_u[:, np.newaxis]
                        * (Cin_std[j] / frac)[:, np.newaxis]
                        * Ni_u[:, j][:, np.newaxis]
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

            for column in drift_df.columns.tolist():
                if type(concentrations) == pd.Series:
                    concentrations.loc[column] = drift_df[column].to_numpy()[0]

                else:
                    concentrations[column] = drift_df[column].to_numpy()

            if type(concentrations) == pd.Series:
                concentrations = pd.DataFrame(concentrations).T
                concentrations["sample"] = sample
                concentrations.set_index("sample", inplace=True)

            unknown_concentrations_list.append(concentrations)

        self.SRM_concentrations = pd.concat(secondary_standards_concentrations_list)
        self.unknown_concentrations = pd.concat(unknown_concentrations_list)

        self.calculate_uncertainties()

        # ADD IN SPOT METADATA NOW

        self.unknown_concentrations[self.unknown_concentrations < 0] = "b.d.l."
        self.SRM_concentrations[self.SRM_concentrations < 0] = "b.d.l."

        self.SRM_concentrations.insert(
            0, "Spot", list(self.data.loc[self.secondary_standards, "Spot"])
        )

        if "timestamp" in self.data.columns.tolist():
            self.SRM_concentrations.insert(
                0,
                "timestamp",
                list(self.data.loc[self.secondary_standards, "timestamp"]),
            )
        else:
            self.SRM_concentrations.insert(
                0, "index", list(self.data.loc[self.secondary_standards, "index"])
            )
        self.unknown_concentrations.insert(
            0, "Spot", list(self.data.loc[self.samples_nostandards, "Spot"])
        )
        if "timestamp" in self.data.columns.tolist():
            self.unknown_concentrations.insert(
                0,
                "timestamp",
                list(self.data.loc[self.samples_nostandards, "timestamp"]),
            )
        else:
            self.unknown_concentrations.insert(
                0, "index", list(self.data.loc[self.samples_nostandards, "index"])
            )

        self.unknown_concentrations.index = [
            "unknown"
        ] * self.unknown_concentrations.shape[0]
        self.unknown_concentrations.index.name = "sample"

    def calculate_uncertainties(self):
        """
        Calculate the uncertainties for each analysis.

        """
        # self.SRM_concentration_uncertainties = 1
        # self.unknown_concentration_uncertainties = 1
        myuncertainties = [analyte + "_se" for analyte in self.analytes]
        srm_rel_uncertainties_list = []
        unk_rel_uncertainties_list = []
        # use RMSE of regression for elements where drift correction is applied rather than the standard error
        # of the mean of all the calibration standard normalized ratios
        rse_i_std = []
        for analyte in self.analytes:
            if "True" in self.calibration_standard_stats.loc[analyte, "drift_correct"]:
                rse_i_std.append(
                    100
                    * self.calibration_standard_stats.loc[analyte, "rmse"]
                    / self.calibration_standard_stats.loc[analyte, "mean"]
                )
            else:
                rse_i_std.append(
                    self.calibration_standard_stats.loc[analyte, "percent_std_err"]
                )

        rse_i_std = np.array(rse_i_std)

        for sample in self.secondary_standards:
            # concentration of internal standard in unknown uncertainties
            int_std_element = re.split(
                "(\d+)", self.calibration_std_data["norm"].unique()[0]
            )[2]
            t1 = (
                self.standards_data.loc[sample, f"{int_std_element}_std"]
                / self.standards_data.loc[sample, f"{int_std_element}"]
            ) ** 2

            # concentration of internal standard in calibration standard uncertainties
            t2 = (
                self.standards_data.loc[
                    self.calibration_standard, f"{int_std_element}_std"
                ]
                / self.standards_data.loc[
                    self.calibration_standard, f"{int_std_element}"
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
                                self.calibration_standard, f"{nomass}_std"
                            ]
                            / self.standards_data.loc[self.calibration_standard, nomass]
                        )
                        ** 2
                    )

            std_conc_stds = np.array(std_conc_stds)

            # Overall uncertainties
            # Need to loop through each row?

            rel_uncertainty = pd.DataFrame(
                np.sqrt(
                    np.array(
                        t1
                        + t2
                        + std_conc_stds
                        + (rse_i_std[np.newaxis, :] / 100) ** 2
                        + (self.data.loc[sample, myuncertainties].to_numpy() / 100) ** 2
                    ).astype(np.float64)
                )
            )
            rel_uncertainty.columns = myuncertainties
            srm_rel_uncertainties_list.append(rel_uncertainty)

        srm_rel_uncertainties = pd.concat(srm_rel_uncertainties_list)

        srm_uncertainties = pd.DataFrame(
            srm_rel_uncertainties.values
            * self.SRM_concentrations.loc[:, self.analytes].values,
            columns=myuncertainties,
            index=self.SRM_concentrations.index,
        )

        self.SRM_concentrations = pd.concat(
            [self.SRM_concentrations, srm_uncertainties], axis="columns"
        )

        ######################################

        for sample in self.samples_nostandards:
            # concentration of internal standard in unknown uncertainties
            int_std_element = re.split(
                "(\d+)", self.calibration_std_data["norm"].unique()[0]
            )[2]
            # concentration of internal standard in unknown uncertainties
            t1 = (self.data.loc[sample, "internal_std_rel_unc"] / 100) ** 2
            t1 = np.array(t1)
            t1 = t1[:, np.newaxis]

            # concentration of internal standard in calibration standard uncertainties
            t2 = (
                self.standards_data.loc[
                    self.calibration_standard, f"{int_std_element}_std"
                ]
                / self.standards_data.loc[
                    self.calibration_standard, f"{int_std_element}"
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
                                self.calibration_standard, f"{nomass}_std"
                            ]
                            / self.standards_data.loc[self.calibration_standard, nomass]
                        )
                        ** 2
                    )

            std_conc_stds = np.array(std_conc_stds)

            # Overall uncertainties
            # Need to loop through each row?

            rel_uncertainty = pd.DataFrame(
                np.sqrt(
                    np.array(
                        t1
                        + t2
                        + std_conc_stds
                        + (rse_i_std[np.newaxis, :] / 100) ** 2
                        + (self.data.loc[sample, myuncertainties].to_numpy() / 100) ** 2
                    ).astype(np.float64)
                )
            )
            rel_uncertainty.columns = myuncertainties
            unk_rel_uncertainties_list.append(rel_uncertainty)

        unk_rel_uncertainties = pd.concat(unk_rel_uncertainties_list)
        overall_uncertainties = (
            unk_rel_uncertainties.values
            * self.unknown_concentrations.loc[:, self.analytes].values
        )
        srm_uncertainties = pd.DataFrame(
            unk_rel_uncertainties.values
            * self.unknown_concentrations.loc[:, self.analytes].values,
            columns=myuncertainties,
            index=self.unknown_concentrations.index,
        )

        self.unknown_concentrations = pd.concat(
            [self.unknown_concentrations, srm_uncertainties], axis="columns"
        )

    # make an accuracy checking function
    # need to use analytes no mass to check SRM vals
    def get_secondary_standard_accuracies(self):
        """
        calculate the accuracy of each secondary standard where accuracy is 100 * measured / accepted value

        Here `accepted` value is the GEOREM preferred value for that SRM analyte pair.

        """
        df_list = []

        nomass = [
            re.split("(\d+)", self.analytes[i])[2] for i in range(len(self.analytes))
        ]

        for standard in self.secondary_standards:
            df = pd.DataFrame(
                100
                * self.SRM_concentrations.loc[standard, self.analytes]
                .replace("b.d.l.", np.nan)
                .values
                / self.standards_data.loc[standard, nomass].values[np.newaxis, :],
                columns=self.analytes,
                index=self.SRM_concentrations.loc[standard, :].index,
            ).fillna("b.d.l.")
            df.insert(0, "Spot", self.SRM_concentrations.loc[standard, "Spot"])
            if "timestamp" in self.data.columns:
                df.insert(
                    0, "timestamp", self.SRM_concentrations.loc[standard, "timestamp"]
                )
            else:
                df.insert(0, "index", self.SRM_concentrations.loc[standard, "index"])

            df_list.append(df)

        self.SRM_accuracies = pd.concat(df_list)
