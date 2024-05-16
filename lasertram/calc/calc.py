"""

Calc module:
This is for taking the output from `tram.LaserTRAM` along with user input
to calculate concentrations for a series of `LaserTRAM` spot objects along
with statistics on calibration standards

"""

import re

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.tools.eval_measures import rmse

from ..helpers import conversions


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
        # all attributes in relative chronological order that they are created in
        # if everything is done correctly. These all will get rewritten throughout the
        # data processing pipeline but this allows us to see what all the potential attributes
        # are going to be from the beginning (PEP convention)
        self.name = name
        self.standards_data = None
        self.database_standards = None
        self.standard_elements = None
        self.standard_element_uncertainties = None
        self.spots = None
        self.analytes = None
        self.calibration_std_stdevs = None
        self.data = None
        self.int_std_element = None
        self.potential_calibration_standards = None
        self.samples_nostandards = None
        self.elements = None
        self.calibration_std = None
        self.calibration_std_data = None
        self.calibration_std_means = None
        self.calibration_std_ses = None
        self.calibration_std_stats = None
        self.calibration_std_conc_ratios = None
        self.secondary_standards = None
        self.SRM_concentrations = None
        self.unknown_concentrations = None
        self.SRM_accuracies = None

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
        self.database_standards = self.standards_data.index.unique().to_list()
        # Get a list of all of the elements supported in the published standard datasheet
        # Get a second list for the same elements but their corresponding uncertainty columns
        self.standard_elements = [
            analyte
            for analyte in self.standards_data.columns.tolist()
            if "_std" not in analyte
        ]
        self.standard_element_uncertainties = [
            analyte + "_std" for analyte in self.standard_elements
        ]

    def get_data(self, df):
        """load in output from `LaserTRAM` for calculation of concentrations

        Args:
            df (pandas DataFrame): a 2D pandas DataFrame representing numerous concatenated calls to `LaserTRAM.make_output_report()`

        """
        # check if first row is nan (output from GUI does this).
        # If so, remove it
        df = df[df.iloc[:, 0].isna() is False]

        data = df.set_index("Spot")
        data.insert(loc=1, column="index", value=np.arange(1, len(data) + 1))

        self.spots = data.index.unique().dropna().tolist()

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
                or "omitted_region" in analyte
            )
        ]
        # elements without isotopes in the front
        self.elements = [re.split(r"(\d+)", analyte)[2] for analyte in self.analytes]

        # internal standard analyte from lasertram
        self.int_std_element = re.split(r"(\d+)", self.data["norm"].unique()[0])[2]

    def set_calibration_standard(self, std):
        """Assign which standard reference material will be the calibration
        standard for calculating concentrations.

        Args:
            std (str): name of standard reference material (e.g., `NIST-612`,`BCR-2G`)
        """
        self.calibration_std = std

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

        f_pvals = []
        f_vals = []
        f_crits = []
        for analyte in self.analytes:
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
                x = self.calibration_std_data["index"].to_numpy()

            y = self.calibration_std_data.loc[:, analyte].astype("float64")

            X = sm.add_constant(x)
            # Note the difference in argument order
            model = sm.OLS(y, X).fit()
            # now generate predictions
            ypred = model.predict(X)

            # calc rmse
            RMSE = rmse(y, ypred)

            calib_std_rmses.append(RMSE)

            if model.params.shape[0] < 2:
                calib_std_slopes.append(model.params.loc["x1"])
                calib_std_intercepts.append(0)

            else:
                calib_std_slopes.append(model.params.loc["x1"])
                calib_std_intercepts.append(model.params.loc["const"])

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

        self.calibration_std_stats = pd.DataFrame(
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

        for element in self.elements:
            if element in self.standard_elements:
                std_conc_ratios.append(
                    self.standards_data.loc[self.calibration_std, element]
                    / self.standards_data.loc[
                        self.calibration_std, self.int_std_element
                    ]
                )

        # make our list an array for easier math going forward
        self.calibration_std_conc_ratios = np.array(std_conc_ratios)

    def set_int_std_concentrations(
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

        self.data["int_std_comp"] = 10.0
        self.data["int_std_rel_unc"] = 1.0
        df = self.data.reset_index().set_index("Spot")

        for spot, concentration, uncertainty in zip(
            spots, concentrations, uncertainties
        ):
            df.loc[spot, "int_std_comp"] = concentration
            df.loc[spot, "int_std_rel_unc"] = uncertainty

        self.data["int_std_comp"] = df["int_std_comp"].to_numpy()
        self.data["int_std_rel_unc"] = df["int_std_rel_unc"].to_numpy()

    def calculate_concentrations(self):
        """
        Calculates the concentration and uncertainty of all spots in the experiment
        using the user specified calibration standard and internal standard
        concentrations/uncertainties.

        """

        secondary_standards = self.potential_calibration_standards.copy()
        secondary_standards.remove(self.calibration_std)
        self.secondary_standards = secondary_standards
        secondary_standards_concentrations_list = []
        unknown_concentrations_list = []

        for sample in secondary_standards:
            Cn_u = self.standards_data.loc[
                sample,
                re.split(
                    r"(\d+)",
                    self.calibration_std_data["norm"].unique()[0],
                )[2],
            ]
            Cin_std = self.calibration_std_conc_ratios
            Ni_std = self.calibration_std_stats["mean"][self.analytes]
            Ni_u = self.data.loc[sample, self.analytes]

            concentrations = Cn_u * (Cin_std / Ni_std) * Ni_u

            drift_concentrations_list = []

            for j, analyte, slope, intercept, drift in zip(
                range(len(self.analytes)),
                self.analytes,
                self.calibration_std_stats["slope"],
                self.calibration_std_stats["intercept"],
                self.calibration_std_stats["drift_correct"],
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

                    if isinstance(drift_concentrations, np.float64):
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
                if isinstance(concentrations, pd.Series):
                    concentrations.loc[column] = drift_df[column].to_numpy()[0]

                else:
                    concentrations[column] = drift_df[column].to_numpy()

            if isinstance(concentrations, pd.Series):
                concentrations = pd.DataFrame(concentrations).T
                concentrations["sample"] = sample
                concentrations.set_index("sample", inplace=True)

            secondary_standards_concentrations_list.append(concentrations)

        ###############################
        for sample in self.samples_nostandards:
            Cn_u = conversions.oxide_to_ppm(
                self.data.loc[sample, "int_std_comp"],
                self.data.loc[sample, "norm"].unique()[0],
            ).to_numpy()
            Cin_std = self.calibration_std_conc_ratios
            Ni_std = self.calibration_std_stats["mean"][self.analytes].to_numpy()
            Ni_u = self.data.loc[sample, self.analytes].to_numpy()

            concentrations = pd.DataFrame(
                Cn_u[:, np.newaxis] * (Cin_std / Ni_std) * Ni_u, columns=self.analytes
            )

            drift_concentrations_list = []

            for j, analyte, slope, intercept, drift in zip(
                range(len(self.analytes)),
                self.analytes,
                self.calibration_std_stats["slope"],
                self.calibration_std_stats["intercept"],
                self.calibration_std_stats["drift_correct"],
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

                    if isinstance(drift_concentrations, np.float64):
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
                if isinstance(concentrations, pd.Series):
                    concentrations.loc[column] = drift_df[column].to_numpy()[0]

                else:
                    concentrations[column] = drift_df[column].to_numpy()

            if isinstance(concentrations, pd.Series):
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

        myuncertainties = [analyte + "_se" for analyte in self.analytes]
        srm_rel_ext_uncertainties_list = []
        unk_rel_ext_uncertainties_list = []
        srm_rel_int_uncertainties_list = []
        unk_rel_int_uncertainties_list = []
        # use RMSE of regression for elements where drift correction is applied rather than the standard error
        # of the mean of all the calibration standard normalized ratios
        rse_i_std = []
        for analyte in self.analytes:
            if "True" in self.calibration_std_stats.loc[analyte, "drift_correct"]:
                rse_i_std.append(
                    100
                    * self.calibration_std_stats.loc[analyte, "rmse"]
                    / self.calibration_std_stats.loc[analyte, "mean"]
                )
            else:
                rse_i_std.append(
                    self.calibration_std_stats.loc[analyte, "percent_std_err"]
                )

        rse_i_std = np.array(rse_i_std)

        for sample in self.secondary_standards:
            t1 = (
                self.standards_data.loc[sample, f"{self.int_std_element}_std"]
                / self.standards_data.loc[sample, f"{self.int_std_element}"]
            ) ** 2

            # concentration of internal standard in calibration standard uncertainties
            t2 = (
                self.standards_data.loc[
                    self.calibration_std, f"{self.int_std_element}_std"
                ]
                / self.standards_data.loc[
                    self.calibration_std, f"{self.int_std_element}"
                ]
            ) ** 2

            # concentration of each analyte in calibration standard uncertainties
            std_conc_stds = []
            for element in self.elements:
                # if our element is in the list of standard elements take the ratio
                if element in self.standard_elements:
                    std_conc_stds.append(
                        (
                            self.standards_data.loc[
                                self.calibration_std, f"{element}_std"
                            ]
                            / self.standards_data.loc[self.calibration_std, element]
                        )
                        ** 2
                    )

            std_conc_stds = np.array(std_conc_stds)

            # Overall uncertainties
            # Need to loop through each row?

            rel_ext_uncertainty = pd.DataFrame(
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
            rel_int_uncertainty = pd.DataFrame(
                np.sqrt(
                    np.array(
                        t1
                        # +t2
                        # + std_conc_stds
                        + (rse_i_std[np.newaxis, :] / 100) ** 2
                        + (self.data.loc[sample, myuncertainties].to_numpy() / 100) ** 2
                    ).astype(np.float64)
                )
            )
            rel_ext_uncertainty.columns = [f"{a}_exterr" for a in self.analytes]
            srm_rel_ext_uncertainties_list.append(rel_ext_uncertainty)
            rel_int_uncertainty.columns = [f"{a}_interr" for a in self.analytes]
            srm_rel_int_uncertainties_list.append(rel_int_uncertainty)

        srm_rel_ext_uncertainties = pd.concat(srm_rel_ext_uncertainties_list)
        srm_rel_int_uncertainties = pd.concat(srm_rel_int_uncertainties_list)

        srm_ext_uncertainties = pd.DataFrame(
            srm_rel_ext_uncertainties.values
            * self.SRM_concentrations.loc[:, self.analytes].values,
            columns=[f"{a}_exterr" for a in self.analytes],
            index=self.SRM_concentrations.index,
        )
        srm_int_uncertainties = pd.DataFrame(
            srm_rel_int_uncertainties.values
            * self.SRM_concentrations.loc[:, self.analytes].values,
            columns=[f"{a}_interr" for a in self.analytes],
            index=self.SRM_concentrations.index,
        )

        self.SRM_concentrations = pd.concat(
            [self.SRM_concentrations, srm_ext_uncertainties, srm_int_uncertainties],
            axis="columns",
        )

        ######################################

        for sample in self.samples_nostandards:
            # concentration of internal standard in unknown uncertainties
            int_std_element = re.split(
                r"(\d+)", self.calibration_std_data["norm"].unique()[0]
            )[2]
            # concentration of internal standard in unknown uncertainties
            t1 = (self.data.loc[sample, "int_std_rel_unc"] / 100) ** 2
            t1 = np.array(t1)
            t1 = t1[:, np.newaxis]

            # concentration of internal standard in calibration standard uncertainties
            t2 = (
                self.standards_data.loc[self.calibration_std, f"{int_std_element}_std"]
                / self.standards_data.loc[self.calibration_std, f"{int_std_element}"]
            ) ** 2

            # concentration of each analyte in calibration standard uncertainties
            std_conc_stds = []
            for element in self.elements:
                # # if our element is in the list of standard elements take the ratio
                if element in self.standard_elements:
                    std_conc_stds.append(
                        (
                            self.standards_data.loc[
                                self.calibration_std, f"{element}_std"
                            ]
                            / self.standards_data.loc[self.calibration_std, element]
                        )
                        ** 2
                    )

            std_conc_stds = np.array(std_conc_stds)

            # Overall uncertainties
            # Need to loop through each row?

            rel_ext_uncertainty = pd.DataFrame(
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
            rel_int_uncertainty = pd.DataFrame(
                np.sqrt(
                    np.array(
                        t1
                        # +t2
                        # + std_conc_stds
                        + (rse_i_std[np.newaxis, :] / 100) ** 2
                        + (self.data.loc[sample, myuncertainties].to_numpy() / 100) ** 2
                    ).astype(np.float64)
                )
            )
            rel_ext_uncertainty.columns = [f"{a}_exterr" for a in self.analytes]
            unk_rel_ext_uncertainties_list.append(rel_ext_uncertainty)
            rel_int_uncertainty.columns = [f"{a}_interr" for a in self.analytes]
            unk_rel_int_uncertainties_list.append(rel_int_uncertainty)

        unk_rel_ext_uncertainties = pd.concat(unk_rel_ext_uncertainties_list)
        unk_rel_int_uncertainties = pd.concat(unk_rel_int_uncertainties_list)

        unknown_ext_uncertainties = pd.DataFrame(
            unk_rel_ext_uncertainties.values
            * self.unknown_concentrations.loc[:, self.analytes].values,
            columns=[f"{a}_exterr" for a in self.analytes],
            index=self.unknown_concentrations.index,
        )

        unknown_int_uncertainties = pd.DataFrame(
            unk_rel_int_uncertainties.values
            * self.unknown_concentrations.loc[:, self.analytes].values,
            columns=[f"{a}_interr" for a in self.analytes],
            index=self.unknown_concentrations.index,
        )

        self.unknown_concentrations = pd.concat(
            [
                self.unknown_concentrations,
                unknown_ext_uncertainties,
                unknown_int_uncertainties,
            ],
            axis="columns",
        )

    # make an accuracy checking function
    # need to use analytes no mass to check SRM vals
    def get_secondary_standard_accuracies(self):
        """
        calculate the accuracy of each secondary standard where accuracy is 100 * measured / accepted value

        Here `accepted` value is the GEOREM preferred value for that SRM analyte pair.

        """
        df_list = []

        for standard in self.secondary_standards:
            df = pd.DataFrame(
                100
                * self.SRM_concentrations.loc[standard, self.analytes]
                .replace("b.d.l.", np.nan)
                .values
                / self.standards_data.loc[standard, self.elements].values[
                    np.newaxis, :
                ],
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
