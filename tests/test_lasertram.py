"""
various tests for the package lasertram
"""
import numpy as np
import pandas as pd
import pytest

from lasertram import lasertram as lt
from lasertram.lasertram import LaserCalc, LaserTRAM

###########LASERTRAM UNIT TESTS##############
spreadsheet_path = r"tests\spot_test_timestamp_raw_data.xlsx"


@pytest.fixture
def load_data():
    data = pd.read_excel(spreadsheet_path).set_index("SampleLabel")
    return data


def test_get_data(load_data):
    """
    checks whether or not data are loaded in properly
    """
    spot = LaserTRAM(name="test")

    samples = load_data.index.unique().dropna().tolist()
    spot.get_data(load_data.loc[samples[0], :])
    df_to_check = spot.data.copy()

    df_to_check["Time"] = df_to_check["Time"] * 1000

    # check to see if input data are the same as the data stored in the lasertram object
    # all other attributes created will be correct if this is correct
    pd.testing.assert_frame_equal(df_to_check, load_data.loc[samples[0], :])


def test_assign_int_std(load_data):
    """
    test that the internal standard is set correctly
    """
    spot = LaserTRAM(name="test")

    samples = load_data.index.unique().dropna().tolist()

    spot.get_data(load_data.loc[samples[0], :])

    spot.assign_int_std("29Si")

    assert spot.int_std == "29Si", "the internal standard should be '29Si'"


def test_assign_intervals(load_data):
    """
    test that the intervals are assigned correctly
    """

    spot = LaserTRAM(name="test")

    samples = load_data.index.unique().dropna().tolist()

    spot.get_data(load_data.loc[samples[0], :])

    spot.assign_int_std("29Si")

    bkgd_interval = (5, 10)
    keep_interval = (20, 50)
    omit_interval = (30, 35)

    spot.assign_intervals(bkgd=bkgd_interval, keep=keep_interval, omit=omit_interval)

    assert spot.bkgd_start == bkgd_interval[0], "the bkgd_start should be 5"
    assert spot.bkgd_stop == bkgd_interval[1], "the bkgd_stop should be 10"
    assert spot.int_start == keep_interval[0], "the int_start should be 20"
    assert spot.int_stop == keep_interval[1], "the int_stop should be 50"
    assert spot.omit_start == omit_interval[0], "the omit_start should be 30"
    assert spot.omit_stop == omit_interval[1], "the omit_stop should be 35"
    assert spot.omitted_region is True, "omittted_region should be True"


def test_get_bkgd_data(load_data):
    """
    test that background signal is being assigned properly
    """

    spot = LaserTRAM(name="test")

    samples = load_data.index.unique().dropna().tolist()

    spot.get_data(load_data.loc[samples[0], :])

    spot.assign_int_std("29Si")

    bkgd_interval = (5, 10)
    keep_interval = (20, 50)
    omit_interval = (30, 35)

    spot.assign_intervals(bkgd=bkgd_interval, keep=keep_interval, omit=omit_interval)
    spot.get_bkgd_data()

    assert np.allclose(
        spot.bkgd_data,
        np.array(
            [
                700.01960055,
                100.0004,
                200.00160001,
                43575.82193016,
                100.0004,
                0.0,
                900.03240117,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ]
        ),
    ), "background values are not correctly assigned"


def test_subtract_bkgd(load_data):
    """
    test that the background signal is correctly subtracted
    from the interval data
    """

    spot = LaserTRAM(name="test")

    samples = load_data.index.unique().dropna().tolist()

    spot.get_data(load_data.loc[samples[0], :])

    spot.assign_int_std("29Si")

    bkgd_interval = (5, 10)
    keep_interval = (20, 50)
    omit_interval = (30, 35)

    spot.assign_intervals(bkgd=bkgd_interval, keep=keep_interval, omit=omit_interval)
    spot.get_bkgd_data()
    spot.subtract_bkgd()

    assert np.allclose(
        spot.bkgd_correct_data,
        spot.data_matrix[spot.int_start_idx : spot.int_stop_idx, 1:] - spot.bkgd_data,
    ), "background not subtracted properly"


def test_get_detection_limits(load_data):
    """
    test to make sure detection limits are generated correctly
    """

    spot = LaserTRAM(name="test")

    samples = load_data.index.unique().dropna().tolist()

    spot.get_data(load_data.loc[samples[0], :])

    spot.assign_int_std("29Si")

    bkgd_interval = (5, 10)
    keep_interval = (20, 50)
    omit_interval = (30, 35)

    spot.assign_intervals(bkgd=bkgd_interval, keep=keep_interval, omit=omit_interval)
    spot.get_bkgd_data()
    spot.get_detection_limits()

    assert np.allclose(
        spot.detection_limits,
        np.array(
            [
                1472.42001196,
                421.41658043,
                727.91181812,
                49692.17439946,
                321.93969037,
                336.49074757,
                1839.41436852,
                0.0,
                71.58217609,
                0.0,
                51.42615343,
                51.42615343,
                287.19571818,
            ]
        ),
    ), "detection limits not calculated correctly"


def test_normalize_interval(load_data):
    """
    check that data are being normalized correctly
    """
    spot = LaserTRAM(name="test")

    samples = load_data.index.unique().dropna().tolist()

    spot.get_data(load_data.loc[samples[0], :])

    spot.assign_int_std("29Si")

    bkgd_interval = (5, 10)
    keep_interval = (20, 50)
    omit_interval = (30, 35)

    spot.assign_intervals(bkgd=bkgd_interval, keep=keep_interval, omit=omit_interval)
    spot.get_bkgd_data()
    spot.subtract_bkgd()
    spot.get_detection_limits()
    spot.normalize_interval()
    assert spot.bkgd_subtract_normal_data.shape[0] == (
        spot.int_stop_idx - spot.int_start_idx
    ) - (
        spot.omit_stop_idx - spot.omit_start_idx
    ), "background subtracted and normalized data is not the right shape. Likely a region omission problem"

    assert np.allclose(
        spot.bkgd_correct_med,
        np.array(
            [
                1.84185261e-03,
                3.20170958e00,
                1.27326069e01,
                1.00000000e00,
                3.62879351e-02,
                4.09905508e00,
                1.26082054e00,
                3.33186264e-01,
                7.27569058e-01,
                2.89757096e-02,
                6.68289728e-02,
                1.82759066e-03,
                9.43237495e-03,
            ]
        ),
    ), "median background and normalized values are incorrect"
    assert np.allclose(
        spot.bkgd_correct_std_err_rel,
        np.array(
            [
                100.99950523,
                3.04429518,
                3.06514828,
                0.0,
                6.68695225,
                3.52178191,
                3.26529918,
                1.96909904,
                4.56046366,
                3.9863298,
                4.35410981,
                20.62679912,
                17.99393343,
            ]
        ),
    ), "standard error values are incorrect"


def test_make_output_report(load_data):
    """
    check to make sure output report is generated correctly
    """

    spot = LaserTRAM(name="test")

    samples = load_data.index.unique().dropna().tolist()

    spot.get_data(load_data.loc[samples[0], :])

    spot.assign_int_std("29Si")

    bkgd_interval = (5, 10)
    keep_interval = (20, 50)
    omit_interval = (30, 35)

    spot.assign_intervals(bkgd=bkgd_interval, keep=keep_interval, omit=omit_interval)
    spot.get_bkgd_data()
    spot.subtract_bkgd()
    spot.get_detection_limits()
    spot.normalize_interval()
    spot.make_output_report()

    pd.testing.assert_frame_equal(
        spot.output_report,
        pd.DataFrame(
            {
                "timestamp": {0: "2021-03-01 22:08:14"},
                "Spot": {0: "test"},
                "despiked": {0: "None"},
                "omitted_region": {0: (30.07824, 35.039379999999994)},
                "bkgd_start": {0: 5.12408},
                "bkgd_stop": {0: 10.084719999999999},
                "int_start": {0: 20.00647},
                "int_stop": {0: 50.07226},
                "norm": {0: "29Si"},
                "norm_cps": {0: 2499024.199695497},
                "7Li": {0: 0.0018418526122406418},
                "24Mg": {0: 3.2017095770096136},
                "27Al": {0: 12.732606932805151},
                "29Si": {0: 1.0},
                "43Ca": {0: 0.0362879351388429},
                "48Ti": {0: 4.099055080072567},
                "57Fe": {0: 1.260820538158862},
                "88Sr": {0: 0.3331862637258016},
                "138Ba": {0: 0.7275690579268954},
                "139La": {0: 0.0289757095722955},
                "140Ce": {0: 0.06682897275601189},
                "153Eu": {0: 0.0018275906598832652},
                "208Pb": {0: 0.009432374952294733},
                "7Li_se": {0: 100.99950522922774},
                "24Mg_se": {0: 3.044295180539666},
                "27Al_se": {0: 3.0651482778196257},
                "29Si_se": {0: 0.0},
                "43Ca_se": {0: 6.6869522521261295},
                "48Ti_se": {0: 3.5217819050704238},
                "57Fe_se": {0: 3.2652991786108005},
                "88Sr_se": {0: 1.9690990404248054},
                "138Ba_se": {0: 4.560463662642775},
                "139La_se": {0: 3.986329799770403},
                "140Ce_se": {0: 4.354109813192227},
                "153Eu_se": {0: 20.62679912145815},
                "208Pb_se": {0: 17.993933433607573},
            }
        ),
    )


def test_process_spot(load_data):
    """
    check to see if the process_spot helper function produces same output
    as doing calculations one by one in LaserTRAM

    """

    spot = LaserTRAM(name="test")

    samples = load_data.index.unique().dropna().tolist()

    spot.get_data(load_data.loc[samples[0], :])

    spot.assign_int_std("29Si")

    bkgd_interval = (5, 10)
    keep_interval = (20, 50)
    omit_interval = (30, 35)

    spot.assign_intervals(bkgd=bkgd_interval, keep=keep_interval, omit=omit_interval)
    spot.get_bkgd_data()
    spot.subtract_bkgd()
    spot.get_detection_limits()
    spot.normalize_interval()
    spot.make_output_report()

    spot2 = LaserTRAM(name="test")
    lt.process_spot(
        spot2,
        raw_data=load_data.loc[samples[0], :],
        bkgd=bkgd_interval,
        keep=keep_interval,
        omit=omit_interval,
        internal_std="29Si",
        despike=False,
        output_report=True,
    )

    pd.testing.assert_frame_equal(spot.output_report, spot2.output_report)


def test_oxide_to_ppm():
    """
    test that oxides wt% are being converted to elemental ppm
    properly. Test
    """

    analytes = ["29Si", "27Al", "CaO"]
    oxide_vals = np.array([65.0, 15.0, 8.0])

    result = {}
    for analyte, oxide in zip(analytes, oxide_vals):
        result[analyte] = lt.oxide_to_ppm(oxide, analyte)

    expected = {
        "29Si": 303833.8631559676,
        "27Al": 107957.04626864659,
        "CaO": 34304.891110317745,
    }

    assert (
        result == expected
    ), "concentrations from oxides not being calculated properly"


###############LASERCALC UNIT TESTS#########################
SRM_path = r"tests\laicpms_stds_tidy.xlsx"


@pytest.fixture
def load_SRM_data():
    data = pd.read_excel(SRM_path)
    return data


LT_complete_path = r"tests\spot_test_timestamp_lasertram_complete.xlsx"


@pytest.fixture
def load_LTcomplete_data():
    data = pd.read_excel(LT_complete_path)
    return data


def test_get_SRM_comps(load_SRM_data):
    concentrations = LaserCalc(name="test")
    concentrations.get_SRM_comps(load_SRM_data)

    assert concentrations.standard_elements == [
        "Ag",
        "Al",
        "As",
        "Au",
        "B",
        "Ba",
        "Be",
        "Bi",
        "Br",
        "Ca",
        "Cd",
        "Ce",
        "Cl",
        "Co",
        "Cr",
        "Cs",
        "Cu",
        "Dy",
        "Er",
        "Eu",
        "F",
        "Fe",
        "Ga",
        "Gd",
        "Ge",
        "Hf",
        "Ho",
        "In",
        "K",
        "La",
        "Li",
        "Lu",
        "Mg",
        "Mn",
        "Mo",
        "Na",
        "Nb",
        "Nd",
        "Ni",
        "P",
        "Pb",
        "Pr",
        "Rb",
        "Re",
        "S",
        "Sb",
        "Sc",
        "Se",
        "Si",
        "Sm",
        "Sn",
        "Sr",
        "Ta",
        "Tb",
        "Th",
        "Ti",
        "Tl",
        "Tm",
        "U",
        "V",
        "W",
        "Y",
        "Yb",
        "Zn",
        "Zr",
        "SiO2",
        "TiO2",
        "Sl2O3",
        "FeO",
        "MgO",
        "MnO",
        "CaO",
        "Na2O",
        "K2O",
        "P2O5",
    ], "standard elements not being accessed properly"
    assert concentrations.database_standards == [
        "BCR-2G",
        "BHVO-2G",
        "BIR-1G",
        "GSA-1G",
        "GSC-1G",
        "GSD-1G",
        "GSE-1G",
        "NIST-610",
        "NIST-612",
        "BM9021-G",
        "GOR128-G",
        "GOR132-G",
        "ATHO-G",
        "KL2-G",
        "ML3B-G",
        "T1-G",
        "StHs680-G",
    ], "standard names not being read in properly"


def test_get_data(load_SRM_data, load_LTcomplete_data):
    concentrations = LaserCalc(name="test")
    concentrations.get_SRM_comps(load_SRM_data)
    concentrations.get_data(load_LTcomplete_data)

    assert concentrations.spots == [
        "BCR-2G_1",
        "BCR-2G_2",
        "ATHO-G_1",
        "ATHO-G_2",
        "BHVO-2G_1",
        "BHVO-2G_2",
        "unknown_nist-612_1",
        "unknown_nist-612_2",
        "BCR-2G_3",
        "ATHO-G_3",
        "BCR-2G_4",
        "ATHO-G_4",
        "BCR-2G_5",
        "ATHO-G_5",
        "BCR-2G_6",
        "ATHO-G_6",
        "BCR-2G_7",
        "ATHO-G_7",
        "BCR-2G_8",
        "ATHO-G_8",
        "BHVO-2G_3",
        "unknown_nist-612_3",
        "BCR-2G_9",
        "ATHO-G_9",
        "BCR-2G_10",
        "ATHO-G_10",
        "BCR-2G_11",
        "ATHO-G_11",
        "BCR-2G_12",
        "ATHO-G_12",
        "BCR-2G_13",
        "ATHO-G_13",
        "BCR-2G_14",
        "BCR-2G_15",
        "ATHO-G_14",
        "ATHO-G_15",
        "BHVO-2G_4",
        "BHVO-2G_5",
        "unknown_nist-612_4",
        "unknown_nist-612_5",
    ], "analysis spots not found correctly"
    assert concentrations.potential_calibration_standards == [
        "ATHO-G",
        "BCR-2G",
        "BHVO-2G",
    ], "potential calibration standards not found correctly"
    assert concentrations.samples_nostandards == [
        "unknown"
    ], "unknown analyses not found correctly"
    assert concentrations.elements == [
        "Li",
        "Mg",
        "Al",
        "Si",
        "Ca",
        "Ti",
        "Fe",
        "Sr",
        "Ba",
        "La",
        "Ce",
        "Eu",
        "Pb",
    ], "analyte to element conversion not correct"


def test_set_calibration_standard(load_SRM_data, load_LTcomplete_data):
    """
    test whether or not calibration standard data is properly assigned
    """
    concentrations = LaserCalc(name="test")
    concentrations.get_SRM_comps(load_SRM_data)
    concentrations.get_data(load_LTcomplete_data)
    concentrations.set_calibration_standard("BCR-2G")
    test_means = pd.Series(
        {
            "7Li": 0.05095309567485837,
            "24Mg": 85.83785949961276,
            "27Al": 337.3903538168668,
            "29Si": 26.146972020239836,
            "43Ca": 1.0,
            "48Ti": 112.35818866678622,
            "57Fe": 33.27331379718277,
            "88Sr": 9.571563094542862,
            "138Ba": 21.393678102682294,
            "139La": 0.8326420933512372,
            "140Ce": 1.9742222315043418,
            "153Eu": 0.05251069362909363,
            "208Pb": 0.2672836117732711,
        }
    )
    test_ses = pd.Series(
        {
            "7Li": 0.558565709206677,
            "24Mg": 0.38138022061850435,
            "27Al": 0.4573559099297444,
            "29Si": 0.6994253098217633,
            "43Ca": 0.0,
            "48Ti": 0.1988979755580142,
            "57Fe": 0.5621273015720097,
            "88Sr": 0.49500279644282247,
            "138Ba": 0.90397731151757,
            "139La": 0.7938673732133003,
            "140Ce": 0.9547499642078431,
            "153Eu": 0.7394818175678296,
            "208Pb": 0.8311922109557456,
        }
    )

    assert concentrations.calibration_std == "BCR-2G"
    pd.testing.assert_series_equal(concentrations.calibration_std_means, test_means)
    pd.testing.assert_series_equal(concentrations.calibration_std_ses, test_ses)


def test_drift_check(load_SRM_data, load_LTcomplete_data):
    """
    test whether or not drift is accounted for properly
    """
    concentrations = LaserCalc(name="test")
    concentrations.get_SRM_comps(load_SRM_data)
    concentrations.get_data(load_LTcomplete_data)
    concentrations.set_calibration_standard("BCR-2G")
    concentrations.drift_check()
    test_ses = pd.Series(
        {
            "7Li": "False",
            "24Mg": "True",
            "27Al": "True",
            "29Si": "True",
            "43Ca": "False",
            "48Ti": "False",
            "57Fe": "True",
            "88Sr": "True",
            "138Ba": "True",
            "139La": "True",
            "140Ce": "True",
            "153Eu": "True",
            "208Pb": "True",
        }
    )
    test_ses.name = "drift_correct"
    pd.testing.assert_series_equal(
        test_ses, concentrations.calibration_std_stats["drift_correct"]
    ), "analytes not being drift corrected properly"


def test_get_calibration_std_ratios(load_SRM_data, load_LTcomplete_data):
    """
    test that the concentration ratio between every analyte and the internal
    standard is accurate
    """

    concentrations = LaserCalc(name="test")
    concentrations.get_SRM_comps(load_SRM_data)
    concentrations.get_data(load_LTcomplete_data)
    concentrations.set_calibration_standard("BCR-2G")
    concentrations.drift_check()
    concentrations.get_calibration_std_ratios()

    test_ratios = np.array(
        [
            1.78368272e-04,
            4.25488849e-01,
            1.40550695e00,
            5.03990796e00,
            1.00000000e00,
            2.79443626e-01,
            1.91023584e00,
            6.77799433e-03,
            1.35361700e-02,
            4.89521813e-04,
            1.05633654e-03,
            3.90428329e-05,
            2.18005666e-04,
        ]
    )
    assert np.allclose(
        concentrations.calibration_std_conc_ratios, test_ratios
    ), "calibration standard concentration ratios are not correct, check again"


def test_set_internal_standard_concentrations(load_SRM_data, load_LTcomplete_data):
    """
    test to make sure concentration of the internal standard is being set correctly
    """

    concentrations = LaserCalc(name="test")
    concentrations.get_SRM_comps(load_SRM_data)
    concentrations.get_data(load_LTcomplete_data)
    concentrations.set_calibration_standard("BCR-2G")
    concentrations.drift_check()
    concentrations.get_calibration_std_ratios()
    concentrations.set_internal_standard_concentrations(
        concentrations.data["Spot"],
        np.full(concentrations.data["Spot"].shape[0], 71.9),
        np.full(concentrations.data["Spot"].shape[0], 1),
    )

    assert np.allclose(
        concentrations.data.loc["unknown", "internal_std_comp"].values,
        np.array([71.9, 71.9, 71.9, 71.9, 71.9]),
    ), "internal standard concentrations for unknowns not set properly"
    assert np.allclose(
        concentrations.data.loc["unknown", "internal_std_rel_unc"].values,
        np.array([1.0, 1.0, 1.0, 1.0, 1.0]),
    ), "internal standard concentration uncertainties for unknowns not set properly"


def test_calculate_concentrations(load_SRM_data, load_LTcomplete_data):
    """
    test to make sure concentrations are calculated correctly
    """

    concentrations = LaserCalc(name="test")
    concentrations.get_SRM_comps(load_SRM_data)
    concentrations.get_data(load_LTcomplete_data)
    concentrations.set_calibration_standard("BCR-2G")
    concentrations.drift_check()
    concentrations.get_calibration_std_ratios()
    concentrations.set_internal_standard_concentrations(
        concentrations.data["Spot"],
        np.full(concentrations.data["Spot"].shape[0], 71.9),
        np.full(concentrations.data["Spot"].shape[0], 1),
    )
    concentrations.calculate_concentrations()

    test_unknown_concentrations = pd.DataFrame(
        {
            "timestamp": [
                "2021-03-01T22:15:56.000000000",
                "2021-03-01T22:17:12.999000000",
                "2021-03-02T00:32:12.000000000",
                "2021-03-02T02:44:04.000000000",
                "2021-03-02T02:45:22.000000000",
            ],
            "Spot": [
                "unknown_nist-612_1",
                "unknown_nist-612_2",
                "unknown_nist-612_3",
                "unknown_nist-612_4",
                "unknown_nist-612_5",
            ],
            "7Li": [
                207.58205449833375,
                209.06241656705257,
                211.661619810332,
                208.1141728704961,
                208.0021823891214,
            ],
            "24Mg": [
                375.9041094947969,
                368.5288462830524,
                369.7103681340728,
                373.2245442918232,
                372.68091572356707,
            ],
            "27Al": [
                62599.11382390356,
                62700.18216934566,
                61134.56508097334,
                61535.64836974937,
                61137.35962927165,
            ],
            "29Si": [
                1898256.5921452672,
                1889902.4124392755,
                1833782.6965243507,
                1848995.4811800483,
                1869701.7686653552,
            ],
            "43Ca": [
                513866.3266579882,
                513866.3266579882,
                513866.3266579882,
                513866.3266579882,
                513866.3266579882,
            ],
            "48Ti": [
                2611.316937677469,
                2584.862677980655,
                2600.625450456699,
                2585.311614429285,
                2556.404072442201,
            ],
            "57Fe": [
                327.57890580947793,
                314.7478981350732,
                300.9029235545256,
                295.7365807741927,
                277.1823289762299,
            ],
            "88Sr": [
                489.9462645872854,
                485.6654338716087,
                485.97056026250135,
                487.3907248754375,
                482.4571478802297,
            ],
            "138Ba": [
                234.4110430915169,
                232.68028900645848,
                235.88527773280674,
                234.27585578953685,
                232.09812648377985,
            ],
            "139La": [
                229.6609961756181,
                227.92641774625443,
                229.8565286717672,
                227.77742247999876,
                222.81901611170542,
            ],
            "140Ce": [
                231.311641090213,
                230.99105056495816,
                230.5293159759712,
                228.57724688795787,
                226.674801608158,
            ],
            "153Eu": [
                226.74379802104772,
                227.00096869676517,
                225.86829618815372,
                224.68524555796338,
                223.83010332189676,
            ],
            "208Pb": [
                214.4793293705701,
                212.19698209171742,
                216.92719783397186,
                219.03746703576175,
                216.22128744737935,
            ],
            "7Li_exterr": [
                23.457200886517615,
                23.624411347548133,
                23.926147326430055,
                23.530038289080682,
                23.515734417054933,
            ],
            "24Mg_exterr": [
                32.46670630987515,
                12.135084373929159,
                12.120561573992989,
                12.223570745692397,
                12.268494815718407,
            ],
            "27Al_exterr": [
                2254.387393435266,
                2258.8756945878913,
                2200.035012020141,
                2217.6046613996673,
                2203.3401755653963,
            ],
            "29Si_exterr": [
                45271.52823314575,
                45110.42262240762,
                43859.95135867492,
                44384.44662899679,
                44615.7297535611,
            ],
            "43Ca_exterr": [
                12434.277214052305,
                12434.277214052305,
                12434.277214052305,
                12434.277214052305,
                12434.277214052305,
            ],
            "48Ti_exterr": [
                192.14934992254877,
                189.7426953800166,
                190.8929253916731,
                189.81761514248365,
                187.67929976424563,
            ],
            "57Fe_exterr": [
                23.129641254856494,
                12.76043772087462,
                12.845814287713937,
                58.23344513283139,
                11.83564142404549,
            ],
            "88Sr_exterr": [
                11.13267716088216,
                11.032762905909145,
                11.048774527275869,
                11.07035361410948,
                10.98119522372554,
            ],
            "138Ba_exterr": [
                5.3603206434413115,
                5.308217164973913,
                5.381494143181964,
                5.34827669044824,
                5.304421473728012,
            ],
            "139La_exterr": [
                5.476858396651021,
                5.418392408993304,
                5.4570269986760085,
                5.417558954516482,
                5.303434339700193,
            ],
            "140Ce_exterr": [
                5.176215023317281,
                5.18273167788242,
                5.1499298169528895,
                5.122494906686536,
                5.064344113079541,
            ],
            "153Eu_exterr": [
                5.605400613215481,
                5.634739949477251,
                5.567809986743265,
                5.586420849604132,
                5.537225826449885,
            ],
            "208Pb_exterr": [
                20.230964078384577,
                20.01912074319983,
                20.45566332189426,
                20.65594898164351,
                20.391685788599492,
            ],
            "7Li_interr": [
                2.7929500734570385,
                2.812250412873164,
                2.9138315363385776,
                2.9048994887221524,
                2.8900067537222123,
            ],
            "24Mg_interr": [
                30.487267565710695,
                5.242959048441498,
                5.134894429888841,
                5.154840966545983,
                5.294361939158086,
            ],
            "27Al_interr": [
                799.4965646295867,
                803.1768786683676,
                776.241916989768,
                790.1787974557354,
                785.3139284723417,
            ],
            "29Si_interr": [
                31303.96498987536,
                31221.319589825514,
                30422.755946824735,
                30906.27978909062,
                30869.511074353813,
            ],
            "43Ca_interr": [
                5138.663266579882,
                5138.663266579882,
                5138.663266579882,
                5138.663266579882,
                5138.663266579882,
            ],
            "48Ti_interr": [
                31.0973133112473,
                27.798381312201826,
                27.92116286749934,
                28.08822251795146,
                27.666646301418236,
            ],
            "57Fe_interr": [
                21.12155250897871,
                8.988508914353561,
                9.488825419337601,
                57.60824183695335,
                8.74416590259588,
            ],
            "88Sr_interr": [
                5.7293693392208525,
                5.674170246570485,
                5.695370249336902,
                5.691211060140753,
                5.678018841436546,
            ],
            "138Ba_interr": [
                3.1019304270067374,
                3.0573309385853595,
                3.099721837214165,
                3.0846450812505357,
                3.0661207613212813,
            ],
            "139La_interr": [
                3.0677446818399963,
                3.013938989272737,
                3.026409825723158,
                3.016835425074584,
                2.957996423473382,
            ],
            "140Ce_interr": [
                3.0159353989302295,
                3.0351925628413006,
                2.990643456578927,
                2.9930845044933188,
                2.9415388242470546,
            ],
            "153Eu_interr": [
                3.6933096100864877,
                3.7322860951727392,
                3.6548005038614386,
                3.7080308239162694,
                3.6516989226262258,
            ],
            "208Pb_interr": [
                4.23649913953211,
                4.207819597621247,
                4.255153097685344,
                4.302757955585895,
                4.253729851429072,
            ],
        }
    )
    test_unknown_concentrations.index = ["unknown"] * test_unknown_concentrations.shape[
        0
    ]
    test_unknown_concentrations.index.name = "sample"

    test_SRM_concentrations = pd.DataFrame(
        {
            "timestamp": [
                "2021-03-01T22:10:47.999000000",
                "2021-03-01T22:12:05.000000000",
                "2021-03-02T02:41:28.000000000",
                "2021-03-02T02:42:45.999000000",
            ],
            "Spot": ["ATHO-G_1", "ATHO-G_2", "BHVO-2G_4", "BHVO-2G_5"],
            "7Li": [
                25.81744018515236,
                25.854562240659416,
                4.024593214720714,
                4.2809620800304256,
            ],
            "24Mg": [
                579.1203858200361,
                576.2972227801915,
                44366.476235376525,
                44295.92560565091,
            ],
            "27Al": [
                66573.81985653027,
                66174.11201124413,
                71763.68995988581,
                71581.55172053768,
            ],
            "29Si": [
                356151.33226708096,
                352671.9080687457,
                244197.30315501872,
                246300.52034322754,
            ],
            "43Ca": [
                12149.799885648943,
                12149.799885648943,
                81475.12864493996,
                81475.12864493996,
            ],
            "48Ti": [
                1546.7969166023204,
                1555.1089547084653,
                16829.215429171712,
                16758.749560341766,
            ],
            "57Fe": [
                24368.052774644304,
                24138.049743047202,
                88638.22463388307,
                88825.16417683096,
            ],
            "88Sr": [
                95.93421857679373,
                95.78763845519435,
                402.5613712527907,
                406.027654866846,
            ],
            "138Ba": [
                548.0963736274922,
                550.408641930586,
                131.94863924234724,
                132.14344009000718,
            ],
            "139La": [
                57.00444892071998,
                57.43836192819776,
                15.284105034707935,
                15.19366825548261,
            ],
            "140Ce": [
                121.08702224217807,
                122.91460266142933,
                37.468464199826016,
                37.35400405818966,
            ],
            "153Eu": [
                2.7926972025277355,
                2.8218407563769374,
                2.076862251038057,
                2.1162913433696207,
            ],
            "208Pb": [
                5.639772926125088,
                5.713687284506415,
                1.8767279716890939,
                1.9512157976817952,
            ],
            "7Li_exterr": [
                2.951486943862664,
                2.9578175450357205,
                0.4683812589848144,
                0.4962596101382965,
            ],
            "24Mg_exterr": [
                21.497848169408744,
                20.729952231225656,
                1419.605890851989,
                1418.6305522177533,
            ],
            "27Al_exterr": [
                2635.9503360618914,
                2603.7293174408505,
                2574.0526852690286,
                2588.192223835008,
            ],
            "29Si_exterr": [
                10400.216595564447,
                10166.547068755182,
                5790.043250888171,
                5804.5536804239755,
            ],
            "43Ca_exterr": [
                342.9898164757125,
                342.9898164757125,
                1932.2930073874898,
                1932.2930073874898,
            ],
            "48Ti_exterr": [
                116.3332500005449,
                116.77746273475508,
                1233.457530277853,
                1228.8521961530887,
            ],
            "57Fe_exterr": [
                874.594068034649,
                855.8745610317629,
                2801.461177043928,
                2815.4884022280294,
            ],
            "88Sr_exterr": [
                2.6949611315303903,
                2.6630775350080635,
                9.040158510227755,
                9.08408517641062,
            ],
            "138Ba_exterr": [
                15.388547429810092,
                15.200668377102298,
                2.991182476077185,
                2.9583172282671333,
            ],
            "139La_exterr": [
                1.6526603214315219,
                1.6425945478652604,
                0.37140881746898735,
                0.3612599667553517,
            ],
            "140Ce_exterr": [
                3.3817793092788118,
                3.3618990976024192,
                0.8259040920441594,
                0.8369507642388556,
            ],
            "153Eu_exterr": [
                0.09079560833981733,
                0.08985730422730541,
                0.057409119254989156,
                0.06020763647917922,
            ],
            "208Pb_exterr": [
                0.6553840387269524,
                0.5490102600658079,
                0.18886744030142424,
                0.1867453391361122,
            ],
            "7Li_interr": [
                0.5661931933907094,
                0.5777866435756921,
                0.12442529803779183,
                0.12477928030756609,
            ],
            "24Mg_interr": [
                12.89926808146804,
                11.69791451129508,
                528.579149937524,
                531.1723087167296,
            ],
            "27Al_interr": [
                1386.734374788498,
                1346.983920288319,
                886.8542210599865,
                942.9235637450245,
            ],
            "29Si_interr": [
                8397.266858416917,
                8151.093394208212,
                3977.968191814399,
                3960.58856575192,
            ],
            "43Ca_interr": [
                214.4082332761578,
                214.4082332761578,
                714.694110920526,
                714.694110920526,
            ],
            "48Ti_interr": [
                30.300080581348855,
                29.760697032883225,
                167.55240268069542,
                170.91949071129307,
            ],
            "57Fe_interr": [
                522.6805815280773,
                500.040065126809,
                1158.5006094673265,
                1180.4415065186247,
            ],
            "88Sr_interr": [
                1.941556554057554,
                1.8998617731255163,
                4.496182949655127,
                4.466316868487512,
            ],
            "138Ba_interr": [
                11.503265573256229,
                11.211361304591513,
                1.7005386960024176,
                1.6365840680590042,
            ],
            "139La_interr": [
                1.2095755937385948,
                1.187654531789101,
                0.2162728062210723,
                0.2010321117813631,
            ],
            "140Ce_interr": [
                2.5664775221808935,
                2.5110245875207173,
                0.46665659694211664,
                0.48884566385644124,
            ],
            "153Eu_interr": [
                0.0744759936234808,
                0.07294230059636406,
                0.04247500954173352,
                0.04556426666627202,
            ],
            "208Pb_interr": [
                0.3986711603106342,
                0.15389476633181678,
                0.07554836727309974,
                0.04984735608154512,
            ],
        }
    )

    test_SRM_concentrations.index = ["ATHO-G", "ATHO-G", "BHVO-2G", "BHVO-2G"]
    test_SRM_concentrations.index.name = "sample"

    pd.testing.assert_frame_equal(
        test_unknown_concentrations,
        concentrations.unknown_concentrations,
        check_index_type=False,
    )
    pd.testing.assert_frame_equal(
        test_SRM_concentrations,
        concentrations.SRM_concentrations.iloc[[0, 1, 18, 19], :],
        check_index_type=False,
    )
