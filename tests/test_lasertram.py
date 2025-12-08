"""
various tests for the package lasertram
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from matplotlib import pyplot as plt

from lasertram import LaserCalc, LaserTRAM, batch, conversions, plotting, preprocessing, formatting

###########LASERTRAM UNIT TESTS##############
spreadsheet_path = Path(__file__).parent / "2022-05-10_LT_ready.xlsx"
raw_thermo_path = Path(__file__).parent / "raw" / "03152024_IODP_tephras4_1.csv"


pytest.bkgd_interval = (5, 10)
pytest.keep_interval = (25, 40)
pytest.omit_interval = (30, 33)
pytest.int_std = "29Si"

# TODO: general pytest notes
# - [ ] add: upload type checking test for lt_ready, lt_complete, srm_database functions
# - [ ] add: test for checking when analytes are not in SRM data and when there are no published values for an analyte

@pytest.fixture
def load_data():
    data = pd.read_excel(spreadsheet_path).set_index("SampleLabel")
    return data


def test_load_preprocess_thermo():

    assert all(
        k in preprocessing.extract_thermo_data(raw_thermo_path).keys()
        for k in ["timestamp", "file", "sample", "data"]
    )

    df = preprocessing.make_lt_ready_file(raw_thermo_path, quad_type="thermo")

    assert df.shape == (194, 35), "dataframe not the right shape"
    assert df.columns.to_list() == [
        "timestamp",
        "SampleLabel",
        "Time",
        "7Li",
        "29Si",
        "31P",
        "43Ca",
        "45Sc",
        "47Ti",
        "51V",
        "55Mn",
        "65Cu",
        "66Zn",
        "85Rb",
        "88Sr",
        "89Y",
        "90Zr",
        "93Nb",
        "133Cs",
        "137Ba",
        "139La",
        "140Ce",
        "141Pr",
        "146Nd",
        "147Sm",
        "153Eu",
        "157Gd",
        "163Dy",
        "166Er",
        "172Yb",
        "178Hf",
        "181Ta",
        "208Pb",
        "232Th",
        "238U",
    ], "columns not correct"


def test_check_included_data():

    raw_data = preprocessing.load_test_rawdata()

    assert raw_data.shape == (27434, 34), "raw data not the right shape"
    assert raw_data.columns.to_list() == [
        "timestamp",
        "Time",
        "7Li",
        "29Si",
        "31P",
        "43Ca",
        "45Sc",
        "47Ti",
        "51V",
        "55Mn",
        "65Cu",
        "66Zn",
        "85Rb",
        "88Sr",
        "89Y",
        "90Zr",
        "93Nb",
        "133Cs",
        "137Ba",
        "139La",
        "140Ce",
        "141Pr",
        "146Nd",
        "147Sm",
        "153Eu",
        "157Gd",
        "163Dy",
        "166Er",
        "172Yb",
        "178Hf",
        "181Ta",
        "208Pb",
        "232Th",
        "238U",
    ]


def test_make_lt_ready_folder():
    df = preprocessing.make_lt_ready_folder(raw_thermo_path.parent, quad_type="thermo")

    assert df.shape == (973, 35)
    assert df.columns.to_list() == [
        "timestamp",
        "SampleLabel",
        "Time",
        "7Li",
        "29Si",
        "31P",
        "43Ca",
        "45Sc",
        "47Ti",
        "51V",
        "55Mn",
        "65Cu",
        "66Zn",
        "85Rb",
        "88Sr",
        "89Y",
        "90Zr",
        "93Nb",
        "133Cs",
        "137Ba",
        "139La",
        "140Ce",
        "141Pr",
        "146Nd",
        "147Sm",
        "153Eu",
        "157Gd",
        "163Dy",
        "166Er",
        "172Yb",
        "178Hf",
        "181Ta",
        "208Pb",
        "232Th",
        "238U",
    ]


def test_get_intervals():

    intervals = preprocessing.load_test_intervals()

    assert intervals.shape == (168, 2)
    assert intervals.columns.to_list() == ["int_start", "int_stop"]


def test_get_intstd_comps():

    int_std_comps = preprocessing.load_test_int_std_comps()

    assert int_std_comps.shape == (131, 3)
    assert int_std_comps.columns.to_list() == ["Spot", "SiO2", "SiO2_std%"]


def test_get_lt_data(load_data):
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

    spot.assign_int_std(pytest.int_std)

    assert (
        spot.int_std == pytest.int_std
    ), f"the internal standard should be {pytest.int_std}"


def test_assign_intervals(load_data):
    """
    test that the intervals are assigned correctly
    """

    spot = LaserTRAM(name="test")

    samples = load_data.index.unique().dropna().tolist()

    spot.get_data(load_data.loc[samples[0], :])

    spot.assign_int_std(pytest.int_std)

    spot.assign_intervals(
        bkgd=pytest.bkgd_interval, keep=pytest.keep_interval, omit=pytest.omit_interval
    )

    assert spot.bkgd_start == pytest.bkgd_interval[0], "the bkgd_start should be 5"
    assert spot.bkgd_stop == pytest.bkgd_interval[1], "the bkgd_stop should be 10"
    assert spot.int_start == pytest.keep_interval[0], "the int_start should be 20"
    assert spot.int_stop == pytest.keep_interval[1], "the int_stop should be 50"
    assert spot.omit_start == pytest.omit_interval[0], "the omit_start should be 30"
    assert spot.omit_stop == pytest.omit_interval[1], "the omit_stop should be 35"
    assert spot.omitted_region is True, "omittted_region should be True"


def test_plot_raw_data():
    """test that the raw data function plots things correctly"""
    raw_data = preprocessing.load_test_rawdata()

    sample = "GSD-1G_-_1"
    ax = plotting.plot_timeseries_data(raw_data.loc[sample, :])

    lines = ax[0].get_lines()
    for line in lines:
        ydata = line.get_ydata()

        np.testing.assert_array_equal(
            ydata, raw_data.loc[sample, line.get_label()].values
        )

    # now test with specified elements
    ax = plotting.plot_timeseries_data(
        raw_data.loc[sample, :], analytes=["7Li", "29Si"]
    )
    lines = ax[0].get_lines()
    for line in lines:
        ydata = line.get_ydata()

        np.testing.assert_array_equal(
            ydata, raw_data.loc[sample, line.get_label()].values
        )


def test_plot_lt_errors(load_data):
    """
    test that the error bars are being plotted correctly for lasertram"""

    spot = LaserTRAM(name="test")

    samples = load_data.index.unique().dropna().tolist()

    spot.get_data(load_data.loc[samples[0], :])

    spot.assign_int_std(pytest.int_std)

    spot.assign_intervals(
        bkgd=pytest.bkgd_interval, keep=pytest.keep_interval, omit=pytest.omit_interval
    )
    spot.get_bkgd_data()
    spot.subtract_bkgd()
    spot.get_detection_limits()
    spot.normalize_interval()
    spot.despike_data()

    fig, ax = plt.subplots()
    plotting.plot_lasertram_uncertainties(spot, ax=ax)
    patches = ax.patches
    heights = []
    for patch in patches:
        heights.append(patch.get_height())
    heights = np.array(heights)

    np.testing.assert_array_equal(heights, spot.bkgd_subtract_std_err_rel)


def test_get_bkgd_data(load_data):
    """
    test that background signal is being assigned properly
    """
    spot = LaserTRAM(name="test")

    samples = load_data.index.unique().dropna().tolist()

    spot.get_data(load_data.loc[samples[0], :])

    spot.assign_int_std(pytest.int_std)

    spot.assign_intervals(
        bkgd=pytest.bkgd_interval, keep=pytest.keep_interval, omit=pytest.omit_interval
    )
    spot.get_bkgd_data()
    assert np.allclose(
        spot.bkgd_data_median,
        np.array(
            [
                0.00000000e00,
                1.82979534e05,
                4.85094118e03,
                1.50001000e02,
                9.00032401e02,
                0.00000000e00,
                0.00000000e00,
                1.40007840e03,
                1.00000400e02,
                0.00000000e00,
                5.00002000e01,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
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

    spot.assign_int_std(pytest.int_std)

    spot.assign_intervals(
        bkgd=pytest.bkgd_interval, keep=pytest.keep_interval, omit=pytest.omit_interval
    )
    spot.get_bkgd_data()
    spot.subtract_bkgd()

    assert np.allclose(
        spot.bkgd_subtract_data,
        spot.data_matrix[spot.int_start_idx : spot.int_stop_idx, 1:]
        - spot.bkgd_data_median,
    ), "background not subtracted properly"


def test_get_detection_limits(load_data):
    """
    test to make sure detection limits are generated correctly
    """

    spot = LaserTRAM(name="test")

    samples = load_data.index.unique().dropna().tolist()

    spot.get_data(load_data.loc[samples[0], :])

    spot.assign_int_std(pytest.int_std)

    spot.assign_intervals(
        bkgd=pytest.bkgd_interval, keep=pytest.keep_interval, omit=pytest.omit_interval
    )
    spot.get_bkgd_data()
    spot.get_detection_limits()

    assert np.allclose(
        spot.detection_limits,
        np.array(
            [
                2.18530954e02,
                1.95562571e05,
                6.43972601e03,
                6.44267522e02,
                1.80721449e03,
                1.83086988e02,
                3.09052197e02,
                2.48611545e03,
                3.11049711e02,
                2.19579097e02,
                3.98836656e02,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                7.72621221e01,
                1.23098263e02,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ]
        ),
    ), "detection limits not calculated correctly"


def test_despike_data(load_data):
    """
    test to make sure data are despiked properly
    """

    spot = LaserTRAM(name="test")

    samples = load_data.index.unique().dropna().tolist()

    spot.get_data(load_data.loc[samples[0], :])

    spot.assign_int_std(pytest.int_std)

    spot.assign_intervals(
        bkgd=pytest.bkgd_interval, keep=pytest.keep_interval, omit=pytest.omit_interval
    )
    spot.get_bkgd_data()
    spot.subtract_bkgd()
    spot.get_detection_limits()
    spot.normalize_interval()
    spot.despike_data()
    assert np.allclose(
        spot.data_matrix[100],
        np.array(
            [
                3.66069400e01,
                2.44238377e04,
                2.58263964e06,
                6.82860111e04,
                8.11626423e04,
                4.99997996e04,
                4.88980604e05,
                4.24720326e04,
                2.59873579e05,
                8.50289098e03,
                1.18055722e04,
                5.99433846e04,
                1.08065112e05,
                4.78915687e04,
                2.41232548e04,
                4.90962286e04,
                7.72378927e04,
                1.81131139e04,
                6.20536448e04,
                8.48872585e04,
                1.12100416e05,
                1.54094922e04,
                1.55096160e04,
                5.74316333e04,
                1.68112972e04,
                2.69289756e04,
                2.73298442e04,
                2.57264468e04,
                2.20193771e04,
                7.73385120e04,
                6.10487147e04,
                7.19062271e04,
                9.42540181e04,
            ]
        ),
    ), "data not despiked properly"

    assert spot.despiked, "spot.despiked should be True"
    assert (
        spot.despiked_elements == spot.analytes
    ), "The list of despiked elements should be the same as all analytes"


def test_normalize_interval(load_data):
    """
    check that data are being normalized correctly
    """
    spot = LaserTRAM(name="test")

    samples = load_data.index.unique().dropna().tolist()

    spot.get_data(load_data.loc[samples[0], :])

    spot.assign_int_std(pytest.int_std)

    spot.assign_intervals(
        bkgd=pytest.bkgd_interval, keep=pytest.keep_interval, omit=pytest.omit_interval
    )
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
        spot.bkgd_subtract_med,
        np.array(
            [
                0.01043726,
                1.0,
                0.02550994,
                0.03587271,
                0.02063795,
                0.21390534,
                0.01682313,
                0.11149184,
                0.0035571,
                0.00530333,
                0.02326246,
                0.04378744,
                0.02122209,
                0.00968445,
                0.0187391,
                0.02857566,
                0.00648418,
                0.02463792,
                0.03045385,
                0.04082149,
                0.00650887,
                0.00582857,
                0.02133758,
                0.00538868,
                0.00981792,
                0.01022141,
                0.0093663,
                0.0086698,
                0.02729486,
                0.02080482,
                0.02554539,
                0.03572571,
            ]
        ),
    ), "median background and normalized values are incorrect"
    assert np.allclose(
        spot.bkgd_subtract_std_err_rel,
        np.array(
            [
                0.97426461,
                0.0,
                0.95142357,
                0.93041269,
                1.08044579,
                0.93726462,
                1.1810172,
                0.88558184,
                1.89780192,
                2.06505148,
                1.33154415,
                1.25263845,
                1.59623616,
                1.44137923,
                1.30504805,
                1.7108269,
                2.07942066,
                1.49831573,
                1.43300245,
                1.65521111,
                1.70766679,
                1.70396693,
                1.864858,
                2.19004078,
                1.74097691,
                1.96200159,
                1.98419571,
                2.30303923,
                1.77776672,
                2.13811969,
                1.87573994,
                1.84468289,
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
    spot.assign_int_std(pytest.int_std)
    spot.assign_intervals(
        bkgd=pytest.bkgd_interval, keep=pytest.keep_interval, omit=pytest.omit_interval
    )
    spot.get_bkgd_data()
    spot.subtract_bkgd()
    spot.get_detection_limits()
    spot.normalize_interval()
    spot.despike_data()
    spot.make_output_report()
    pd.testing.assert_frame_equal(
        spot.output_report,
        pd.DataFrame(
            {
                "timestamp": {0: pd.Timestamp("2022-05-10 23:08:59")},
                "Spot": {0: "test"},
                "despiked": {
                    0: [
                        "7Li",
                        "29Si",
                        "31P",
                        "43Ca",
                        "45Sc",
                        "47Ti",
                        "51V",
                        "55Mn",
                        "65Cu",
                        "66Zn",
                        "85Rb",
                        "88Sr",
                        "89Y",
                        "90Zr",
                        "93Nb",
                        "133Cs",
                        "137Ba",
                        "139La",
                        "140Ce",
                        "141Pr",
                        "146Nd",
                        "147Sm",
                        "153Eu",
                        "157Gd",
                        "163Dy",
                        "166Er",
                        "172Yb",
                        "178Hf",
                        "181Ta",
                        "208Pb",
                        "232Th",
                        "238U",
                    ]
                },
                "omitted_region": {
                    0: (np.float64(30.019650000000002), np.float64(33.31339))
                },
                "bkgd_start": {0: 5.1363},
                "bkgd_stop": {0: 10.25943},
                "int_start": {0: 25.2623},
                "int_stop": {0: 40.26601},
                "norm": {0: "29Si"},
                "norm_cps": {0: 2819243.2722727386},
                "7Li": {0: 0.010437259034101},
                "29Si": {0: 1.0},
                "31P": {0: 0.02550993641024596},
                "43Ca": {0: 0.03587271241690164},
                "45Sc": {0: 0.020637954074539257},
                "47Ti": {0: 0.2139053355651085},
                "51V": {0: 0.016823127327213718},
                "55Mn": {0: 0.11149183525132685},
                "65Cu": {0: 0.003557100318428378},
                "66Zn": {0: 0.005303326725585216},
                "85Rb": {0: 0.02326245724103047},
                "88Sr": {0: 0.04378744384175842},
                "89Y": {0: 0.021222094964942637},
                "90Zr": {0: 0.00968444848099825},
                "93Nb": {0: 0.018739101968643035},
                "133Cs": {0: 0.028575660636982726},
                "137Ba": {0: 0.0064841776030295445},
                "139La": {0: 0.02463791668966695},
                "140Ce": {0: 0.030453845935454353},
                "141Pr": {0: 0.04082149351465142},
                "146Nd": {0: 0.006508870567873512},
                "147Sm": {0: 0.005828571501401572},
                "153Eu": {0: 0.021337580959465384},
                "157Gd": {0: 0.005388680603868104},
                "163Dy": {0: 0.009817921222589381},
                "166Er": {0: 0.010221411134162577},
                "172Yb": {0: 0.009366301785549105},
                "178Hf": {0: 0.008669797278608683},
                "181Ta": {0: 0.02729485532766075},
                "208Pb": {0: 0.020804823707881108},
                "232Th": {0: 0.025545385313949856},
                "238U": {0: 0.03572571386095365},
                "7Li_se": {0: 0.9742646111044662},
                "29Si_se": {0: 0.0},
                "31P_se": {0: 0.9514235658823135},
                "43Ca_se": {0: 0.930412694285101},
                "45Sc_se": {0: 1.080445788344373},
                "47Ti_se": {0: 0.9372646210445522},
                "51V_se": {0: 1.1810172014404927},
                "55Mn_se": {0: 0.8855818382733421},
                "65Cu_se": {0: 1.8978019212448682},
                "66Zn_se": {0: 2.0650514767520787},
                "85Rb_se": {0: 1.331544147368859},
                "88Sr_se": {0: 1.2526384454267696},
                "89Y_se": {0: 1.5962361610371518},
                "90Zr_se": {0: 1.4413792299360928},
                "93Nb_se": {0: 1.3050480549054508},
                "133Cs_se": {0: 1.710826904827876},
                "137Ba_se": {0: 2.079420660182808},
                "139La_se": {0: 1.4983157335779338},
                "140Ce_se": {0: 1.4330024474613523},
                "141Pr_se": {0: 1.6552111073601967},
                "146Nd_se": {0: 1.7076667903943026},
                "147Sm_se": {0: 1.7039669327734615},
                "153Eu_se": {0: 1.8648580005076827},
                "157Gd_se": {0: 2.190040781302208},
                "163Dy_se": {0: 1.740976911125997},
                "166Er_se": {0: 1.9620015928084438},
                "172Yb_se": {0: 1.9841957089474223},
                "178Hf_se": {0: 2.3030392323015816},
                "181Ta_se": {0: 1.777766722012861},
                "208Pb_se": {0: 2.1381196925657986},
                "232Th_se": {0: 1.875739939915299},
                "238U_se": {0: 1.844682892581869},
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

    spot.assign_int_std(pytest.int_std)

    spot.assign_intervals(
        bkgd=pytest.bkgd_interval, keep=pytest.keep_interval, omit=pytest.omit_interval
    )
    spot.get_bkgd_data()
    spot.subtract_bkgd()
    spot.get_detection_limits()
    spot.normalize_interval()
    spot.make_output_report()

    spot2 = LaserTRAM(name="test")
    batch.process_spot(
        spot2,
        raw_data=load_data.loc[samples[0], :],
        bkgd=pytest.bkgd_interval,
        keep=pytest.keep_interval,
        omit=pytest.omit_interval,
        int_std=pytest.int_std,
        despike=False,
        output_report=True,
    )

    pd.testing.assert_frame_equal(spot.output_report, spot2.output_report)


def test_oxide_to_ppm():
    """
    test that oxides wt% are being converted to elemental ppm
    properly. Test
    """

    elements = ["Si", "Al", "Ca"]
    oxide_vals = np.array([65.0, 15.0, 8.0])

    result = []
    for element, val in zip(elements, oxide_vals):
        result.append(conversions.oxide_to_ppm(val, element))

    result = np.array(result)

    expected = np.array([303833.86315597, 79388.5390063, 57175.66916918])

    np.testing.assert_allclose(
        result,
        expected,
        err_msg="conversions from wt percent oxide to ppm not correct, please check!",
    )

def test_wt_percent_to_oxide():

    elements = ["Si", "Al", "Ca"]

    wt_percent_values = expected = np.array([303833.86315597, 79388.5390063, 57175.66916918]) /1e4
    

    result = []
    for element, val in zip(elements, wt_percent_values):
        result.append(conversions.wt_percent_to_oxide(val, element))

    result = np.array(result)

    expected = np.array([65.0, 15.0, 8.0])

    np.testing.assert_allclose(
        result,
        expected,
        err_msg="conversions from wt percent oxide to ppm not correct, please check!",
    )

def test_supported_oxides():

    
    result = conversions.supported_internal_standard_oxides

    expected = [
        "Al2O3",
        "As2O3",
        "Au2O",
        "B2O3",
        "BaO",
        "BeO",
        "CO2",
        "CaO",
        "Ce2O3",
        "CoO",
        "Cr2O3",
        "Cs2O",
        "CuO",
        "Dy2O3",
        "Er2O3",
        "Eu2O3",
        "FeOT",
        "Ga2O3",
        "Gd2O3",
        "GeO2",
        "H2O",
        "HfO2",
        "Ho2O3",
        "K2O",
        "La2O3",
        "Li2O",
        "Lu2O3",
        "MgO",
        "MnO",
        "MoO3",
        "Na2O",
        "Nb2O5",
        "Nd2O3",
        "NiO",
        "P2O5",
        "PbO",
        "Pr2O3",
        "Rb2O",
        "SO3",
        "Sb2O3",
        "Sc2O3",
        "SiO2",
        "Sm2O3",
        "SnO2",
        "SrO",
        "Ta2O5",
        "Tb2O3",
        "ThO2",
        "TiO2",
        "Tm2O3",
        "V2O5",
        "WO3",
        "Y2O3",
        "Yb2O3",
        "ZnO",
        "ZrO2",
    ]

    assert sorted(result) == sorted(
        expected
    ), "supported internal standard oxides do not match expected - have you changed this list recently??"


###############LASERCALC UNIT TESTS#########################
SRM_path = Path(__file__).parent / "laicpms_stds_tidy.xlsx"
# SRM_path = r"./tests/laicpms_stds_tidy.xlsx"
intervals_path = Path(__file__).parent / "example_intervals.xlsx"
internal_std_path = Path(__file__).parent / "example_internal_std.xlsx"


@pytest.fixture
def load_SRM_data():
    data = pd.read_excel(SRM_path)
    return data


# new for intervals and concentrations
@pytest.fixture
def load_intervals():
    intervals = pd.read_excel(intervals_path).set_index("Spot")
    return intervals


@pytest.fixture
def load_internal_std_comps():
    internal_std_comps = pd.read_excel(internal_std_path)
    return internal_std_comps


# LT_complete_path = r"./tests/spot_test_timestamp_lasertram_complete.xlsx"
LT_complete_path = Path(__file__).parent / "example_lt_complete.xlsx"


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

def test_duplicate_samples(load_LTcomplete_data):
    """test check for finding duplicate sample names in "Spot" column of lasertram output

    Args:
        load_LTcomplete_data (_type_): _description_
    """
    lt_complete = load_LTcomplete_data

    lt_complete.loc[1, "Spot"] = 'GSD-1G_-_1'

    lt_complete.loc[164,"Spot"] = 'ATHO-G_-_4'

    duplicates = formatting.check_duplicate_values(lt_complete, 'Spot',print_output = True)

    expected = pd.Series(
        {0: 'GSD-1G_-_1', 1: 'GSD-1G_-_1', 164: 'ATHO-G_-_4', 165: 'ATHO-G_-_4'}, name = "Spot"
    )

    pd.testing.assert_series_equal(duplicates,expected)


def test_replace_duplicate_samples(load_LTcomplete_data):
    """test replacing duplicate sample names in "Spot" column from lasertram output

    Args:
        load_LTcomplete_data (_type_): _description_
    """
    lt_complete = load_LTcomplete_data

    lt_complete.loc[1, "Spot"] = 'GSD-1G_-_1'

    lt_complete.loc[164,"Spot"] = 'ATHO-G_-_4'
    result = formatting.rename_duplicate_values(lt_complete, 'Spot', print_output = True)
    expected = pd.Series({0: 'GSD-1G_-_1-a', 1: 'GSD-1G_-_1-b', 164: 'ATHO-G_-_4-a', 165: 'ATHO-G_-_4-b'},name = 'Spot')
    
    pd.testing.assert_series_equal(result.loc[[0,1,164,165],'Spot'], expected)


def test_get_lc_data(load_SRM_data, load_LTcomplete_data):
    """_summary_

    Parameters
    ----------
    load_SRM_data : _type_
        _description_
    load_LTcomplete_data : _type_
        _description_
    """
    # potential calibration standards due to nans in the SRM data
    concentrations = LaserCalc(name="test")
    concentrations.get_SRM_comps(load_SRM_data)
    concentrations.get_data(load_LTcomplete_data)
    assert concentrations.spots == [
        "GSD-1G_-_1",
        "GSD-1G_-_2",
        "GSE-1G_-_1",
        "GSE-1G_-_2",
        "BCR-2G_-_1",
        "BCR-2G_-_2",
        "ATHO-G_-_1",
        "ATHO-G_-_2",
        "NIST-612_-_1",
        "NIST-612_-_2",
        "AT-3214-2_shard1_-_1",
        "AT-3214-2_shard1_-_2",
        "AT-3214-2_shard1_-_3",
        "AT-3214-2_shard2_-_1",
        "AT-3214-2_shard2_-_2",
        "AT-3214-2_shard2_-_3",
        "AT-3214-2_shard3_-_2",
        "AT-3214-2_shard3_-_3",
        "AT-3214-2_shard4_-_1",
        "AT-3214-2_shard4_-_2",
        "AT-3214-2_shard5_-_1",
        "AT-3214-2_shard5_-_2",
        "AT-3214-2_shard5_-_3",
        "GSD-1G_-_3",
        "GSE-1G_-_3",
        "AT-3214-2_shard6_-_1",
        "AT-3214-2_shard6_-_2",
        "AT-3214-2_shard6_-_3",
        "AT-3214-2_shard7_-_1",
        "AT-3214-2_shard7_-_3",
        "AT-3214-2_shard8_-_1",
        "AT-3214-2_shard8_-_3",
        "AT-5846_shard1_-_1",
        "AT-5846_shard1_-_2",
        "AT-5846_shard1_-_3",
        "AT-5846_shard1_-_4",
        "AT-5846_shard2_-_1",
        "AT-5846_shard2_-_2",
        "GSD-1G_-_4",
        "GSE-1G_-_4",
        "AT-5846_shard2_-_3",
        "AT-5846_shard2_-_4",
        "AT-5846_shard2_-_5",
        "AT-5846_shard3_-_1",
        "AT-5846_shard3_-_2",
        "AT-5846_shard3_-_3",
        "AT-5846_shard3_-_4",
        "AT-5846_shard3_-_5",
        "AT-5846_shard4_-_1",
        "AT-5846_shard4_-_2",
        "AT-5846_shard4_-_4",
        "AT-5846_shard4_-_5",
        "AT-5846_shard5_-_1",
        "AT-5846_shard5_-_2",
        "GSD-1G_-_5",
        "GSE-1G_-_5",
        "AT-5846_shard5_-_3",
        "AT-5846_shard5_-_4",
        "AT-5846_shard5_-_5",
        "AT-5846_shard6_-_1",
        "AT-5846_shard6_-_2",
        "AT-5846_shard6_-_3",
        "AT-5846_shard6_-_4",
        "AT-5846_shard6_-_5",
        "AT-5846_shard7_-_1",
        "AT-5846_shard7_-_2",
        "AT-5846_shard7_-_3",
        "AT-5846_shard7_-_4",
        "AT-5846_shard7_-_5",
        "AT-5846_shard8_-_1",
        "AT-5846_shard8_-_2",
        "AT-5846_shard8_-_3",
        "AT-5846_shard8_-_4",
        "AT-5846_shard8_-_5",
        "AT-5846_shard9_-_1",
        "AT-5846_shard9_-_2",
        "AT-5846_shard9_-_3",
        "AT-5846_shard9_-_5",
        "AT-5846_shard10_-_1",
        "AT-5846_shard10_-_2",
        "AT-5846_shard11_-_1",
        "AT-5846_shard11_-_2",
        "AT-5846_shard11_-_3",
        "AT-5846_shard11_-_4",
        "GSD-1G_-_7",
        "GSE-1G_-_7",
        "AT-5844_shard1_-_1",
        "AT-5844_shard1_-_2",
        "AT-5844_shard1_-_3",
        "AT-5844_shard1_-_4",
        "AT-5844_shard1_-_5",
        "AT-5844_shard2_-_1",
        "AT-5844_shard2_-_2",
        "AT-5844_shard2_-_3",
        "AT-5844_shard3_-_1",
        "AT-5844_shard3_-_2",
        "AT-5844_shard3_-_4",
        "AT-5844_shard4_-_1",
        "AT-5844_shard4_-_2",
        "GSD-1G_-_8",
        "GSE-1G_-_8",
        "AT-5844_shard4_-_3",
        "AT-5844_shard4_-_4",
        "AT-5844_shard5_-_1",
        "AT-5844_shard5_-_2",
        "AT-5844_shard5_-_3",
        "AT-5844_shard6_-_2",
        "AT-5844_shard7_-_1",
        "AT-5844_shard7_-_2",
        "AT-5844_shard8_-_1",
        "GSD-1G_-_9",
        "GSE-1G_-_9",
        "AT-5844_shard8_-_3",
        "AT-5844_shard9_-_1",
        "AT-5844_shard9_-_3",
        "AT-5844_shard10_-_1",
        "AT-5844_shard10_-_3",
        "AT-5844_shard10_-_4",
        "AT-5843_shard1_-_1",
        "AT-5843_shard1_-_2",
        "AT-5843_shard1_-_3",
        "AT-5843_shard1_-_4",
        "AT-5843_shard2_-_1",
        "AT-5843_shard2_-_2",
        "GSD-1G_-_10",
        "AT-5843_shard2_-_4",
        "AT-5843_shard2_-_5",
        "AT-5843_shard3_-_1",
        "AT-5843_shard3_-_2",
        "AT-5843_shard3_-_3",
        "AT-5843_shard3_-_4",
        "AT-5843_shard3_-_5",
        "AT-5843_shard4_-_1",
        "AT-5843_shard4_-_2",
        "AT-5843_shard4_-_3",
        "AT-5843_shard5_-_2",
        "AT-5843_shard5_-_3",
        "AT-5843_shard5_-_4",
        "GSD-1G_-_11",
        "GSE-1G_-_11",
        "AT-5843_shard5_-_5",
        "AT-5843_shard6_-_1",
        "AT-5843_shard6_-_2",
        "AT-5843_shard6_-_3",
        "AT-5843_shard7_-_1",
        "AT-5843_shard7_-_2",
        "AT-5843_shard7_-_3",
        "AT-5843_shard7_-_4",
        "AT-5843_shard7_-_5",
        "AT-5843_shard8_-_2",
        "AT-5843_shard8_-_3",
        "AT-5843_shard9_-_1",
        "GSD-1G_-_12",
        "GSE-1G_-_12",
        "AT-5843_shard9_-_3",
        "AT-5843_shard10_-_1",
        "AT-5843_shard10_-_2",
        "AT-5843_shard10_-_3",
        "GSD-1G_-_13",
        "GSD-1G_-_14",
        "GSE-1G_-_13",
        "GSE-1G_-_14",
        "BCR-2G_-_3",
        "BCR-2G_-_4",
        "ATHO-G_-_3",
        "ATHO-G_-_4",
        "NIST-612_-_3",
        "NIST-612_-_4",
    ], "analysis spots not found correctly"
    assert concentrations.potential_calibration_standards == [
        "ATHO-G",
        "BCR-2G",
        "GSD-1G",
        "GSE-1G",
        "NIST-612",
    ], "potential calibration standards not found correctly"
    assert concentrations.samples_nostandards == [
        "unknown"
    ], "unknown analyses not found correctly"
    assert concentrations.elements == [
        "Li",
        "Si",
        "P",
        "Ca",
        "Sc",
        "Ti",
        "V",
        "Mn",
        "Cu",
        "Zn",
        "Rb",
        "Sr",
        "Y",
        "Zr",
        "Nb",
        "Cs",
        "Ba",
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Sm",
        "Eu",
        "Gd",
        "Dy",
        "Er",
        "Yb",
        "Hf",
        "Ta",
        "Pb",
        "Th",
        "U",
    ], "analyte to element conversion not correct"


def test_set_calibration_standard(load_SRM_data, load_LTcomplete_data):
    """
    test whether or not calibration standard data is properly assigned
    """
    concentrations = LaserCalc(name="test")
    concentrations.get_SRM_comps(load_SRM_data)
    concentrations.get_data(load_LTcomplete_data)
    concentrations.set_calibration_standard("GSD-1G")
    test_means = pd.Series(
        {
            "7Li": 0.010514340000220277,
            "29Si": 1.0,
            "31P": 0.025310021780240298,
            "43Ca": 0.0370916922740531,
            "45Sc": 0.020889371447489306,
            "47Ti": 0.21808824754316664,
            "51V": 0.01680311011458345,
            "55Mn": 0.11086931664495428,
            "65Cu": 0.0033192819392361688,
            "66Zn": 0.00538551006241637,
            "85Rb": 0.023780278610644006,
            "88Sr": 0.0458124521713978,
            "89Y": 0.022219077036608106,
            "90Zr": 0.010176602015108708,
            "93Nb": 0.018954572715224556,
            "133Cs": 0.029569786588847926,
            "137Ba": 0.006920249495778854,
            "139La": 0.025776100916924973,
            "140Ce": 0.03170470023214307,
            "141Pr": 0.042701015762156756,
            "146Nd": 0.006765706468540337,
            "147Sm": 0.006111752387191597,
            "153Eu": 0.021745128675798656,
            "157Gd": 0.005784290157112075,
            "163Dy": 0.010391424206551071,
            "166Er": 0.010458064791942844,
            "172Yb": 0.009639013782208355,
            "178Hf": 0.008836468906653999,
            "181Ta": 0.027866368422299996,
            "208Pb": 0.0221058848225871,
            "232Th": 0.026274269275601274,
            "238U": 0.03797459547232614,
        }
    )
    test_ses = pd.Series(
        {
            "7Li": 0.6105826456771032,
            "29Si": 0.0,
            "31P": 0.37609726278010913,
            "43Ca": 0.8181052440219647,
            "45Sc": 0.3997294135301341,
            "47Ti": 0.66784873601474,
            "51V": 0.5534010798860425,
            "55Mn": 0.21344500508247022,
            "65Cu": 0.8648118581479076,
            "66Zn": 0.5815308508092042,
            "85Rb": 0.7383146590763143,
            "88Sr": 0.7558710199245325,
            "89Y": 0.8545119820624132,
            "90Zr": 0.674839265226237,
            "93Nb": 0.5243827009967068,
            "133Cs": 0.5353765848118395,
            "137Ba": 0.8017829253581518,
            "139La": 0.6303908414537106,
            "140Ce": 0.7100211634636429,
            "141Pr": 0.5932059593301126,
            "146Nd": 0.7282897866007346,
            "147Sm": 0.7229524547080699,
            "153Eu": 0.442190296432719,
            "157Gd": 0.9088055996443535,
            "163Dy": 0.5559032996500961,
            "166Er": 0.504096058802396,
            "172Yb": 0.37575150143491476,
            "178Hf": 0.47892875794985806,
            "181Ta": 0.5089185148379584,
            "208Pb": 0.9049237861922058,
            "232Th": 0.8068702643170635,
            "238U": 1.2949064704079674,
        }
    )
    assert concentrations.calibration_std == "GSD-1G"
    pd.testing.assert_series_equal(concentrations.calibration_std_means, test_means)
    pd.testing.assert_series_equal(concentrations.calibration_std_ses, test_ses)


def test_drift_check(load_SRM_data, load_LTcomplete_data):
    """
    test whether or not drift is accounted for properly
    """
    concentrations = LaserCalc(name="test")
    concentrations.get_SRM_comps(load_SRM_data)
    concentrations.get_data(load_LTcomplete_data)
    concentrations.set_calibration_standard("GSD-1G")
    concentrations.drift_check()
    test_ses = pd.Series(
        {
            "7Li": "False",
            "29Si": "False",
            "31P": "False",
            "43Ca": "False",
            "45Sc": "False",
            "47Ti": "False",
            "51V": "True",
            "55Mn": "False",
            "65Cu": "False",
            "66Zn": "False",
            "85Rb": "True",
            "88Sr": "False",
            "89Y": "False",
            "90Zr": "False",
            "93Nb": "True",
            "133Cs": "True",
            "137Ba": "False",
            "139La": "False",
            "140Ce": "True",
            "141Pr": "False",
            "146Nd": "True",
            "147Sm": "False",
            "153Eu": "True",
            "157Gd": "False",
            "163Dy": "False",
            "166Er": "False",
            "172Yb": "False",
            "178Hf": "False",
            "181Ta": "False",
            "208Pb": "False",
            "232Th": "True",
            "238U": "True",
        }
    )
    test_ses.name = "drift_correct"
    (
        pd.testing.assert_series_equal(
            test_ses, concentrations.calibration_std_stats["drift_correct"]
        ),
        "analytes not being drift corrected properly",
    )


def test_get_calibration_std_ratios(load_SRM_data, load_LTcomplete_data):
    """
    test that the concentration ratio between every analyte and the internal
    standard is accurate
    """

    concentrations = LaserCalc(name="test")
    concentrations.get_SRM_comps(load_SRM_data)
    concentrations.get_data(load_LTcomplete_data)
    concentrations.set_calibration_standard("GSD-1G")
    concentrations.drift_check()
    concentrations.get_calibration_std_ratios()
    test_ratios = np.array(
        [
            1.72905263e-04,
            1.00000000e00,
            3.45810526e-03,
            2.06915230e-01,
            2.09094737e-04,
            2.98909254e-02,
            1.76926316e-04,
            8.84631579e-04,
            1.68884211e-04,
            2.17136842e-04,
            1.49985263e-04,
            2.79061053e-04,
            1.68884211e-04,
            1.68884211e-04,
            1.68884211e-04,
            1.28673684e-04,
            2.69410526e-04,
            1.57223158e-04,
            1.66471579e-04,
            1.80947368e-04,
            1.79741053e-04,
            1.92206316e-04,
            1.64863158e-04,
            2.03867368e-04,
            2.05877895e-04,
            1.61244211e-04,
            2.04671579e-04,
            1.56821053e-04,
            1.60842105e-04,
            2.01052632e-04,
            1.64863158e-04,
            1.64863158e-04,
        ]
    )
    assert np.allclose(
        concentrations.calibration_std_conc_ratios, test_ratios, rtol = 1e-3,
    ), "calibration standard concentration ratios are not correct, check again"


def test_set_int_std_concentrations(
    load_SRM_data, load_LTcomplete_data, load_internal_std_comps
):
    """
    test to make sure concentration of the internal standard is being set correctly
    """
    concentrations = LaserCalc(name="test")
    concentrations.get_SRM_comps(load_SRM_data)
    concentrations.get_data(load_LTcomplete_data)
    concentrations.set_calibration_standard("GSD-1G")
    concentrations.drift_check()
    concentrations.get_calibration_std_ratios()
    concentrations.set_int_std_concentrations(
        spots=load_internal_std_comps["Spot"],
        concentrations=load_internal_std_comps["SiO2"],
        uncertainties=load_internal_std_comps["SiO2_std%"],
    )
    # TODO check that all different units work for this, oxide and element (wt percent and ppm)

    assert np.allclose(
        concentrations.data.loc["unknown", "int_std_comp"].values,
        load_internal_std_comps["SiO2"].values,
    ), "internal standard concentrations for unknowns not set properly"
    assert np.allclose(
        concentrations.data.loc["unknown", "int_std_rel_unc"].values,
        load_internal_std_comps["SiO2_std%"].values,
    ), "internal standard concentration uncertainties for unknowns not set properly"

def test_calculate_concentrations_units(load_SRM_data,load_LTcomplete_data,load_internal_std_comps):
    """test to make sure that no matter which internal standard units are used (wt percent oxide, wt percent element, ppm element), the output concentration
    is the same

    Args:
        load_SRM_data (_type_): _description_
        load_LTcomplete_data (_type_): _description_
        load_internal_std_comps (_type_): _description_
    """


    out = []

    for units in ["wt_per_ox","wt_per_el","ppm_el"]:
        concentrations = LaserCalc(name="test")
        concentrations.get_SRM_comps(load_SRM_data)
        concentrations.get_data(load_LTcomplete_data)
        concentrations.set_calibration_standard("GSD-1G")
        concentrations.drift_check()
        concentrations.get_calibration_std_ratios()

        int_std_data = load_internal_std_comps

        if units == 'ppm_el':
            values= conversions.oxide_to_ppm(int_std_data['SiO2'],'Si')

        elif units == 'wt_per_el':
            values = conversions.oxide_to_ppm(int_std_data['SiO2'],'Si') / 1e4

        else:
            values = int_std_data['SiO2']

        concentrations.set_int_std_concentrations(
            spots=int_std_data["Spot"],
            concentrations=values,
            uncertainties=int_std_data["SiO2_std%"],
            units = units
        )
        concentrations.calculate_concentrations()


        out.append(concentrations.unknown_concentrations)

    pd.testing.assert_frame_equal(out[0],out[1])
    pd.testing.assert_frame_equal(out[0],out[2])


def test_calculate_concentrations(
    load_SRM_data, load_LTcomplete_data, load_internal_std_comps
):
    """
    test to make sure concentrations are calculated correctly
    """

    concentrations = LaserCalc(name="test")
    concentrations.get_SRM_comps(load_SRM_data)
    concentrations.get_data(load_LTcomplete_data)
    concentrations.set_calibration_standard("GSD-1G")
    concentrations.drift_check()
    concentrations.get_calibration_std_ratios()
    concentrations.set_int_std_concentrations(
        spots=load_internal_std_comps["Spot"],
        concentrations=load_internal_std_comps["SiO2"],
        uncertainties=load_internal_std_comps["SiO2_std%"],
    )
    concentrations.calculate_concentrations()
    test_unknown_concentrations = pd.DataFrame(
        {
            "sample": {0: "unknown", 1: "unknown", 2: "unknown"},
            "timestamp": {
                0: pd.Timestamp("2022-05-10 23:20:00"),
                1: pd.Timestamp("2022-05-10 23:23:54"),
                2: pd.Timestamp("2022-05-10 23:28:38"),
            },
            "Spot": {
                0: "AT-3214-2_shard1_-_1",
                1: "AT-3214-2_shard2_-_2",
                2: "AT-3214-2_shard4_-_1",
            },
            "7Li": {0: 25.83466730151736, 1: 23.918060502303252, 2: 22.898910852608875},
            "29Si": {0: 300468.3188256246, 1: 296869.0561390077, 2: 295840.6953714029},
            "31P": {0: 1421.3751008530235, 1: 1636.101981644864, 2: 1532.6844041862823},
            "43Ca": {
                0: 25099.494041325357,
                1: 22648.37666159304,
                2: 28832.776841374318,
            },
            "45Sc": {
                0: 24.208874127662884,
                1: 20.326108676842324,
                2: 21.46088134459972,
            },
            "47Ti": {0: 4732.74988863527, 1: 4680.5058504286435, 2: 4576.317766288286},
            "51V": {0: 43.73936324577931, 1: 48.60322423359829, 2: 42.34470309581641},
            "55Mn": {
                0: 1505.9792608077726,
                1: 1617.0033414643483,
                2: 1279.8000674172301,
            },
            "65Cu": {0: 12.644458185332347, 1: 5.982796702104342, 2: 73.59571667001988},
            "66Zn": {0: 56.7048398895819, 1: 124.96016798388698, 2: 56.46187108380628},
            "85Rb": {0: 79.84142247722356, 1: 72.6255196165883, 2: 65.61281239603193},
            "88Sr": {
                0: 237.64030939645255,
                1: 219.53128092884606,
                2: 338.14561833843857,
            },
            "89Y": {0: 50.071248735557745, 1: 43.670889832171106, 2: 51.94861169093933},
            "90Zr": {
                0: 276.3598102099743,
                1: 246.43134646221145,
                2: 272.00146173823344,
            },
            "93Nb": {0: 5.408509197734714, 1: 5.525396860821138, 2: 5.22902288360957},
            "133Cs": {0: 4.4101696764603, 1: 4.015918052902527, 2: 4.068127725556567},
            "137Ba": {0: 787.3578800883071, 1: 790.5669866567802, 2: 791.4237350674512},
            "139La": {
                0: 21.635706627513255,
                1: 20.726410440578963,
                2: 23.768057076480215,
            },
            "140Ce": {
                0: 46.269236366419335,
                1: 47.68409912316797,
                2: 43.81208805810481,
            },
            "141Pr": {0: 6.743673359933416, 1: 6.951852128842739, 2: 6.588177637867477},
            "146Nd": {0: 33.18163680025176, 1: 34.05870726570113, 2: 33.01610824946759},
            "147Sm": {0: 8.367230094612692, 1: 7.262884090611597, 2: 7.869364509578179},
            "153Eu": {
                0: 1.8940535609036937,
                1: 1.9087940064353053,
                2: 2.0813461448657553,
            },
            "157Gd": {0: 9.218809116864849, 1: 9.543399686837528, 2: 8.32700367389847},
            "163Dy": {0: 8.49714954834481, 1: 7.37994516103108, 2: 8.062732138972914},
            "166Er": {
                0: 5.593944248493936,
                1: 4.850702045502252,
                2: 4.9775984135422275,
            },
            "172Yb": {0: 5.174813678097828, 1: 4.56357817921105, 2: 5.469729611360127},
            "178Hf": {0: 7.462626444510423, 1: 6.497386801216015, 2: 8.29006322633164},
            "181Ta": {
                0: 0.3764927934496587,
                1: 0.36281882414153954,
                2: 0.37921256753111593,
            },
            "208Pb": {
                0: 13.72550793674585,
                1: 17.071974536655834,
                2: 26.187512800220095,
            },
            "232Th": {0: 6.723318480334518, 1: 6.333449325128908, 2: 6.545985789687301},
            "238U": {0: 2.6739828189203343, 1: 3.244106257133943, 2: 2.395583056221252},
            "7Li_exterr": {
                0: 3.746690556360006,
                1: 3.4044568049955157,
                2: 3.419132330086863,
            },
            "29Si_exterr": {
                0: 7061.059613016031,
                1: 6358.337511664369,
                2: 7323.397435305808,
            },
            "31P_exterr": {
                0: 266.0660996201221,
                1: 306.1288489676574,
                2: 288.3887918741892,
            },
            "43Ca_exterr": {
                0: 890.3019876389949,
                1: 697.8173032747455,
                2: 934.0891172020392,
            },
            "45Sc_exterr": {
                0: 1.3701550661981472,
                1: 1.0325753078279531,
                2: 1.1969194135557393,
            },
            "47Ti_exterr": {
                0: 257.2956393385891,
                1: 244.60114393055707,
                2: 289.76957547009584,
            },
            "51V_exterr": {
                0: 3.1995411310923414,
                1: 3.747898237902943,
                2: 2.4296515339846247,
            },
            "55Mn_exterr": {
                0: 142.51590098502598,
                1: 156.84241915882797,
                2: 128.5474618014899,
            },
            "65Cu_exterr": {
                0: 1.4372234255853766,
                1: 0.6898042855723011,
                2: 12.563694996492702,
            },
            "66Zn_exterr": {
                0: 3.1559717769109317,
                1: 6.137024115949706,
                2: 4.545971320250709,
            },
            "85Rb_exterr": {
                0: 2.861506670590474,
                1: 2.1673168640685834,
                2: 3.580654795729462,
            },
            "88Sr_exterr": {
                0: 10.026184320034725,
                1: 5.1459921536137445,
                2: 12.60555855984909,
            },
            "89Y_exterr": {
                0: 2.6955632373129133,
                1: 2.297741137824117,
                2: 3.5834133825285437,
            },
            "90Zr_exterr": {
                0: 14.599976063966599,
                1: 12.638039756585423,
                2: 18.941670886441518,
            },
            "93Nb_exterr": {
                0: 0.43060120260574236,
                1: 0.45641766993158317,
                2: 0.5799714581511027,
            },
            "133Cs_exterr": {
                0: 0.32284388055490015,
                1: 0.2806951645233029,
                2: 0.38017330898016133,
            },
            "137Ba_exterr": {
                0: 21.125133893905943,
                1: 20.034138832767766,
                2: 37.77781805152602,
            },
            "139La_exterr": {
                0: 0.5482688386277118,
                1: 0.5099509218283487,
                2: 1.2794450618257014,
            },
            "140Ce_exterr": {
                0: 1.3863420330062375,
                1: 1.2973740940711758,
                2: 2.4884248737272983,
            },
            "141Pr_exterr": {
                0: 0.24920998686175078,
                1: 0.23686495671396804,
                2: 0.4633869769543705,
            },
            "146Nd_exterr": {
                0: 1.129532394241278,
                1: 1.2390478522291153,
                2: 1.9820093502265685,
            },
            "147Sm_exterr": {
                0: 0.47435935285628267,
                1: 0.5029018266829225,
                2: 0.8765505832775088,
            },
            "153Eu_exterr": {
                0: 0.130052266984548,
                1: 0.1788383903410621,
                2: 0.19168398089800237,
            },
            "157Gd_exterr": {
                0: 0.4357194349720907,
                1: 0.4936984158394564,
                2: 0.9299565578075277,
            },
            "163Dy_exterr": {
                0: 0.4322460515029835,
                1: 0.40638975079932266,
                2: 0.6302386175454217,
            },
            "166Er_exterr": {
                0: 0.288433405335399,
                1: 0.29360115340903575,
                2: 0.4155726365816195,
            },
            "172Yb_exterr": {
                0: 0.22725980991277314,
                1: 0.3820566856802195,
                2: 0.6482127500580166,
            },
            "178Hf_exterr": {
                0: 0.5078356054143585,
                1: 0.4172465735354423,
                2: 0.704785202709277,
            },
            "181Ta_exterr": {
                0: 0.048110041620031574,
                1: 0.062408119676709,
                2: 0.06065476478327006,
            },
            "208Pb_exterr": {
                0: 0.7490711246165848,
                1: 0.9081482281180053,
                2: 1.6890265146203731,
            },
            "232Th_exterr": {
                0: 0.40930404796875053,
                1: 0.4103324307664915,
                2: 0.6034921726657597,
            },
            "238U_exterr": {
                0: 0.16531990027572227,
                1: 0.20296775789543933,
                2: 0.22354624918415952,
            },
            "7Li_interr": {
                0: 0.9444113188548456,
                1: 0.5680664612109076,
                2: 1.167308761662562,
            },
            "29Si_interr": {
                0: 3004.683188256246,
                1: 755.2205734563631,
                2: 3748.294616399762,
            },
            "31P_interr": {
                0: 20.119924578317118,
                1: 21.345759526459894,
                2: 36.414737030625496,
            },
            "43Ca_interr": {
                0: 727.0865748396102,
                1: 521.5436165949247,
                2: 723.9960572268593,
            },
            "45Sc_interr": {
                0: 0.9369261419387447,
                1: 0.6013457002961587,
                2: 0.8044576738859164,
            },
            "47Ti_interr": {
                0: 93.2370675360729,
                1: 59.87033359682083,
                2: 173.7737431433709,
            },
            "51V_interr": {
                0: 2.419028017135937,
                1: 2.9379979494750663,
                2: 1.3390452030664584,
            },
            "55Mn_interr": {
                0: 32.4706635705944,
                1: 48.98157896517682,
                2: 51.1645317684098,
            },
            "65Cu_interr": {
                0: 1.2910895237696989,
                1: 0.6217481064724331,
                2: 12.014146461399141,
            },
            "66Zn_interr": {
                0: 2.195973640220475,
                1: 3.565422037143328,
                2: 3.946128818978103,
            },
            "85Rb_interr": {
                0: 2.4522708793699337,
                1: 1.7023458989720235,
                2: 3.369348345826858,
            },
            "88Sr_interr": {
                0: 9.055872858668758,
                1: 3.2680385061997006,
                2: 11.018671568927077,
            },
            "89Y_interr": {
                0: 1.006987520679556,
                1: 0.7237107345324737,
                2: 2.4720840340058663,
            },
            "90Zr_interr": {
                0: 4.7648060306809406,
                1: 2.877748892997167,
                2: 13.201923932102424,
            },
            "93Nb_interr": {
                0: 0.17192415761124952,
                1: 0.21365393055255255,
                2: 0.43666971722209363,
            },
            "133Cs_interr": {
                0: 0.15445084084289062,
                1: 0.11020174740292064,
                2: 0.2759388707265419,
            },
            "137Ba_interr": {
                0: 12.960945018911605,
                1: 10.991289497632488,
                2: 33.85257513837871,
            },
            "139La_interr": {
                0: 0.3817808677598914,
                1: 0.34343795317400455,
                2: 1.2042059528043174,
            },
            "140Ce_interr": {
                0: 1.112648779775906,
                1: 0.9781379015440885,
                2: 2.361994117496722,
            },
            "141Pr_interr": {
                0: 0.17135946655898776,
                1: 0.14598211046821571,
                2: 0.4283438364377244,
            },
            "146Nd_interr": {
                0: 0.9429269342019645,
                1: 1.0619758137475164,
                2: 1.8829438884062888,
            },
            "147Sm_interr": {
                0: 0.4489153965201001,
                1: 0.48498491111443637,
                2: 0.8646164449908144,
            },
            "153Eu_interr": {
                0: 0.08698233951325585,
                1: 0.14996480519096567,
                2: 0.15954623831668893,
            },
            "157Gd_interr": {
                0: 0.40294902846536484,
                1: 0.46290947400379007,
                2: 0.9178212650454481,
            },
            "163Dy_interr": {
                0: 0.40450474830303695,
                1: 0.38424310904267395,
                2: 0.6134337556330752,
            },
            "166Er_interr": {
                0: 0.27019273077634526,
                1: 0.28024946393967726,
                2: 0.4057494569998901,
            },
            "172Yb_interr": {
                0: 0.20738265822754995,
                1: 0.37315983330858205,
                2: 0.640724214214251,
            },
            "178Hf_interr": {
                0: 0.3143967561985685,
                1: 0.2313582385014198,
                2: 0.5481283328187775,
            },
            "181Ta_interr": {
                0: 0.029411813426301503,
                1: 0.050483970954845296,
                2: 0.04699426457007916,
            },
            "208Pb_interr": {
                0: 0.4659225817903032,
                1: 0.540836031154741,
                2: 1.265099609337153,
            },
            "232Th_interr": {
                0: 0.22303786431837658,
                1: 0.25269101308122194,
                2: 0.5025437594406454,
            },
            "238U_interr": {
                0: 0.09327209301798912,
                1: 0.11735971454467317,
                2: 0.18713494214486212,
            },
        }
    )

    test_SRM_concentrations = pd.DataFrame(
        {
            "timestamp": {
                "ATHO-G": pd.Timestamp("2022-05-10 23:16:07"),
                "BCR-2G": pd.Timestamp("2022-05-10 23:14:11"),
                "GSE-1G": pd.Timestamp("2022-05-10 23:12:17"),
            },
            "Spot": {
                "ATHO-G": "ATHO-G_-_1",
                "BCR-2G": "BCR-2G_-_1",
                "GSE-1G": "GSE-1G_-_1",
            },
            "7Li": {
                "ATHO-G": 28.768884913887472,
                "BCR-2G": 8.93044517930292,
                "GSE-1G": 425.1463454584108,
            },
            "29Si": {
                "ATHO-G": 353403.14136125654,
                "BCR-2G": 254300.67314884067,
                "GSE-1G": 251028.42183994013,
            },
            "31P": {
                "ATHO-G": 90.5138207140824,
                "BCR-2G": 1362.5115588875708,
                "GSE-1G": 38.741230678889735,
            },
            "43Ca": {
                "ATHO-G": 12314.232847586807,
                "BCR-2G": 50195.18211362366,
                "GSE-1G": 51837.93575615514,
            },
            "45Sc": {
                "ATHO-G": 13.40091718112487,
                "BCR-2G": 36.0287545030922,
                "GSE-1G": 481.0944179115579,
            },
            "47Ti": {
                "ATHO-G": 1487.1275901878787,
                "BCR-2G": 12918.136938672484,
                "GSE-1G": 420.88703728273146,
            },
            "51V": {
                "ATHO-G": 3.67275536338785,
                "BCR-2G": 445.6477405836571,
                "GSE-1G": 445.0274485739689,
            },
            "55Mn": {
                "ATHO-G": 852.9665571969798,
                "BCR-2G": 1617.1683982500567,
                "GSE-1G": 607.2705548074669,
            },
            "65Cu": {
                "ATHO-G": 19.94748681604397,
                "BCR-2G": 20.70791849953658,
                "GSE-1G": 405.30390739772366,
            },
            "66Zn": {
                "ATHO-G": 147.95036905087176,
                "BCR-2G": 155.09947014221015,
                "GSE-1G": 478.09954740623766,
            },
            "85Rb": {
                "ATHO-G": 62.95251128875639,
                "BCR-2G": 46.75983618383034,
                "GSE-1G": 367.24230816666363,
            },
            "88Sr": {
                "ATHO-G": 91.34711316762585,
                "BCR-2G": 339.45977489177625,
                "GSE-1G": 448.0706851877401,
            },
            "89Y": {
                "ATHO-G": 100.09074660479311,
                "BCR-2G": 34.57154420366904,
                "GSE-1G": 426.45034092551947,
            },
            "90Zr": {
                "ATHO-G": 527.1370594354174,
                "BCR-2G": 180.45038267036054,
                "GSE-1G": 401.96893763534064,
            },
            "93Nb": {
                "ATHO-G": 58.88236322079589,
                "BCR-2G": 12.01723667292893,
                "GSE-1G": 428.7729779935295,
            },
            "133Cs": {
                "ATHO-G": 0.935433084641077,
                "BCR-2G": 1.1545306854016135,
                "GSE-1G": 315.151930434475,
            },
            "137Ba": {
                "ATHO-G": 526.3706435247419,
                "BCR-2G": 657.0181651712675,
                "GSE-1G": 406.1303491087627,
            },
            "139La": {
                "ATHO-G": 56.97793780511031,
                "BCR-2G": 24.917418402196525,
                "GSE-1G": 391.94856989741885,
            },
            "140Ce": {
                "ATHO-G": 124.5443748695831,
                "BCR-2G": 52.677882888380935,
                "GSE-1G": 422.1221633819843,
            },
            "141Pr": {
                "ATHO-G": 14.625403615824704,
                "BCR-2G": 6.364232022491421,
                "GSE-1G": 458.5375880248435,
            },
            "146Nd": {
                "ATHO-G": 64.42828750162171,
                "BCR-2G": 29.060065737892124,
                "GSE-1G": 458.29928213565273,
            },
            "147Sm": {
                "ATHO-G": 14.29985935201175,
                "BCR-2G": 6.77264433258067,
                "GSE-1G": 489.5535583954849,
            },
            "153Eu": {
                "ATHO-G": 2.7039917033718504,
                "BCR-2G": 1.9969259239569765,
                "GSE-1G": 420.16971325887806,
            },
            "157Gd": {
                "ATHO-G": 17.1064685478974,
                "BCR-2G": 7.492722619083715,
                "GSE-1G": 532.3865049456698,
            },
            "163Dy": {
                "ATHO-G": 17.038258787642523,
                "BCR-2G": 6.4420347971498,
                "GSE-1G": 523.679861401997,
            },
            "166Er": {
                "ATHO-G": 11.128831792428361,
                "BCR-2G": 3.7723208670236037,
                "GSE-1G": 468.70898178570565,
            },
            "172Yb": {
                "ATHO-G": 11.216787179349128,
                "BCR-2G": 3.292952484516414,
                "GSE-1G": 521.8250988466178,
            },
            "178Hf": {
                "ATHO-G": 14.263375846312426,
                "BCR-2G": 4.740101033249867,
                "GSE-1G": 396.72613241290765,
            },
            "181Ta": {
                "ATHO-G": 3.606627921019415,
                "BCR-2G": 0.7470490392547179,
                "GSE-1G": 401.6728562420093,
            },
            "208Pb": {
                "ATHO-G": 5.883343836373958,
                "BCR-2G": 10.67022569729979,
                "GSE-1G": 388.7885182337956,
            },
            "232Th": {
                "ATHO-G": 7.801102569693415,
                "BCR-2G": 5.917926043680966,
                "GSE-1G": 337.80201034251223,
            },
            "238U": {
                "ATHO-G": 2.341115376530018,
                "BCR-2G": 1.7037244957678466,
                "GSE-1G": 430.57282360029575,
            },
            "7Li_exterr": {
                "ATHO-G": 4.061209153889923,
                "BCR-2G": 1.2719416169375029,
                "GSE-1G": 60.928533901103336,
            },
            "29Si_exterr": {
                "ATHO-G": 8197.066432270218,
                "BCR-2G": 5722.184498613356,
                "GSE-1G": 8812.880979503669,
            },
            "31P_exterr": {
                "ATHO-G": 17.076612641272447,
                "BCR-2G": 254.73425031041637,
                "GSE-1G": 7.500389691926532,
            },
            "43Ca_exterr": {
                "ATHO-G": 338.9794656030722,
                "BCR-2G": 1325.7375006802085,
                "GSE-1G": 1918.7107465966951,
            },
            "45Sc_exterr": {
                "ATHO-G": 0.5997430271294923,
                "BCR-2G": 1.564402322704271,
                "GSE-1G": 24.479483450972626,
            },
            "47Ti_exterr": {
                "ATHO-G": 78.38016570613868,
                "BCR-2G": 678.3289227076508,
                "GSE-1G": 25.139982306372662,
            },
            "51V_exterr": {
                "ATHO-G": 0.20604657505775073,
                "BCR-2G": 22.459111494286066,
                "GSE-1G": 25.417830867447837,
            },
            "55Mn_exterr": {
                "ATHO-G": 79.30859413525295,
                "BCR-2G": 149.98538207887424,
                "GSE-1G": 58.65113000711994,
            },
            "65Cu_exterr": {
                "ATHO-G": 1.1729209397826947,
                "BCR-2G": 1.1374500840223054,
                "GSE-1G": 23.721088416768907,
            },
            "66Zn_exterr": {
                "ATHO-G": 6.543605883380328,
                "BCR-2G": 6.657775018733137,
                "GSE-1G": 24.009982080029534,
            },
            "85Rb_exterr": {
                "ATHO-G": 1.73403766996681,
                "BCR-2G": 1.2922727662526472,
                "GSE-1G": 14.221214810091542,
            },
            "88Sr_exterr": {
                "ATHO-G": 2.208202409779133,
                "BCR-2G": 8.086345731794701,
                "GSE-1G": 16.163232866407522,
            },
            "89Y_exterr": {
                "ATHO-G": 5.262923944185317,
                "BCR-2G": 1.816319628606558,
                "GSE-1G": 25.205274422239714,
            },
            "90Zr_exterr": {
                "ATHO-G": 27.52543360533757,
                "BCR-2G": 9.381041561423181,
                "GSE-1G": 23.69045871745045,
            },
            "93Nb_exterr": {
                "ATHO-G": 4.442186758411581,
                "BCR-2G": 0.911304044656405,
                "GSE-1G": 34.25123828322263,
            },
            "133Cs_exterr": {
                "ATHO-G": 0.06677605402596191,
                "BCR-2G": 0.08305991655862681,
                "GSE-1G": 22.749572255147147,
            },
            "137Ba_exterr": {
                "ATHO-G": 14.247811873023075,
                "BCR-2G": 17.729796673771332,
                "GSE-1G": 15.392696835663328,
            },
            "139La_exterr": {
                "ATHO-G": 1.358802002235542,
                "BCR-2G": 0.5986173854301686,
                "GSE-1G": 14.047973546867828,
            },
            "140Ce_exterr": {
                "ATHO-G": 3.52459603122835,
                "BCR-2G": 1.5210447520356156,
                "GSE-1G": 16.445466212944893,
            },
            "141Pr_exterr": {
                "ATHO-G": 0.45353222738654114,
                "BCR-2G": 0.20393235589189235,
                "GSE-1G": 18.768012322983697,
            },
            "146Nd_exterr": {
                "ATHO-G": 1.885847311498854,
                "BCR-2G": 0.9312988980981447,
                "GSE-1G": 18.11396532636725,
            },
            "147Sm_exterr": {
                "ATHO-G": 0.4763475997846337,
                "BCR-2G": 0.2829049777087526,
                "GSE-1G": 17.66215109837514,
            },
            "153Eu_exterr": {
                "ATHO-G": 0.16127043704503938,
                "BCR-2G": 0.11830645908585853,
                "GSE-1G": 25.45312708056567,
            },
            "157Gd_exterr": {
                "ATHO-G": 0.5750804202706984,
                "BCR-2G": 0.2953548088129469,
                "GSE-1G": 19.399432450575144,
            },
            "163Dy_exterr": {
                "ATHO-G": 0.481475391791309,
                "BCR-2G": 0.2091516022542951,
                "GSE-1G": 18.7358883826142,
            },
            "166Er_exterr": {
                "ATHO-G": 0.3091683027267639,
                "BCR-2G": 0.12055058178613962,
                "GSE-1G": 16.818325556582188,
            },
            "172Yb_exterr": {
                "ATHO-G": 0.31697068189410404,
                "BCR-2G": 0.12601922214244402,
                "GSE-1G": 18.789078843633817,
            },
            "178Hf_exterr": {
                "ATHO-G": 0.8111003429080107,
                "BCR-2G": 0.2951100472024632,
                "GSE-1G": 24.48152858092754,
            },
            "181Ta_exterr": {
                "ATHO-G": 0.3712299915608881,
                "BCR-2G": 0.08210814196911717,
                "GSE-1G": 42.42203579049489,
            },
            "208Pb_exterr": {
                "ATHO-G": 0.2921855681616748,
                "BCR-2G": 0.5125645966199847,
                "GSE-1G": 20.623468361477382,
            },
            "232Th_exterr": {
                "ATHO-G": 0.453840559375883,
                "BCR-2G": 0.3407891313374696,
                "GSE-1G": 21.166985298634533,
            },
            "238U_exterr": {
                "ATHO-G": 0.1363318712373106,
                "BCR-2G": 0.1028261655024203,
                "GSE-1G": 26.287139035257212,
            },
            "7Li_interr": {
                "ATHO-G": 0.43812862695942767,
                "BCR-2G": 0.21682783244228468,
                "GSE-1G": 12.33796329528715,
            },
            "29Si_interr": {
                "ATHO-G": 3272.2513089005233,
                "BCR-2G": 1869.8578908002992,
                "GSE-1G": 7011.967090501121,
            },
            "31P_interr": {
                "ATHO-G": 2.4859340390141753,
                "BCR-2G": 14.57046319136426,
                "GSE-1G": 1.9914651171438318,
            },
            "43Ca_interr": {
                "ATHO-G": 226.63883047377755,
                "BCR-2G": 837.7404511597699,
                "GSE-1G": 1598.5756674111283,
            },
            "45Sc_interr": {
                "ATHO-G": 0.23113786005572054,
                "BCR-2G": 0.48332100873498257,
                "GSE-1G": 14.301143524343253,
            },
            "47Ti_interr": {
                "ATHO-G": 21.573862930699224,
                "BCR-2G": 177.98916039845187,
                "GSE-1G": 13.31191874159885,
            },
            "51V_interr": {
                "ATHO-G": 0.10740017075060583,
                "BCR-2G": 7.011951923410519,
                "GSE-1G": 13.859593182290983,
            },
            "55Mn_interr": {
                "ATHO-G": 10.606899919107574,
                "BCR-2G": 17.050894490276548,
                "GSE-1G": 17.573147613678117,
            },
            "65Cu_interr": {
                "ATHO-G": 0.6192685797652988,
                "BCR-2G": 0.4737597163895182,
                "GSE-1G": 12.371189112895982,
            },
            "66Zn_interr": {
                "ATHO-G": 2.8004417542092788,
                "BCR-2G": 2.426483415679609,
                "GSE-1G": 14.5340676359435,
            },
            "85Rb_interr": {
                "ATHO-G": 1.286460422885337,
                "BCR-2G": 0.9612977096249774,
                "GSE-1G": 12.499429754652445,
            },
            "88Sr_interr": {
                "ATHO-G": 1.4629930903756208,
                "BCR-2G": 5.254335926048232,
                "GSE-1G": 13.979461709620583,
            },
            "89Y_interr": {
                "ATHO-G": 1.6480472086437783,
                "BCR-2G": 0.5644147663640385,
                "GSE-1G": 13.483371617343364,
            },
            "90Zr_interr": {
                "ATHO-G": 8.044599779654426,
                "BCR-2G": 2.6082823146883816,
                "GSE-1G": 12.582040446118338,
            },
            "93Nb_interr": {
                "ATHO-G": 1.1222970177714287,
                "BCR-2G": 0.24700776004667602,
                "GSE-1G": 13.913416269047437,
            },
            "133Cs_interr": {
                "ATHO-G": 0.02903557642787282,
                "BCR-2G": 0.03729251404587188,
                "GSE-1G": 10.349496485159197,
            },
            "137Ba_interr": {
                "ATHO-G": 8.867147270789388,
                "BCR-2G": 10.980412898863303,
                "GSE-1G": 12.762968899474957,
            },
            "139La_interr": {
                "ATHO-G": 0.8788922720240755,
                "BCR-2G": 0.3911077049514785,
                "GSE-1G": 12.104917556615588,
            },
            "140Ce_interr": {
                "ATHO-G": 2.7326218073555384,
                "BCR-2G": 1.1945845660503913,
                "GSE-1G": 14.61253173240063,
            },
            "141Pr_interr": {
                "ATHO-G": 0.22735667700771367,
                "BCR-2G": 0.11147954318770736,
                "GSE-1G": 14.172600926152699,
            },
            "146Nd_interr": {
                "ATHO-G": 1.448580391789084,
                "BCR-2G": 0.755441930503341,
                "GSE-1G": 15.948043117829451,
            },
            "147Sm_interr": {
                "ATHO-G": 0.39785998376025555,
                "BCR-2G": 0.2542521581729974,
                "GSE-1G": 15.216231638644281,
            },
            "153Eu_interr": {
                "ATHO-G": 0.08340651796071483,
                "BCR-2G": 0.06004808737043845,
                "GSE-1G": 13.705867223414787,
            },
            "157Gd_interr": {
                "ATHO-G": 0.48588512625980496,
                "BCR-2G": 0.2628292809833573,
                "GSE-1G": 16.872423409966935,
            },
            "163Dy_interr": {
                "ATHO-G": 0.372138877042072,
                "BCR-2G": 0.17436267228087265,
                "GSE-1G": 16.213145166932765,
            },
            "166Er_interr": {
                "ATHO-G": 0.23506453446355863,
                "BCR-2G": 0.09949160970748813,
                "GSE-1G": 14.536808225392745,
            },
            "172Yb_interr": {
                "ATHO-G": 0.24470192364868007,
                "BCR-2G": 0.11127649225141403,
                "GSE-1G": 16.28430117575458,
            },
            "178Hf_interr": {
                "ATHO-G": 0.2772228116534601,
                "BCR-2G": 0.15139429581048855,
                "GSE-1G": 12.240842419032335,
            },
            "181Ta_interr": {
                "ATHO-G": 0.06922875931137475,
                "BCR-2G": 0.03216715082483274,
                "GSE-1G": 12.236598823488052,
            },
            "208Pb_interr": {
                "ATHO-G": 0.14887327985044205,
                "BCR-2G": 0.2341172053985284,
                "GSE-1G": 12.218679229793542,
            },
            "232Th_interr": {
                "ATHO-G": 0.2177102668915752,
                "BCR-2G": 0.15774035289465613,
                "GSE-1G": 12.276330542019133,
            },
            "238U_interr": {
                "ATHO-G": 0.06561397642781303,
                "BCR-2G": 0.05486178440506828,
                "GSE-1G": 14.420183752480144,
            },
        }
    )

    test_SRM_concentrations.index = ["ATHO-G", "BCR-2G", "GSE-1G"]
    test_SRM_concentrations.index.name = "sample"

    pd.testing.assert_frame_equal(
        test_unknown_concentrations,
        concentrations.unknown_concentrations.iloc[[0, 4, 8], :].reset_index(),
        check_index_type=False,
        check_dtype=False,
    )
    pd.testing.assert_frame_equal(
        test_SRM_concentrations,
        concentrations.SRM_concentrations.iloc[[0, 4, 8], :].head(),
        check_index_type=False,
        check_dtype=False,
    )


def test_SRM_accuracies(load_SRM_data, load_LTcomplete_data, load_internal_std_comps):
    concentrations = LaserCalc(name="test")
    concentrations.get_SRM_comps(load_SRM_data)
    concentrations.get_data(load_LTcomplete_data)
    concentrations.set_calibration_standard("GSD-1G")
    concentrations.drift_check()
    concentrations.get_calibration_std_ratios()
    concentrations.set_int_std_concentrations(
        spots=load_internal_std_comps["Spot"],
        concentrations=load_internal_std_comps["SiO2"],
        uncertainties=load_internal_std_comps["SiO2_std%"],
    )
    concentrations.calculate_concentrations()
    concentrations.get_secondary_standard_accuracies()
    test_SRM_accuracies = pd.DataFrame(
        {
            "timestamp": {
                "ATHO-G": pd.Timestamp("2022-05-10 23:16:07"),
                "BCR-2G": pd.Timestamp("2022-05-10 23:14:11"),
                "GSE-1G": pd.Timestamp("2022-05-10 23:12:17"),
            },
            "Spot": {
                "ATHO-G": "ATHO-G_-_1",
                "BCR-2G": "BCR-2G_-_1",
                "GSE-1G": "GSE-1G_-_1",
            },
            "7Li": {
                "ATHO-G": 100.59050669191424,
                "BCR-2G": 99.22716865892133,
                "GSE-1G": 98.87124312986298,
            },
            "29Si": {"ATHO-G": 100.0, "BCR-2G": 100.0, "GSE-1G": 100.0},
            "31P": {
                "ATHO-G": 82.96858861935648,
                "BCR-2G": 84.38733752288533,
                "GSE-1G": 55.34461525555676,
            },
            "43Ca": {
                "ATHO-G": 101.35338000202034,
                "BCR-2G": 99.48030993396917,
                "GSE-1G": 98.01572933785442,
            },
            "45Sc": {
                "ATHO-G": 191.44167401606956,
                "BCR-2G": 109.17804394876424,
                "GSE-1G": 90.77253168142602,
            },
            "47Ti": {
                "ATHO-G": 97.28147188989806,
                "BCR-2G": 91.6179924728545,
                "GSE-1G": 93.61088891884324,
            },
            "51V": {
                "ATHO-G": 93.93236223498337,
                "BCR-2G": 104.85829190203695,
                "GSE-1G": 101.1426019486293,
            },
            "55Mn": {
                "ATHO-G": 103.90098289176795,
                "BCR-2G": 104.33344504839076,
                "GSE-1G": 102.92721267923167,
            },
            "65Cu": {
                "ATHO-G": 107.24455277442993,
                "BCR-2G": 98.60913571207894,
                "GSE-1G": 106.65892299940097,
            },
            "66Zn": {
                "ATHO-G": 104.92933975239131,
                "BCR-2G": 124.07957611376813,
                "GSE-1G": 103.93468421874732,
            },
            "85Rb": {
                "ATHO-G": 96.40507088630382,
                "BCR-2G": 99.48901315708582,
                "GSE-1G": 103.15795173220889,
            },
            "88Sr": {
                "ATHO-G": 97.07450921107954,
                "BCR-2G": 99.25724412040242,
                "GSE-1G": 100.23952688763761,
            },
            "89Y": {
                "ATHO-G": 105.91613397332604,
                "BCR-2G": 98.77584058191155,
                "GSE-1G": 104.01227827451694,
            },
            "90Zr": {
                "ATHO-G": 102.95645692097996,
                "BCR-2G": 98.07086014693506,
                "GSE-1G": 98.04120430130259,
            },
            "93Nb": {
                "ATHO-G": 94.36276157178828,
                "BCR-2G": 96.13789338343143,
                "GSE-1G": 102.08880428417369,
            },
            "133Cs": {
                "ATHO-G": 86.61417450380341,
                "BCR-2G": 99.52850736220807,
                "GSE-1G": 101.66191304337903,
            },
            "137Ba": {
                "ATHO-G": 96.22863684181753,
                "BCR-2G": 96.19592462244033,
                "GSE-1G": 95.11249393647837,
            },
            "139La": {
                "ATHO-G": 102.47830540487465,
                "BCR-2G": 100.88023644613979,
                "GSE-1G": 99.98688007587215,
            },
            "140Ce": {
                "ATHO-G": 102.92923542940751,
                "BCR-2G": 98.83280091628694,
                "GSE-1G": 101.96187521303968,
            },
            "141Pr": {
                "ATHO-G": 100.17399736866236,
                "BCR-2G": 94.98853764912567,
                "GSE-1G": 99.68208435322684,
            },
            "146Nd": {
                "ATHO-G": 105.79357553632464,
                "BCR-2G": 100.55386068474783,
                "GSE-1G": 101.16981945599399,
            },
            "147Sm": {
                "ATHO-G": 100.70323487332219,
                "BCR-2G": 102.77153767193734,
                "GSE-1G": 100.31835213022232,
            },
            "153Eu": {
                "ATHO-G": 97.97071389028443,
                "BCR-2G": 101.3667981704049,
                "GSE-1G": 102.48041786801903,
            },
            "157Gd": {
                "ATHO-G": 111.80698397318562,
                "BCR-2G": 111.66501667784969,
                "GSE-1G": 103.57714104001359,
            },
            "163Dy": {
                "ATHO-G": 105.17443696075632,
                "BCR-2G": 100.03159622903415,
                "GSE-1G": 99.93890484770934,
            },
            "166Er": {
                "ATHO-G": 108.0469106061006,
                "BCR-2G": 101.95461802766496,
                "GSE-1G": 78.77461878751355,
            },
            "172Yb": {
                "ATHO-G": 106.8265445652298,
                "BCR-2G": 97.13724143116266,
                "GSE-1G": 100.3509805474265,
            },
            "178Hf": {
                "ATHO-G": 104.11223245483524,
                "BCR-2G": 97.93597176136089,
                "GSE-1G": 100.43699554757154,
            },
            "181Ta": {
                "ATHO-G": 92.47763900049782,
                "BCR-2G": 95.77551785316895,
                "GSE-1G": 102.99304006205367,
            },
            "208Pb": {
                "ATHO-G": 103.76267789019325,
                "BCR-2G": 97.00205179363446,
                "GSE-1G": 102.85410535285598,
            },
            "232Th": {
                "ATHO-G": 105.42030499585697,
                "BCR-2G": 100.30383124882992,
                "GSE-1G": 88.89526587960847,
            },
            "238U": {
                "ATHO-G": 98.78123951603452,
                "BCR-2G": 100.81210034129269,
                "GSE-1G": 102.51733895245137,
            },
        }
    )
    test_SRM_accuracies.index = ["ATHO-G", "BCR-2G", "GSE-1G"]
    test_SRM_accuracies.index.name = "sample"

    pd.testing.assert_frame_equal(
        test_SRM_accuracies,
        concentrations.SRM_accuracies.iloc[[0, 4, 8], :].head(),
        check_index_type=False,
        check_dtype=False,
    )
