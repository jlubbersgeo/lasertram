"""
various tests for the package lasertram
"""
import numpy as np
import pandas as pd
import pytest

from lasertram import lasertram as lt
from lasertram.lasertram import LaserCalc, LaserTRAM

###########UNIT TESTS##############
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
