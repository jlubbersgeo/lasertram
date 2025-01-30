# for preprocessing a folder of data to make a dataframe that
# is ready for LaserTRAM

import re
from pathlib import Path

import numpy as np
import pandas as pd
from dateutil.parser import parse


def extract_agilent_data(file):
    """
    read raw output from an Agilent quadrupole .csv file and
    return a pandas dataframe and metadata ready for processing with LaserTRAM

    Parameters
    ----------
    file : path-like
        path to the csv file for data to be extracted

    Returns
    -------
    dict
        dictionary that contains timestamp, filename, and data
        for preprocessing

    """
    # import data
    # extract sample name
    # extract timestamp
    # extract data and make headers ready for lasertram

    df = pd.read_csv(file, sep="\t", header=None)

    sample = df.iloc[0, 0].split("\\")[-1].split(".")[0].replace("_", "-")

    timestamp = parse(df.iloc[2, 0].split(" ")[7] + " " + df.iloc[2, 0].split(" ")[8])

    data = pd.DataFrame([sub.split(",") for sub in df.iloc[3:-1, 0]])

    header = data.iloc[0, :]
    data = data[1:]
    data.columns = header
    newcols = []
    for s in data.columns.tolist():
        l = re.findall("(\d+|[A-Za-z]+)", s)
        if "Time" in l:
            newcols.append(l[0])
        else:

            newcols.append(l[1] + l[0])
    data.columns = newcols

    return {"timestamp": timestamp, "file": file, "sample": sample, "data": data}


def extract_thermo_data(file):
    """
    read raw output from an ThermoFisher quadrupole .csv file and
    return a pandas dataframe and metadata ready for processing with LaserTRAM

    Parameters
    ----------
    file : path-like
        path to the csv file for data to be extracted

    Returns
    -------
    dict
        dictionary that contains timestamp, filename, and data
        for preprocessing
    """

    # gets the top row in your csv and turns it into a pandas series
    top = pd.read_csv(file, nrows=0)
    # since it is only 1 long it is also the column name
    # extract that as a list
    sample = list(top.columns)

    # turn that list value to a string
    sample = str(sample[0])

    # because its a string it can be split
    # split at : removes the time stamp
    sample = sample.split(":")[0]

    # .strip() removes leading and trailing spaces
    sample = sample.strip()

    # replace middle spaces with _ because spaces are bad
    nospace = sample.replace(" ", "_")

    # get the timestamp by splitting the string by the previously
    # designated sample. Also drops the colon in front of the date
    timestamp = top.columns.tolist()[0].split(sample)[1:][0][1:]

    timestamp = parse(timestamp)

    # import data
    # remove the top rows. Double check that your header is the specified
    # amount of rows to be skipped in 'skiprows' argument
    data = pd.read_csv(file, skiprows=13)
    # drop empty column at the end
    data.drop(data.columns[len(data.columns) - 1], axis=1, inplace=True)

    # remove dwell time row beneath header row
    data = data.dropna()

    return {"timestamp": timestamp, "file": file, "sample": nospace, "data": data}


def make_lt_ready_folder(folder, quad_type):
    """
    Take a folder of csv files from either an Agilent or ThermoFisher
    quadrupole mass spectrometer, and combine their data such that it is
    a pandas.DataFrame ready for processing in LaserTRAM

    Parameters
    ----------
    folder : path-like
        path to the folder where the csv files are. This looks at all csv
        files so make sure ONLY the data are in there.
    quad_type : str
        "agilent" or "thermo"

    Returns
    -------
    pandas.DataFrame
        dataframe ready to be processed using LaserTRAM.
    """

    if isinstance(folder, Path):
        pass
    else:
        folder = Path(folder)
    assert (
        folder.is_dir()
    ), f"{folder} is not a directory, please choose a directory to your data .csv files"
    my_dict = {}
    for i in folder.glob("*.csv"):
        if quad_type == "thermo":
            temp = extract_thermo_data(i)

        elif quad_type == "agilent":
            temp = extract_agilent_data(i)

        my_dict[temp["timestamp"]] = temp

    my_dict = dict(sorted(my_dict.items()))

    outdf = pd.DataFrame()
    for timestamp in my_dict:
        samplelabel = pd.DataFrame(
            np.repeat(
                my_dict[timestamp]["sample"], my_dict[timestamp]["data"].shape[0]
            ),
            columns=["SampleLabel"],
            index=my_dict[timestamp]["data"].index,
        )
        ts = pd.DataFrame(
            np.repeat(
                my_dict[timestamp]["timestamp"], my_dict[timestamp]["data"].shape[0]
            ),
            columns=["timestamp"],
            index=my_dict[timestamp]["data"].index,
        )
        df = pd.concat([ts, samplelabel, my_dict[timestamp]["data"]], axis="columns")

        outdf = pd.concat([outdf, df])
        outdf.index = np.arange(outdf.shape[0], dtype=int)

    return outdf


def make_lt_ready_file(file, quad_type):
    """
    Take an individual csv file from either an Agilent or ThermoFisher
    quadrupole mass spectrometer and convert it to a pandas.DataFrame
    object ready for processing in LaserTRAM

    Parameters
    ----------
    folder : path-like
        path to the csv file.
    quad_type : str
        "agilent" or "thermo"

    Returns
    -------
    pandas.DataFrame
        dataframe ready to be processed using LaserTRAM.
    """

    if isinstance(file, Path):
        pass
    else:
        file = Path(file)

    assert file.name.endswith(".csv"), f"File '{file}' does not have a CSV extension."

    if quad_type == "thermo":
        temp = extract_thermo_data(file)

    elif quad_type == "agilent":
        temp = extract_agilent_data(file)
    else:
        temp = None

    if temp:
        outdf = temp["data"]
        outdf.insert(0, "SampleLabel", temp["sample"])
        outdf.insert(0, "timestamp", temp["timestamp"])

    else:
        raise ValueError("please choose either 'thermo' or 'agilent' for quad_type")

    return outdf


def load_test_rawdata():
    """
    Load in raw data used as examples in the following manuscript:

    Lubbers, J., Kent, A., Russo, C. (2025) "lasertram: a Python
    library for time resolved analysis of laser ablation inductively
    coupled plasma mass spectrometry data "

    """

    current_path = Path(__file__).parent

    lt_ready = pd.read_excel(
        current_path.parents[1]
        / "test_data"
        / "computers_and_geosciences_examples"
        / "2022-05-10_LT_ready.xlsx"
    ).set_index("SampleLabel")

    return lt_ready


def load_test_intervals():
    """
    Load in interval regions used as examples in the following manuscript:

    Lubbers, J., Kent, A., Russo, C. (2025) "lasertram: a Python
    library for time resolved analysis of laser ablation inductively
    coupled plasma mass spectrometry data "

    """

    current_path = Path(__file__).parent

    intervals = pd.read_excel(
        current_path.parents[1]
        / "test_data"
        / "computers_and_geosciences_examples"
        / "example_intervals.xlsx"
    ).set_index("Spot")

    return intervals


def load_test_int_std_comps():
    """
    Load in internal standard comps used as examples in the following manuscript:

    Lubbers, J., Kent, A., Russo, C. (2025) "lasertram: a Python
    library for time resolved analysis of laser ablation inductively
    coupled plasma mass spectrometry data "

    """

    current_path = Path(__file__).parent

    concentrations = pd.read_excel(
        current_path.parents[1]
        / "test_data"
        / "computers_and_geosciences_examples"
        / "example_internal_std.xlsx"
    )

    return concentrations
