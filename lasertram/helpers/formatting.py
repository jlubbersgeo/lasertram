"""
# TODO: add objects to represent all the various formatting things in other files
- [x] add: srm formatting
- [x] add: lt_ready formatting for upload into lasertram
- [x] add: lt_complete formatting for upload into lasercalc
- [x] add: checks for the correct types in each column for all uploads
- [ ] add: duplicate name checks for lasercalc upload
"""

import pandas as pd
import string
from tabulate import tabulate


def _check_cols(in_cols: list, correct_cols: list) -> None | list:
    """check to see if all elements of an input list are in a reference list
    if not, return a list of the missing elements that are supposed to be in the
    input list but are not.

    Args:
        in_cols (list): the input list to check against the reference
        correct_cols (list): the reference list to use for checking that all elements of the 
        input list are preset

    Returns:
        None | list:
        if there are no missing elements will return None, if there are missing elements will return
        a list of the missing elements
    """
    missing_cols = []
    for col in correct_cols:
        if col in in_cols:
            pass
        else:
            missing_cols.append(col)

    if len(missing_cols) > 0:

        pass
    else:
        missing_cols = None

    return missing_cols

def _check_col_types(in_df, correct_types: list) -> None | list:
    """item by item check of two lists to see if the item types match at each index

    Args:
        in_types (list): item types to be checked
        correct_types (list): reference item types to check against

    Returns:
        None | list: if there are no incorrect types will return none, if there are incorrect types will return 
        a list of the indices where the types do not match
    """

    incorrect_types = []
    
    for i, correct_type in enumerate(correct_types):
        # print(f"checking column {in_df.columns[i]} of type {in_df[in_df.columns[i]].dtype} for type {correct_type}")
        # special datetime check
        if correct_type == 'datetime':
            if pd.api.types.is_datetime64_any_dtype(in_df[in_df.columns[i]]):
                pass
            else:
                incorrect_types.append(i)
        # special string check
        elif correct_type == 'str':
            if pd.api.types.is_string_dtype(in_df[in_df.columns[i]].dropna()):
                pass
            else:
                incorrect_types.append(i)
        # special float check
        elif correct_type == 'float':
            if pd.api.types.is_numeric_dtype(in_df[in_df.columns[i]]):
                pass
            else:
                incorrect_types.append(i)
        # special bool check
        elif correct_type == 'bool':
            if pd.api.types.is_bool_dtype(in_df[in_df.columns[i]]):
                pass
            else:
                incorrect_types.append(i)
        # general type check

        # if they're not datetimes, strings, or the same type
        # append to incorrect types output
        else:
            incorrect_types.append(i)

    if len(incorrect_types) > 0:
        pass
    else:
        incorrect_types = None

    return incorrect_types





def check_srm_format(df: pd.DataFrame) -> None | list:
    """check to make sure uploaded standards database is in the correct format 

    Args:
        df (pd.DataFrame): input dataframe to be used as the standards database

    Returns:
        None | list: _description_
    """

    in_cols = df.columns.to_list()
    correct_cols = [
        "Standard",
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
        "Ag_std",
        "Al_std",
        "As_std",
        "Au_std",
        "B_std",
        "Ba_std",
        "Be_std",
        "Bi_std",
        "Br_std",
        "Ca_std",
        "Cd_std",
        "Ce_std",
        "Cl_std",
        "Co_std",
        "Cr_std",
        "Cs_std",
        "Cu_std",
        "Dy_std",
        "Er_std",
        "Eu_std",
        "F_std",
        "Fe_std",
        "Ga_std",
        "Gd_std",
        "Ge_std",
        "Hf_std",
        "Ho_std",
        "In_std",
        "K_std",
        "La_std",
        "Li_std",
        "Lu_std",
        "Mg_std",
        "Mn_std",
        "Mo_std",
        "Na_std",
        "Nb_std",
        "Nd_std",
        "Ni_std",
        "P_std",
        "Pb_std",
        "Pr_std",
        "Rb_std",
        "Re_std",
        "S_std",
        "Sb_std",
        "Sc_std",
        "Se_std",
        "Si_std",
        "Sm_std",
        "Sn_std",
        "Sr_std",
        "Ta_std",
        "Tb_std",
        "Th_std",
        "Ti_std",
        "Tl_std",
        "Tm_std",
        "U_std",
        "V_std",
        "W_std",
        "Y_std",
        "Yb_std",
        "Zn_std",
        "Zr_std",
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
        "SiO2_std",
        "TiO2_std",
        "Sl2O3_std",
        "FeO_std",
        "MgO_std",
        "MnO_std",
        "CaO_std",
        "Na2O_std",
        "K2O_std",
        "P2O5_std",
    ]

    correct_types = ['str'] + ['float'] * (len(correct_cols) - 1)
    
    col_names_check = _check_cols(in_cols, correct_cols)
    col_types_check = _check_col_types(df, correct_types)

    return col_names_check, col_types_check


def check_lt_input_format(df: pd.DataFrame) -> None | list:
    """check to make sure that data uploaded into lasertram have the correct format


    Args:
        df (pd.DataFrame): input dataframe to be used as the input to lasertram for initial processing.

    Returns:
        None | list: _description_
    """
    in_cols = df.columns.to_list()
    in_cols = in_cols[:3]
    
    # columns and types to check for
    correct_cols = ['SampleLabel','timestamp','Time']
    correct_types = ['str', 'datetime', 'float']

    col_names_check = _check_cols(in_cols, correct_cols)

    col_types_check = _check_col_types(df[in_cols], correct_types)

    return col_names_check, col_types_check


def check_lt_complete_format(df: pd.DataFrame) -> None | list:
    """check to see that the data uploaded for processing using lasercalc have the correct format

    Args:
        df (pd.DataFrame): input dataframe to be used as the input to lasercalc after it has been processed in lasertram

    Returns:
        None | list: _description_
    """

    in_cols = df.columns.to_list()
    # columns and types to check for
    # omit despike and omitted region columns
    # because they have mixed data types
    indices_to_get = [0,1,4,5,6,7,8,9,10]
    in_cols = [in_cols[i] for i in indices_to_get]


    correct_cols = [
            "timestamp",
            "Spot",
            "bkgd_start",
            "bkgd_stop",
            "int_start",
            "int_stop",
            "norm",
            "norm_cps",
        ]
    correct_types = ['datetime', 'str', 'float', 'float', 'float', 'float', 'str', 'float']
    
    col_names_check = _check_cols(in_cols, correct_cols)
    col_types_check = _check_col_types(df[in_cols], correct_types)

    return col_names_check, col_types_check


def check_duplicate_values(df: pd.DataFrame, col: str, print_output: bool = True) -> pd.Series:
    """check a column in a pandas dataframe for duplicate values and return them as a pandas series

    Args:
        df (pd.DataFrame): input dataframe
        col (str): column to check for duplicate values
        print_output (bool, optional): whether or not to print a nicely formatted table of the duplicates. Defaults to True.

    Returns:
        pd.Series: the duplicate values in the specified column and their indices 
    """

    assert isinstance(print_output, bool), "print_output must be boolean"
    assert isinstance(df, pd.core.frame.DataFrame), "df must be a pandas dataframe"
    assert col in df.columns, f"'{col}' is not in the input dataframe columns - please choose a column that exists in the dataframe"

    duplicates = df[col][df[col].duplicated(keep = False).values]
    if duplicates.shape[0] > 0:
        if print_output:

            print("duplicate sample names found:\n")
            print(tabulate(pd.DataFrame(duplicates),headers = 'keys',tablefmt='pipe'))
    else:
        duplicates = None
        if print_output:

            print( f"No duplicate values in column {col} found")


    return duplicates

def rename_duplicate_values(df: pd.DataFrame, col: str, print_output: bool = True) -> pd.DataFrame:
    """rename duplicate values in a specified column from a pandas dataframe. Renaming will append 
    -a, -b, -c ... -z. For use with columns that have values that are solely strings e.g. sample names

    Args:
        df (pd.DataFrame): input dataframe
        col (str): column to check for and rename duplicate values
        print_output (bool, optional): Whether or not to print a nicely formatted table of the duplicates as well as messages that show the results of the renaming. Defaults to True.

    Returns:
        pd.DataFrame: copy of the input dataframe with duplicate values in the specified column renamed. All other values are left alone.
    """

    assert isinstance(print_output, bool), "print_output must be boolean"
    assert isinstance(df, pd.core.frame.DataFrame), "df must be a pandas dataframe"
    assert col in df.columns, f"'{col}' is not in the input dataframe columns - please choose a column that exists in the dataframe"
    
    df_copy = df.copy()

    duplicates = check_duplicate_values(df_copy, col, print_output)

    if duplicates is not None:
        print("Renaming columns:")

        unique_duplicates = duplicates.unique()

        alphabet = list(string.ascii_lowercase)
        for duplicate in unique_duplicates:
            subset = df_copy[col][df_copy[col] == duplicate]
            old_vals = subset.values.copy()
            
            for val, letter in zip(range(len(subset)),alphabet):
                subset.iloc[val] = f"{subset.iloc[val]}-{letter}"

            new_vals = subset.values
            if print_output:

                print(f"Rename {col} indices {subset.index.values} from {old_vals} ---> {new_vals} complete")

            df_copy.loc[subset.index,col] = subset.values
    else:
        print("No duplicates to rename, returning copy of unmodified DataFrame")
    
    return df_copy
    
    
