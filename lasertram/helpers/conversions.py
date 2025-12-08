"""
conversions module:
For converting wt% oxide to ppm
"""

import pandas as pd
import numpy as np

# TODO: oxide_to_ppm function
# - [ ] add: element_to_oxide and oxide_dict to conversions (take from oxide converter tool) so that we can utilize many more things as internal standards
# - [ ] fix: Make oxide_to_ppm reference oxide_dict as well so we can seamlessly use both functions back and forth
# - [ ] add: in error handling for if element is not found in oxides list
# list of supported oxides - not exhaustive but common geologic ones
# if not in this list they can just use the option to set the internal std concentration as ppm rather than wt
# % oxide
oxide_dict = {
    "Al": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Al2O3",
        "cation_atomic_weight": 26.9815385,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 101.960077,
        "cation_atomic_number": 13,
    },
    "As": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "As2O3",
        "cation_atomic_weight": 74.921595,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 197.84019,
        "cation_atomic_number": 33,
    },
    "Au": {
        "num_oxygens": 1,
        "num_cations": 2,
        "label": "Au2O",
        "cation_atomic_weight": 196.966569,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 409.932138,
        "cation_atomic_number": 79,
    },
    "B": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "B2O3",
        "cation_atomic_weight": 10.81,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 69.617,
        "cation_atomic_number": 5,
    },
    "Ba": {
        "num_oxygens": 1,
        "num_cations": 1,
        "label": "BaO",
        "cation_atomic_weight": 137.327,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 153.326,
        "cation_atomic_number": 56,
    },
    "Be": {
        "num_oxygens": 1,
        "num_cations": 1,
        "label": "BeO",
        "cation_atomic_weight": 9.0121831,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 25.0111831,
        "cation_atomic_number": 4,
    },
    "C": {
        "num_oxygens": 2,
        "num_cations": 1,
        "label": "CO2",
        "cation_atomic_weight": 12.011,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 44.009,
        "cation_atomic_number": 6,
    },
    "Ca": {
        "num_oxygens": 1,
        "num_cations": 1,
        "label": "CaO",
        "cation_atomic_weight": 40.078,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 56.077000000000005,
        "cation_atomic_number": 20,
    },
    "Ce": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Ce2O3",
        "cation_atomic_weight": 140.116,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 328.22900000000004,
        "cation_atomic_number": 58,
    },
    "Co": {
        "num_oxygens": 1,
        "num_cations": 1,
        "label": "CoO",
        "cation_atomic_weight": 58.933194,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 74.932194,
        "cation_atomic_number": 27,
    },
    "Cr": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Cr2O3",
        "cation_atomic_weight": 51.9961,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 151.98919999999998,
        "cation_atomic_number": 24,
    },
    "Cs": {
        "num_oxygens": 1,
        "num_cations": 2,
        "label": "Cs2O",
        "cation_atomic_weight": 132.90545196,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 281.80990392,
        "cation_atomic_number": 55,
    },
    "Cu": {
        "num_oxygens": 1,
        "num_cations": 1,
        "label": "CuO",
        "cation_atomic_weight": 63.546,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 79.545,
        "cation_atomic_number": 29,
    },
    "Dy": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Dy2O3",
        "cation_atomic_weight": 162.5,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 372.997,
        "cation_atomic_number": 66,
    },
    "Er": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Er2O3",
        "cation_atomic_weight": 167.259,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 382.515,
        "cation_atomic_number": 68,
    },
    "Eu": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Eu2O3",
        "cation_atomic_weight": 151.964,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 351.925,
        "cation_atomic_number": 63,
    },
    "Fe": {
        "num_oxygens": 1,
        "num_cations": 1,
        "label": "FeOT",
        "cation_atomic_weight": 55.845,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 71.844,
        "cation_atomic_number": 26,
    },
    "Ga": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Ga2O3",
        "cation_atomic_weight": 69.723,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 187.44299999999998,
        "cation_atomic_number": 31,
    },
    "Gd": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Gd2O3",
        "cation_atomic_weight": 157.25,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 362.497,
        "cation_atomic_number": 64,
    },
    "Ge": {
        "num_oxygens": 2,
        "num_cations": 1,
        "label": "GeO2",
        "cation_atomic_weight": 72.63,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 104.628,
        "cation_atomic_number": 32,
    },
    "H": {
        "num_oxygens": 1,
        "num_cations": 2,
        "label": "H2O",
        "cation_atomic_weight": 1.008,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 18.015,
        "cation_atomic_number": 1,
    },
    "Hf": {
        "num_oxygens": 2,
        "num_cations": 1,
        "label": "HfO2",
        "cation_atomic_weight": 178.49,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 210.488,
        "cation_atomic_number": 72,
    },
    "Ho": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Ho2O3",
        "cation_atomic_weight": 164.93033,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 377.85766,
        "cation_atomic_number": 67,
    },
    "K": {
        "num_oxygens": 1,
        "num_cations": 2,
        "label": "K2O",
        "cation_atomic_weight": 39.0983,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 94.1956,
        "cation_atomic_number": 19,
    },
    "La": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "La2O3",
        "cation_atomic_weight": 138.90547,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 325.80794000000003,
        "cation_atomic_number": 57,
    },
    "Li": {
        "num_oxygens": 1,
        "num_cations": 2,
        "label": "Li2O",
        "cation_atomic_weight": 6.94,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 29.879,
        "cation_atomic_number": 3,
    },
    "Lu": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Lu2O3",
        "cation_atomic_weight": 174.9668,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 397.9306,
        "cation_atomic_number": 71,
    },
    "Mg": {
        "num_oxygens": 1,
        "num_cations": 1,
        "label": "MgO",
        "cation_atomic_weight": 24.305,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 40.304,
        "cation_atomic_number": 12,
    },
    "Mn": {
        "num_oxygens": 1,
        "num_cations": 1,
        "label": "MnO",
        "cation_atomic_weight": 54.938044,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 70.937044,
        "cation_atomic_number": 25,
    },
    "Mo": {
        "num_oxygens": 3,
        "num_cations": 1,
        "label": "MoO3",
        "cation_atomic_weight": 95.95,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 143.947,
        "cation_atomic_number": 42,
    },
    "Na": {
        "num_oxygens": 1,
        "num_cations": 2,
        "label": "Na2O",
        "cation_atomic_weight": 22.98976928,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 61.978538560000004,
        "cation_atomic_number": 11,
    },
    "Nb": {
        "num_oxygens": 5,
        "num_cations": 2,
        "label": "Nb2O5",
        "cation_atomic_weight": 92.90637,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 265.80773999999997,
        "cation_atomic_number": 41,
    },
    "Nd": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Nd2O3",
        "cation_atomic_weight": 144.242,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 336.481,
        "cation_atomic_number": 60,
    },
    "Ni": {
        "num_oxygens": 1,
        "num_cations": 1,
        "label": "NiO",
        "cation_atomic_weight": 58.6934,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 74.69239999999999,
        "cation_atomic_number": 28,
    },
    "P": {
        "num_oxygens": 5,
        "num_cations": 2,
        "label": "P2O5",
        "cation_atomic_weight": 30.973761998,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 141.942523996,
        "cation_atomic_number": 15,
    },
    "Pb": {
        "num_oxygens": 1,
        "num_cations": 1,
        "label": "PbO",
        "cation_atomic_weight": 207.2,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 223.19899999999998,
        "cation_atomic_number": 82,
    },
    "Pr": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Pr2O3",
        "cation_atomic_weight": 140.90766,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 329.81232,
        "cation_atomic_number": 59,
    },
    "Rb": {
        "num_oxygens": 1,
        "num_cations": 2,
        "label": "Rb2O",
        "cation_atomic_weight": 85.4678,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 186.9346,
        "cation_atomic_number": 37,
    },
    "S": {
        "num_oxygens": 3,
        "num_cations": 1,
        "label": "SO3",
        "cation_atomic_weight": 32.06,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 80.057,
        "cation_atomic_number": 16,
    },
    "Sb": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Sb2O3",
        "cation_atomic_weight": 121.76,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 291.517,
        "cation_atomic_number": 51,
    },
    "Sc": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Sc2O3",
        "cation_atomic_weight": 44.955908,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 137.908816,
        "cation_atomic_number": 21,
    },
    "Si": {
        "num_oxygens": 2,
        "num_cations": 1,
        "label": "SiO2",
        "cation_atomic_weight": 28.085,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 60.083,
        "cation_atomic_number": 14,
    },
    "Sm": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Sm2O3",
        "cation_atomic_weight": 150.36,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 348.71700000000004,
        "cation_atomic_number": 62,
    },
    "Sn": {
        "num_oxygens": 2,
        "num_cations": 1,
        "label": "SnO2",
        "cation_atomic_weight": 118.71,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 150.708,
        "cation_atomic_number": 50,
    },
    "Sr": {
        "num_oxygens": 1,
        "num_cations": 1,
        "label": "SrO",
        "cation_atomic_weight": 87.62,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 103.619,
        "cation_atomic_number": 38,
    },
    "Ta": {
        "num_oxygens": 5,
        "num_cations": 2,
        "label": "Ta2O5",
        "cation_atomic_weight": 180.94788,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 441.89076,
        "cation_atomic_number": 73,
    },
    "Tb": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Tb2O3",
        "cation_atomic_weight": 158.92535,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 365.84770000000003,
        "cation_atomic_number": 65,
    },
    "Th": {
        "num_oxygens": 2,
        "num_cations": 1,
        "label": "ThO2",
        "cation_atomic_weight": 232.0377,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 264.0357,
        "cation_atomic_number": 90,
    },
    "Ti": {
        "num_oxygens": 2,
        "num_cations": 1,
        "label": "TiO2",
        "cation_atomic_weight": 47.867,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 79.865,
        "cation_atomic_number": 22,
    },
    "Tm": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Tm2O3",
        "cation_atomic_weight": 168.93422,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 385.86544000000004,
        "cation_atomic_number": 69,
    },
    "V": {
        "num_oxygens": 5,
        "num_cations": 2,
        "label": "V2O5",
        "cation_atomic_weight": 50.9415,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 181.878,
        "cation_atomic_number": 23,
    },
    "W": {
        "num_oxygens": 3,
        "num_cations": 1,
        "label": "WO3",
        "cation_atomic_weight": 183.84,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 231.837,
        "cation_atomic_number": 74,
    },
    "Y": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Y2O3",
        "cation_atomic_weight": 88.90584,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 225.80867999999998,
        "cation_atomic_number": 39,
    },
    "Yb": {
        "num_oxygens": 3,
        "num_cations": 2,
        "label": "Yb2O3",
        "cation_atomic_weight": 173.045,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 394.087,
        "cation_atomic_number": 70,
    },
    "Zn": {
        "num_oxygens": 1,
        "num_cations": 1,
        "label": "ZnO",
        "cation_atomic_weight": 65.38,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 81.37899999999999,
        "cation_atomic_number": 30,
    },
    "Zr": {
        "num_oxygens": 2,
        "num_cations": 1,
        "label": "ZrO2",
        "cation_atomic_weight": 91.224,
        "oxygen_atomic_weight": 15.999,
        "molecular_weight": 123.22200000000001,
        "cation_atomic_number": 40,
    },
}


def wt_percent_to_oxide(y: int | float | list | np.ndarray | pd.Series, element: str):
    """convert weight percent element to weight percent oxide

    Parameters
    ----------
    element : str
        element to convert chemistry for
    y : array-like
        data to convert

    Returns
    -------
    weight percent oxide in same shape as input data

    """
    assert (
        element in oxide_dict
    ), f"{element} is not a supported oxide for conversion to parts per million. please use conversions.supported_internal_standard_oxides for supported oxides"
    assert isinstance(
        y, (int, float, list, pd.core.series.Series, np.ndarray)
    ), "y should be something you can do math on - float, integer, list with numbers, numpy array, pandas Series"

    el_dict = oxide_dict[element]

    if isinstance(y, pd.Series):
        y = y.values

    wt_per_el = y

    return (el_dict["molecular_weight"] * wt_per_el) / (
        el_dict["cation_atomic_weight"] * el_dict["num_cations"]
    )


supported_internal_standard_oxides = [
    oxide_dict[e]["label"] for e, i in oxide_dict.items()
]


# def oxide_to_ppm(wt_percent, int_std):
#     """
#     convert concentration internal standard analyte oxide in weight percent to
#     concentration ppm for a 1D series of data

#     Args:
#     wt_percent (array-like): the oxide values to be converted to ppm
#     int_std (str): the internal standard used in the experiment (e.g., '29Si', '43Ca','47Ti')

#     Returns:
#     ppm (array-like): concentrations in ppm the same shape as the wt_percent input

#     """

#     el = [i for i in int_std if not i.isdigit()]

#     if len(el) == 2:
#         element = el[0] + el[1]

#     else:
#         element = el[0]

#     oxides = [
#         "SiO2",
#         "TiO2",
#         "Al2O3",
#         "Cr2O3",
#         "MnO",
#         "FeO",
#         "K2O",
#         "CaO",
#         "Na2O",
#         "NiO",
#         "MgO",
#     ]

#     for o in oxides:
#         if element in o:
#             oxide = o


#     s = oxide.split("O")
#     cat_subscript = s[0]
#     an_subscript = s[1]

#     cat_subscript = [i for i in cat_subscript if i.isdigit()]
#     if cat_subscript:
#         cat_subscript = int(cat_subscript[0])
#     else:
#         cat_subscript = 1

#     an_subscript = [i for i in an_subscript if i.isdigit()]
#     if an_subscript:
#         an_subscript = int(an_subscript[0])
#     else:
#         an_subscript = 1

#     ppm = 1e4 * (
#         (wt_percent * mendeleev.element(element).atomic_weight * cat_subscript)
#         / (
#             mendeleev.element(element).atomic_weight
#             + mendeleev.element("O").atomic_weight * an_subscript
#         )
#     )
#     return ppm


# this now removes the need to have the mendeleev library as a dependency
def oxide_to_ppm(y: int | float | list | np.ndarray | pd.Series, element: str):
    """convert weight percent oxide to parts per million element for select oxides
    see `conversions.supported_internal_standard_oxides` for list of supported oxides

    Args:
        y (int | float | list | np.ndarray | pd.Series):oxide weight percent concentration
        element (str): element to have converted from weight percent oxide to ppm element (this is the cation in the oxide)

    Returns:
        ppm converted values in same type and shape as input
    """

    assert (
        element in oxide_dict
    ), f"{element} is not a supported oxide for conversion to parts per million. please use conversions.supported_internal_standard_oxides for supported oxides"
    assert isinstance(
        y, (int, float, list, pd.core.series.Series, np.ndarray)
    ), "y should be something you can do math on - float, integer, list with numbers, numpy array, pandas Series"

    el_dict = oxide_dict[element]

    if isinstance(y, pd.Series):
        y = y.values

    wt_per_oxide = y

    return (
        1e4 * wt_per_oxide * el_dict["cation_atomic_weight"] * el_dict["num_cations"]
    ) / el_dict["molecular_weight"]
