"""
conversions module:
For converting wt% oxide to ppm
"""

import mendeleev


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
