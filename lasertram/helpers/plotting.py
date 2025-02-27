import re

import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1.axes_size as Size
from mpl_toolkits.axes_grid1 import Divider


def plot_timeseries_data(
    df,
    analytes="all",
    marker="",
    fig=None,
    ax=None,
    **kwargs,
):
    """Plot time-series data related to laser ablation ICP-MS analyses,
        typically where the x-axis is analysis time and y-axis is either
        counts per second data or data derived from it.

    Parameters
    ----------
    df : pandas DataFrame
        the dataframe to be plotted
    analytes : str, optional
        list of columns to be plotted from the dataframe, by default 'all'.
        Meant to be utilized when the input dataframe is either a LaserTRAM spot
        object so columns reflect only 'Time' and analytes.
    marker : str, optional
        matplotlib marker to use for plotting symbol, by default ''
    fig : matplotlib.Figure, optional
        The figure to apply the plot to, by default None
    ax : matplotlib.Axes, optional
        the axis to apply the plot to, by default None

    Returns
    -------
    ax


    Ex:
    ```python
    from lasertram import preprocessing, plotting, LaserTRAM
    import matplotlib.pyplot as plt
    plt.style.use("lasertram.lasertram")

    raw_data  = preprocessing.load_test_rawdata()

    sample = 'GSD-1G_-_1'

    ax = plotting.plot_timeseries_data(raw_data.loc[sample,:])
    ax[0].set_title(sample)
    ax[0].set_ylabel("cps")
    ax[0].set_xlabel("Time (ms)")
    ```
    """

    if fig is None:
        fig = plt.figure(figsize=(8, 4))
    else:
        fig = plt.gcf()

    if ax is None:
        # setting up default axes
        rect = (0.1, 0.1, 0.8, 0.8)
        ax = [fig.add_axes(rect, label=f"{i}") for i in range(2)]

        horiz = [Size.AxesX(ax[0]), Size.Fixed(0.5), Size.AxesX(ax[1])]
        vert = [Size.AxesY(ax[0]), Size.Fixed(0.5), Size.AxesY(ax[1])]

        # divide the Axes rectangle into grid whose size is specified by horiz * vert
        divider = Divider(fig, rect, horiz, vert, aspect=False)
        ax[0].set_axes_locator(divider.new_locator(nx=0, ny=0))
        ax[1].set_axes_locator(divider.new_locator(nx=2, ny=0))

    if analytes == "all":
        analytes = [
            column
            for column in df.columns
            if ("timestamp" not in column) and ("Time" not in column)
        ]

        df.loc[:, ["Time"] + analytes].plot(
            x="Time",
            y=analytes,
            kind="line",
            marker=marker,
            ax=ax[0],
            lw=1,
            legend=False,
            **kwargs,
        )

    else:
        if isinstance(analytes, list):
            pass
        else:
            analytes = [analytes]

        df.loc[:, ["Time"] + analytes].plot(
            x="Time",
            y=analytes,
            kind="line",
            marker=marker,
            ax=ax[0],
            lw=1,
            legend=False,
            **kwargs,
        )

    ax[0].set_yscale("log")

    handles, labels = ax[0].get_legend_handles_labels()
    cols = 2
    ax[1].legend(
        handles, labels, loc="upper left", bbox_to_anchor=(0.15, 1.1), ncol=cols
    )
    ax[1].axis("off")

    return ax


def plot_lasertram_uncertainties(spot, fig=None, ax=None, **kwargs):
    """plot a bar chart of analyte uncertainties related to the output from
    processing using the `LaserTRAM` module

    Parameters
    ----------
    spot : LaserTRAM.spot
        the `LaserTRAM.spot` object to plot the uncertainties for
    fig : matplotlib.Figure, optional
        The figure to apply the plot to, by default None
    ax : matplotlib.Axes, optional
        the axis to apply the plot to, by default None

    Returns
    -------
    ax
    """

    if fig is None:
        fig = plt.figure(figsize=(12, 3))
    else:
        fig = plt.gcf()

    if ax is None:
        ax = fig.add_subplot()

    ax.bar(x=spot.analytes, height=spot.bkgd_subtract_std_err_rel, **kwargs)

    labels = [analyte for analyte in spot.analytes]
    labels = [
        "$^{{{}}}${}".format(
            re.findall(r"\d+", label)[0],
            label.replace(re.findall(r"\d+", label)[0], ""),
        )
        for label in labels
    ]
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(labels, rotation=90)
    ax.set_ylabel("% SE")

    return ax
