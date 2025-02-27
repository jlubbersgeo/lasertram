# Welcome to the `lasertram` Documentation!

This site contains the project documentation for `lasertram` (laser <b><i>T</i></b>ime <b><i>R</i></b>esolved <b><i>A</i></b>nalysis <b><i>M</i></b>odule): a lightweight, but powerful package for processing laser ablation inductively coupled plasma mass spectrometry (LA-ICP-MS) data. It is designed to give the user complete control over which portion of their LA-ICP-MS data is converted into concentrations. The underlying logic is largely inspired from two references:

1. [Longerich, H. P., Jackson, S. E., & Günther, D. (1996). Inter-laboratory note. Laser ablation inductively coupled plasma mass spectrometric transient signal data acquisition and analyte concentration calculation. Journal of analytical atomic spectrometry, 11(9), 899-904.](https://pubs.rsc.org/en/content/articlepdf/1996/ja/ja9961100899)
2. [Kent, A. J., & Ungerer, C. A. (2006). Analysis of light lithophile elements (Li, Be, B) by laser ablation ICP-MS: comparison between magnetic sector and quadrupole ICP-MS. American Mineralogist, 91(8-9), 1401-1411.](https://pubs.geoscienceworld.org/msa/ammin/article/91/8-9/1401/134304/Analysis-of-light-lithophile-elements-Li-Be-B-by)

## Motivation

With a wide array of applications in the natural sciences, laser ablation inductively coupled plasma mass spectrometry (LA-ICP-MS) is a now a commonplace tool for the gathering of <i>in situ</i> trace element (i.e., $<$ 0.1 wt%) data from solid materials. The last two decades have seen significant advances in both instrument capabilities and operating software, allowing users to generate large volumes of <i>in situ</i> geochemical data in comparatively little time to previous methodologies (i.e., micro-drilling) while still maintaining high degrees of accuracy and precision.

Raw data output from LA-ICP-MS, however, is in the form of counts per second (cps) for the selected analyte isotopes, not elemental concentrations. In order to be converted into accurate concentrations, a modest amount of user input and interpretation is required and should not be automated. That does not mean, however, we cannot streamline our workflow! This documuentation is designed to outline the theory, workflow, and structure, behind `lasertram` in an effort to maximize its effectiveness in the petrology and volcanology communities.

Overall the goals of `lasertram` are three-fold:

1. provide an intuitive and efficient workflow for processing LA-ICP-MS data that takes advantage of the scientific python ecosystem
2. make LA-ICP-MS data processing as reproducible as possible
3. Allow processing by experts and non-experts alike. Defaults will get you there but the API allows for flexibility for advanced users to get into the `weeds`
   > More time doing science less time processing data!!

## Projects

- `lasertram` is now published in _Applied Computing and Geosciences_:

```tex
@article{lubbers2025lasertram,
  title={lasertram: A Python library for time resolved analysis of laser ablation inductively coupled plasma mass spectrometry data},
  author={Lubbers, Jordan and Kent, Adam JR and Russo, Chris},
  journal={Applied Computing and Geosciences},
  pages={100225},
  year={2025},
  publisher={Elsevier},
  doi={10.1016/j.acags.2025.100225}
}
```

- `lasertram` is the backend engine to [LaserTRAM-DB](https://github.com/jlubbersgeo/laserTRAM-DB), a dashboard for interactively and transparently processing LA-ICP-MS data.

## Get started

```
pip install lasertram
```

## Table of Contents

1. [Background](explanation.md)
2. [Tutorials - Basic](lasertram_tutorial.ipynb)
3. [Tutorials - Signal despiking](lasertram_despiking.ipynb)
4. [Tutorials - Region omission](lasertram_region_omission.ipynb)
5. [API reference](reference.md)
