# Background Theory

## Calculating concentrations

We calculate the concentration of analyte ($i$) in an unknown material ($u$) using the following relationship from Longerich et al., (1996):

$$
{C_i}^u = \frac{{R_i}^u}{S}
$$

Where ${C_i}^u$ and ${R_i}^u$ are the concentration of analyte and count rate of analyte ($i$) in the unknown material, respectively, and $S$ is the normalized sensitivity. When using naturally occuring internal standards, $S$ can be defined as:

$$
S = \frac{{R_i}^{std}}{{C_i}^{std}}\left[\frac{{R_{n}}^u}{{R_{n}}^{std}} \frac{{C_{n}}^{std}}{{C_{n}}^{u}} \right]
$$

${R_i}^{std}$ and ${C_i}^{std}$ are the count rate and and concentration of analyte ($i$) in the calibration standard, ${R_{n}}^u$ and ${R_{n}}^{std}$ are the mean count rates of the internal standard in the unknown material and calibration standard, ${C_{n}}^{u}$ and ${C_{n}}^{std}$ are the concentrations of the internal standard in the unknown material and calibration standard.

Kent and Ungerer (2006) re-arrange this relationship such that the count rate expressions always contain unknown analytes in the numerator:

$$
{C_i}^u = {C_n}^u \frac{\left[\frac{{C_i}^{std}}{{C_n}^{std}}\right]}{\left[\frac{{R_i}^{std}}{{R_n}^{std}}\right]}\frac{{R_i}^u}{{R*{n}}^u}
$$

## Normalizing to an internal standard

The purpose of `lasertram` is to give the user complete control over which portion of the analytical spectra gets used in calculating concentrations (e.g., filtering out portions of the signal not reflective of the material under investigation). In complex natural materials, selection of this interval and an overall judgement about data quality require an operator to make a decision. This software is optimized to allow that decision to be made as efficient as possible.

When a given interval from the analytical spectra has been chosen, every analyte is normalized to a chosen internal standard. `lasertram` allows for any analyte in the experiment to be used as the internal standard (see caveats on this in the [How-To](how-to-guides.md) section). Prior to normalization to an internal standard, raw data first has the background analyte levels subtracted from it. Background is determined by taking the median counts per second value for each analyte over the user specified background range. Once data have been background subtracted, each normalized ratio is calculated the following way:

$$
N_i = median\left[\frac{cps_{i}}{cps_{is}}\right]
$$

Where $cps_i$ is the background subtracted counts per second data for analyte ($i$), and $cps_{is}$ is the background subtracted counts per second data for the internal standard. Since counts per second is analogous to count rate above in Equation 3, we can simplify the above relationship to now reflect our $N_i$ values:

$$
{C_i}^u = {C_n}^u \frac{\left[\frac{{C_i}^{std}}{{C_n}^{std}}\right]}{{N_{i}}^{std}}{N_{i}}^u
$$

Here, ${N_{i}}^{std}$ and ${N_{i}}^{u}$ are the normalized counts per second value of analyte $i$ in the calibration standard and unknown, respectively. The uncertainty for any given normalized ratio is expressed as:

$$
SE = \frac{\sigma _{N_i}}{\sqrt{n}}
$$

$\sigma_N$ is the standard deviation of a given analyte's normalized ratio for the interval and $n$ is the number of time steps in the interval (i.e., cycles through the mass range). The relative standard error is then:

$$
{RSE_i}^u = \left[\frac{SE}{N_i}\right]100
$$

## Detection Limits

Count rate detection limits for each analyte are determined similar to Longerich et al., (1996):

$$
LOD = 3\sigma\ _{b} + b
$$

where $b$ is the background median counts per second value for a given analyte. This is standard practice in LA-ICP-MS data reduction. To reflect this in data output, measurements that are below detection limit will have values that say "b.d.l." rather than concentrations.

## Drift correction

To check for drift in calibration standard normalized ratios over time, a linear regression is applied to the calibration standard for each analyte, where the dependent variable is the count rate normalized to the internal standard and the independent variable is the timestamp associated with each analysis.

We determine the significance of each regression by evaluating following null hypothesis: there is no relationship between a given analyte's internal standard normalized ratio and time. We reject this if both the following conditions are true: The p-value for the coefficient (i.e., slope) is significant; The F-statisic comparing the regression and observed data is greater than the critical F value. By default, we set the threshold for p-value significance at .01 (i.e., we have 99\% confidence that we can reject the null hypothesis) in an effort to mitigate drift correcting all but the most linear of changes in normalized count rates, but this may be changed by the user. If the null hypothesis for a given analyte is rejected, the analyte is linearly corrected for drift and the regression parameters (e.g., slope and intercept) are used to calculate a normalized count rate for the calibration standard at the point in time where an unknown was analyzed:

$$
{C_i}^u = {C_n}^u \frac{\left[\frac{{C_i}^{std}}{{C_n}^{std}}\right]}{\left[m_ix+b_i\right]}{N_i}^u
$$

Where $m$ is the regression slope, $x$ is the analysis time, and $b$ is the intercept for analyte $i$.

## Uncertainties

Calculating concentrations of a given analyte in an unknown material can be considered a series of nested quotients and products. Therefore, we quantify the overall uncertainty of a given analyte as (Taylor, 1997):

$$
\sigma_{{C_i}} = {C_i}^u \sqrt{ \left( \frac{\sigma_{{C_n}^{u}}}{{C_n}^{u}}\right)^2 + \left( \frac{\sigma_{{C_i}^{std}}}{{C_i}^{std}}\right)^2 + \left( \frac{\sigma_{{C_n}^{std}}}{{C_n}^{std}}\right)^2 + \left({RSE_i}^{std}\right)^2 + \left({RSE_i}^{u}\right)^2}
$$

Where ${RSE_i}^{std}$ is defined as:

$$
{RSE\_{i}}^{std} = \left[\frac{\frac{\sigma_i}{\sqrt{n_i}}}{\mu_i}\right]100
$$

$\sigma_i$ and $\mu_i$ are the standard deviation and mean of all of the calibration standard normalized ratios respectively and $n_i$ is the total number of calibration standard analyses for analyte ($i$).

For analytes where drift correction has been applied, ${RSE_i}^{std}$ is replaced with:

$$
100\left[\frac{RMSE_i}{\mu_i}\right]
$$

Where $RMSE_i$ is the Root Mean Squared Error as specified in the Drift Correction section.

## Concentrations of internal standard in unknown

To calculate concentrations of a given analyte list in an unknown sample, the concentration of the internal standard must be known. LaserCalc takes these concentrations in the form of wt\% oxide and utilizes user interaction to input concentrations of the internal standard and its relative uncertainty. A default value of 10\% is used for this, but may and <b>should</b> be updated by the user.

## References

1. [Longerich, H. P., Jackson, S. E., & Günther, D. (1996). Inter-laboratory note. Laser ablation inductively coupled plasma mass spectrometric transient signal data acquisition and analyte concentration calculation. Journal of analytical atomic spectrometry, 11(9), 899-904.](https://pubs.rsc.org/en/content/articlepdf/1996/ja/ja9961100899)
2. [Kent, A. J., & Ungerer, C. A. (2006). Analysis of light lithophile elements (Li, Be, B) by laser ablation ICP-MS: comparison between magnetic sector and quadrupole ICP-MS. American Mineralogist, 91(8-9), 1401-1411.](https://pubs.geoscienceworld.org/msa/ammin/article/91/8-9/1401/134304/Analysis-of-light-lithophile-elements-Li-Be-B-by)
