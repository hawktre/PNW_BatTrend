# README - Posteriors-for-UpdatingPriors-BatModels.csv

**Project Name**: `BatHub\_Trend`  
**GitHub Location**: *TBD*

## Overview
This file describes the posterior means and sd contained in the file `Posteriors-for-UpdatingPriors-BatModels.csv`.

These data were generated from the posterior means and sd estimated during previous analyses to be used as
Gaussian priors in subsequent models with new information.

## Source documentation

Rodhouse, T. J., P. C. Ormsbee, K. M. Irvine, L. A. Vierling, J. M. Szewczak, and K. E. Vierling.
2015. Establishing conservation baselines with dynamic distribution models for bat populations
facing imminent decline. Diversity and Distributions 21:1401–1413.

Rodhouse, T. J., R. M. Rodriguez, K. M. Banner, P. C. Ormsbee, J. Barnett, and K. M. Irvine.
2019. Evidence of region-wide bat population decline from long-term monitoring and Bayesian occupancy models
with empirically informed priors. Ecology and Evolution 9:11078–11088.

Banner, K. M., K. M. Irvine, and T. J. Rodhouse.
2020. The use of Bayesian priors in ecology: the good, the bad and the not great.
Methods in Ecology and Evolution 11:882–889.

## Details about the posteriors

1. No posteriors are included for the detection model parameters because methods changed significantly as did
model structure.

2. Informative priors on occupancy parameters are described in detail by Rodhouse et al. 2019 Table 1. Additional
information can be found in Banner et al. 2020.

3. The .csv includes posteriors from the Rodhouse et al. 2015 analyses as a starting baseline for updating priors
with new data. It also incudes more recent updated posteriors for MYLU and LACI from Rodhouse et al. 2019. Note
that no information is available for TABR (see Rodhouse et al. 2015).

4. Columns include posterior means and sd. Priors from these posteriors are Gaussian.

### Data Dictionary for `Posteriors-for-UpdatingPriors-BatModels.csv`:

| Column | Original Model Notation | Source | Description |
|-|-|-|-|
|`Intercept`|Beta0|Rodhouse et al. 2019 (Table 1)|Notated as **alpha01** in the `.jags` models used in this analysis repository.|
|`Gamma`|alpha|Rodhouse et al. 2019 (Table 1)|Colonization parameter; notated as **gamma** in the `.jags` model.|
|`Phi`|beta|Rodhouse et al. 2019 (Table 1)|Survival parameter; notated as **phi** in the `.jags` model.|
|`Elevation`||*see previous analyses*|Covariate representing grid cell % cover|
|`Precipitation`||*see previous analyses*|Covariate representing grid cell % cover|
|`Forest`||*see previous analyses*|Covariate representing grid cell % cover|
|`Cliffs`||*see previous analyses*|Covariate representing grid cell % cover|

\###########

