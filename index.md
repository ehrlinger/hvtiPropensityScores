# hvtiPropensityScores

Propensity score methods for cardiac surgery comparative-effectiveness
research. Ports the balancing logic from SAS programs into a tidy R API
whose outputs are compatible with
[hvtiPlotR](https://ehrlinger.github.io/hvtiPlotR/) and
[hvtiRutilities](https://github.com/ehrlinger/hvtiRutilities).

**Full documentation:**
<https://ehrlinger.github.io/hvtiPropensityScores/>

## Installation

``` r
# Install the development version from GitHub:
# install.packages("pak")
pak::pak("ehrlinger/hvtiPropensityScores")
```

## Quick start

``` r
library(hvtiPropensityScores)

# Reproducible synthetic dataset
dta <- sample_ps_data(n = 500, seed = 42)

# --- Propensity score estimation ---
obj <- ps_logistic(tavr ~ age + female + ef + diabetes + hypertension, data = dta)

# --- Nearest-neighbour 1:1 matching ---
m <- ps_match(obj$data)
print(m)
summary(m)
matched <- m$data[m$data$match == 1L, ]

# --- IPTW weighting (ATE) ---
w <- ps_weight(dta, estimand = "ATE")
print(w)
summary(w)
```

## Functions

**Score estimation**

| Function                                                                                         | Description                                                         |
|--------------------------------------------------------------------------------------------------|---------------------------------------------------------------------|
| [`ps_logistic()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_logistic.md)     | Binary propensity score via logistic regression (with optional MI)  |
| [`ps_ordinal()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_ordinal.md)       | Ordered treatment propensity via cumulative logit (MASS::polr)      |
| [`ps_nominal()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_nominal.md)       | Nominal treatment propensity via generalised logit (nnet::multinom) |
| [`bs_continuous()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/bs_continuous.md) | Balancing score for continuous exposures via linear regression      |
| [`bs_count()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/bs_count.md)           | Balancing score for count exposures via negative-binomial / Poisson |

**Balancing methods**

| Function                                                                                 | Description                                                 |
|------------------------------------------------------------------------------------------|-------------------------------------------------------------|
| [`ps_match()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_match.md)   | Greedy nearest-neighbour 1:1 matching with optional caliper |
| [`ps_weight()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_weight.md) | IPTW weighting — ATE, ATT, or ATC estimand                  |

**Sensitivity analysis**

| Function                                                                                         | Description                                                   |
|--------------------------------------------------------------------------------------------------|---------------------------------------------------------------|
| [`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/sa_rosenbaum.md)   | Rosenbaum gamma bounds for matched analyses                   |
| [`sa_evalue()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/sa_evalue.md)         | E-values (VanderWeele & Ding 2017) for unmeasured confounding |
| [`sa_overlap()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/sa_overlap.md)       | PS overlap / positivity diagnostics                           |
| [`sa_trim_sweep()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/sa_trim_sweep.md) | IPTW trim-threshold sensitivity sweep                         |

**Utilities**

| Function                                                                                                           | Description                                          |
|--------------------------------------------------------------------------------------------------------------------|------------------------------------------------------|
| [`sample_ps_data()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/sample_ps_data.md)                 | Synthetic cardiac-surgery dataset (binary treatment) |
| [`sample_ps_data_ordinal()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/sample_ps_data_ordinal.md) | Synthetic dataset with ordered treatment             |
| [`sample_ps_data_nominal()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/sample_ps_data_nominal.md) | Synthetic dataset with nominal treatment             |
| [`sample_ps_data_count()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/sample_ps_data_count.md)     | Synthetic dataset with count exposure                |
| [`is_ps_data()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/is_ps_data.md)                         | Predicate for `ps_data` objects                      |

## Related packages

- **hvtiPlotR** — mirror histograms, balance plots, survival and hazard
  curves
- **hvtiRutilities** — data wrangling utilities for cardiac surgery
  datasets
