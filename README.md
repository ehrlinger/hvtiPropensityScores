# hvtiPropensityScores

<!-- badges: start -->
[![R-CMD-check](https://github.com/ehrlinger/hvtiPropensityScores/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ehrlinger/hvtiPropensityScores/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/ehrlinger/hvtiPropensityScores/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ehrlinger/hvtiPropensityScores?branch=main)
[![pkgdown site](https://img.shields.io/badge/docs-pkgdown-blue)](https://ehrlinger.github.io/hvtiPropensityScores/)
<!-- badges: end -->

Propensity score methods for cardiac surgery comparative-effectiveness
research.  Ports the balancing logic from SAS programs into a tidy R API
whose outputs are compatible with
[hvtiPlotR](https://ehrlinger.github.io/hvtiPlotR/) and
[hvtiRutilities](https://github.com/ehrlinger/hvtiRutilities).

**Full documentation:** <https://ehrlinger.github.io/hvtiPropensityScores/>

## Installation

```r
# Install the development version from GitHub:
# install.packages("pak")
pak::pak("ehrlinger/hvtiPropensityScores")
```

## Quick start

```r
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

| Function | Description |
|---|---|
| `ps_logistic()` | Binary propensity score via logistic regression (with optional MI) |
| `ps_ordinal()` | Ordered treatment propensity via cumulative logit (MASS::polr) |
| `ps_nominal()` | Nominal treatment propensity via generalised logit (nnet::multinom) |
| `bs_continuous()` | Balancing score for continuous exposures via linear regression |
| `bs_count()` | Balancing score for count exposures via negative-binomial / Poisson |

**Balancing methods**

| Function | Description |
|---|---|
| `ps_match()` | Greedy nearest-neighbour 1:1 matching with optional caliper |
| `ps_weight()` | IPTW weighting — ATE, ATT, or ATC estimand |

**Sensitivity analysis**

| Function | Description |
|---|---|
| `sa_rosenbaum()` | Rosenbaum gamma bounds for matched analyses |
| `sa_evalue()` | E-values (VanderWeele & Ding 2017) for unmeasured confounding |
| `sa_overlap()` | PS overlap / positivity diagnostics |
| `sa_trim_sweep()` | IPTW trim-threshold sensitivity sweep |

**Utilities**

| Function | Description |
|---|---|
| `sample_ps_data()` | Synthetic cardiac-surgery dataset (binary treatment) |
| `sample_ps_data_ordinal()` | Synthetic dataset with ordered treatment |
| `sample_ps_data_nominal()` | Synthetic dataset with nominal treatment |
| `sample_ps_data_count()` | Synthetic dataset with count exposure |
| `is_ps_data()` | Predicate for `ps_data` objects |

## Related packages

- **hvtiPlotR** — mirror histograms, balance plots, survival and hazard curves
- **hvtiRutilities** — data wrangling utilities for cardiac surgery datasets
