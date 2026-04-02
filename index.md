# hvtiPropensityScores

Propensity score methods for cardiac surgery comparative-effectiveness
research. Ports the balancing logic from SAS programs into a tidy R API
whose outputs are compatible with
[hvtiPlotR](https://ehrlinger.github.io/hvtiPlotR/) and
[hvtiRutilities](https://github.com/ehrlinger/hvtiRutilities).

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

# --- Nearest-neighbour 1:1 matching ---
m <- ps_match(dta)
print(m)
summary(m)
matched <- m$data[m$data$match == 1L, ]

# --- IPTW weighting (ATE) ---
w <- ps_weight(dta, estimand = "ATE")
print(w)
summary(w)
```

## Functions

| Function                                                                                           | Description                                                 |
|----------------------------------------------------------------------------------------------------|-------------------------------------------------------------|
| [`sample_ps_data()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/sample_ps_data.md) | Synthetic cardiac-surgery dataset for examples and tests    |
| [`ps_match()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_match.md)             | Greedy nearest-neighbour 1:1 matching with optional caliper |
| [`ps_weight()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_weight.md)           | IPTW weighting — ATE, ATT, or ATC estimand                  |
| [`is_ps_data()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/is_ps_data.md)         | Predicate for `ps_data` objects                             |

## Related packages

- **hvtiPlotR** — mirror histograms, balance plots, survival and hazard
  curves
- **hvtiRutilities** — data wrangling utilities for cardiac surgery
  datasets
