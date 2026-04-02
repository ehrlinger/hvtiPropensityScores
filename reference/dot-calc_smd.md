# Compute standardised mean difference between two numeric vectors

Uses the pooled standard deviation: `SMD = (mean1 - mean0) / sd_pooled`.
Returns `NA_real_` if either group is empty or has zero variance.

## Usage

``` r
.calc_smd(x, treatment)
```

## Arguments

- x:

  Numeric vector of covariate values.

- treatment:

  Binary vector (0/1 or logical) of treatment indicators.

## Value

A single numeric SMD value.
