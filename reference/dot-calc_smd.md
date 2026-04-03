# Compute standardised mean difference between two numeric vectors

Uses the sample-size-weighted pooled standard deviation:
`sp = sqrt(((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2))`. This is correct
for unequal group sizes (e.g. pre-match diagnostics). The equal-weight
formula `sqrt((s0^2 + s1^2) / 2)` is only unbiased when n0 == n1.

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

## Details

Returns `NA_real_` if either group has fewer than 2 observations or the
pooled SD is zero (constant covariate).
