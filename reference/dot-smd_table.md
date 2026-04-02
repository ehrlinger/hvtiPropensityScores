# Build a tidy SMD table across all numeric covariates

Build a tidy SMD table across all numeric covariates

## Usage

``` r
.smd_table(data, treatment, covariates = NULL)
```

## Arguments

- data:

  A data frame.

- treatment:

  Name of the binary treatment column.

- covariates:

  Character vector of covariate column names. If `NULL`, all numeric
  columns other than `treatment` are used.

## Value

A data frame with columns `variable` and `smd`.
