# Assert that a column contains only 0/1 (or logical) values

Emits a warning (not an error) if NA values are present, because NA
treatment assignments silently drop patients from downstream matching
and weighting steps.

## Usage

``` r
.check_binary(data, col, call_env = rlang::caller_env())
```

## Arguments

- data:

  A data frame.

- col:

  Name of the column to check.

- call_env:

  Environment for the error call.

## Value

Invisible NULL.
