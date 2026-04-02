# Assert that a column contains valid probabilities (0 \<= x \<= 1)

Assert that a column contains valid probabilities (0 \<= x \<= 1)

## Usage

``` r
.check_probability(data, col, call_env = rlang::caller_env())
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
