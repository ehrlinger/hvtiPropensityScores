# Assert that required columns exist in a data frame

Assert that required columns exist in a data frame

## Usage

``` r
.check_cols(data, cols, call_env = rlang::caller_env())
```

## Arguments

- data:

  A data frame.

- cols:

  Character vector of required column names.

- call_env:

  Environment for the error call (use
  [`rlang::caller_env()`](https://rlang.r-lib.org/reference/stack.html)).

## Value

Invisible NULL. Signals a classed error on failure.
