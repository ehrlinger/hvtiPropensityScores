# Assert that a column contains only 0/1 (or logical) values

Assert that a column contains only 0/1 (or logical) values

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
