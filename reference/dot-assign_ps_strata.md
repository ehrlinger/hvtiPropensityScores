# Assign rank-based quintile and decile columns

Mirrors the SAS logic: proc sort by \_propen; quintile = int(*n* /
(nobs/5)) + 1; if quintile \> 5 then quintile = 5;

## Usage

``` r
.assign_ps_strata(score)
```
