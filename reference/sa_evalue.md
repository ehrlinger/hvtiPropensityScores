# E-Value Sensitivity Analysis

Computes the E-value for a point estimate and optionally for a
confidence interval bound, following VanderWeele & Ding (2017). The
E-value is the minimum strength of association (on the risk-ratio scale)
that an unmeasured confounder would need with **both** the treatment and
the outcome, conditional on measured covariates, to fully explain away
the observed association.

## Usage

``` r
sa_evalue(
  estimate,
  ci_lo = NULL,
  ci_hi = NULL,
  type = c("RR", "OR", "HR", "IRR", "RD"),
  p0 = NULL
)
```

## Arguments

- estimate:

  Numeric. The point estimate (RR, OR, HR, or RD on the original scale).
  Must be \> 0 for RR/OR/HR; any finite value for RD.

- ci_lo:

  Optional numeric. Lower confidence interval bound.

- ci_hi:

  Optional numeric. Upper confidence interval bound. Both `ci_lo` and
  `ci_hi` must be supplied together to compute the CI E-value.

- type:

  One of `"RR"` (default), `"OR"`, `"HR"`, `"IRR"`, or `"RD"`. `"IRR"`
  (incidence rate ratio) uses the same formula as `"RR"`; it is provided
  as a distinct value for clarity in count-outcome studies.

- p0:

  Baseline outcome probability in the unexposed group. Required when
  `type = "RD"`. Must be in (0, 1).

## Value

A list with:

- `$evalue_estimate`:

  E-value for the point estimate.

- `$evalue_ci`:

  E-value for the confidence interval bound closest to the null. `NA` if
  `ci_lo`/`ci_hi` are not supplied. Returns `1` if the CI already
  crosses the null.

- `$rr_estimate`:

  The point estimate converted to the risk-ratio scale (same as
  `estimate` when `type = "RR"`).

- `$type`:

  The measure type supplied.

## Details

For odds ratios and hazard ratios the rare-outcome approximation
`RR_approx = sqrt(OR)` / `sqrt(HR)` is used to convert to the risk-ratio
scale before applying the E-value formula. For risk differences a
baseline outcome probability `p0` must be supplied.

## References

VanderWeele, T. J., & Ding, P. (2017). Sensitivity analysis in
observational research: Introducing the E-value. *Annals of Internal
Medicine*, 167(4), 268-274.
[doi:10.7326/M16-2607](https://doi.org/10.7326/M16-2607)

## See also

[`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/sa_rosenbaum.md),
[`sa_overlap()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/sa_overlap.md),
[`sa_trim_sweep()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/sa_trim_sweep.md)

## Examples

``` r
# Risk ratio of 2.5 -- how robust is the result?
sa_evalue(2.5, ci_lo = 1.8, ci_hi = 3.5)
#> $evalue_estimate
#> [1] 4.436492
#> 
#> $evalue_ci
#> [1] 3
#> 
#> $rr_estimate
#> [1] 2.5
#> 
#> $type
#> [1] "RR"
#> 

# Odds ratio of 0.6 (protective effect)
sa_evalue(0.6, ci_lo = 0.4, ci_hi = 0.9, type = "OR")
#> $evalue_estimate
#> [1] 1.903916
#> 
#> $evalue_ci
#> [1] 1.292878
#> 
#> $rr_estimate
#> [1] 0.7745967
#> 
#> $type
#> [1] "OR"
#> 

# Risk difference of 0.10 with baseline probability 0.20
sa_evalue(0.10, ci_lo = 0.02, ci_hi = 0.18, type = "RD", p0 = 0.20)
#> $evalue_estimate
#> [1] 2.366025
#> 
#> $evalue_ci
#> [1] 1.431662
#> 
#> $rr_estimate
#> [1] 1.5
#> 
#> $type
#> [1] "RD"
#> 
```
