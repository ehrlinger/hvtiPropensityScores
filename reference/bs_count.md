# Estimate a balancing score via negative-binomial or Poisson regression

Fits a GLM with a count distribution
([`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) for
negative binomial, or [`stats::glm()`](https://rdrr.io/r/stats/glm.html)
with Poisson family) on a (possibly multiply-imputed) dataset, averages
the per-patient **linear predictors** (log-scale, i.e. xbeta) across
imputations, and assigns rank-based strata. Mirrors
`tp.pm.count.balncing_score.sas` (PROC GENMOD `dist=nb link=log`

- PROC SUMMARY + the `%grp` macro).

## Usage

``` r
bs_count(
  formula,
  data,
  outcome_col = NULL,
  id_col = "id",
  imputation_col = NULL,
  dist = c("negbin", "poisson"),
  score_col = "bs",
  n_strata = 10L,
  strata_col = "cluster",
  create_binary_strata = TRUE,
  covariates = NULL
)
```

## Arguments

- formula:

  A formula with the count outcome on the LHS, e.g.
  `rbc_tot ~ female + agee + in_hct + hx_csurg + iv_cpb`.

- data:

  A data frame, or a stacked MI data frame when `imputation_col` is set.

- outcome_col:

  Name of the count outcome column. `NULL` (default) extracts from
  `formula`.

- id_col:

  Patient identifier column. Required when `imputation_col` is set.
  Default `"id"`.

- imputation_col:

  Imputation-index column. `NULL` (default) means a single complete
  dataset.

- dist:

  Distribution for the GLM. `"negbin"` (default) calls
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) (negative
  binomial, log link); `"poisson"` calls
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) with
  `family = poisson(link = "log")`.

- score_col:

  Name of the output balancing score column (averaged linear predictor
  on the log scale). Default `"bs"`.

- n_strata:

  Number of strata. Default `10L`.

- strata_col:

  Name of the cluster/stratum column. Default `"cluster"`.

- create_binary_strata:

  Logical. If `TRUE` (default), binary indicator columns `stra_1`
  through `stra_<n_strata>` are appended.

- covariates:

  Covariate columns for diagnostics.

## Value

An object of class `c("bs_count", "ps_data")` with:

- `$data`:

  Base data frame with `score_col`, `"quintile"`, `"decile"`,
  `strata_col`, and (if `create_binary_strata`) `stra_1`...`stra_k`
  appended.

- `$meta`:

  Named list: `formula`, `outcome_col`, `id_col`, `imputation_col`,
  `score_col`, `dist`, `n_strata`, `strata_col`, `method`,
  `n_imputations`, `n_total`.

- `$tables`:

  Named list: `strata_counts`.

## Details

The linear predictor (xbeta, log-scale) rather than the fitted count is
used as the balancing score because it is on an unbounded continuous
scale more suitable for rank-based stratification.

## Package dependency

`dist = "negbin"` (the default) requires the MASS package (a recommended
R package, usually pre-installed). Install with
`install.packages("MASS")` if missing, or use `dist = "poisson"`.

## See also

[`bs_continuous()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/bs_continuous.md),
[`ps_logistic()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_logistic.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dta <- sample_ps_data_count(n = 300, seed = 42)
obj <- bs_count(
  rbc_tot ~ female + age + diabetes + hypertension,
  data = dta,
  dist = "poisson"   # or "negbin" if MASS is available
)
print(obj)
table(obj$data$cluster)
} # }
```
