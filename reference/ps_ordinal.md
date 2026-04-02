# Estimate propensity scores via ordinal (cumulative logit) logistic regression

Fits a proportional-odds cumulative logit model via
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html) for an ordered
treatment variable and appends one probability column per treatment
level to the dataset. The model is equivalent to the default
`PROC LOGISTIC` behaviour for ordinal responses in SAS.

## Usage

``` r
ps_ordinal(
  formula,
  data,
  treatment_col = NULL,
  id_col = "id",
  imputation_col = NULL,
  score_col_prefix = "prob",
  covariates = NULL
)
```

## Arguments

- formula:

  A formula with the ordered factor treatment on the LHS, e.g.
  `nyha_grp ~ age + female + ef`. The response is coerced to an ordered
  factor if it is not already one.

- data:

  A data frame (single or stacked MI).

- treatment_col:

  Name of the ordinal treatment column. `NULL` (default) extracts from
  `formula`.

- id_col:

  Patient identifier column. Default `"id"`.

- imputation_col:

  Imputation-index column for stacked MI data. `NULL` (default) assumes
  a single complete dataset.

- score_col_prefix:

  Prefix for the per-level output columns. Columns are named
  `<prefix>_<level>` for each level of the treatment. Default `"prob"`.

- covariates:

  Covariate columns for diagnostics.

## Value

An object of class `c("ps_ordinal", "ps_data")` with:

- `$data`:

  Base data frame with one `prob_<level>` column per treatment level.

- `$meta`:

  Named list: `formula`, `treatment_col`, `id_col`, `imputation_col`,
  `score_cols`, `levels`, `method`, `n_imputations`, `n_total`.

- `$tables`:

  Named list: `group_counts`.

## Details

When `imputation_col` is supplied, models are fit on each imputed
dataset separately and per-level probabilities are averaged per patient
across imputations.

## Package dependency

Requires the MASS package (a recommended package shipped with base R).
Install with `install.packages("MASS")` if missing.

## See also

[`ps_logistic()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_logistic.md),
[`ps_nominal()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_nominal.md),
[`ps_match()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_match.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dta <- sample_ps_data_ordinal(n = 300, seed = 42)
obj <- ps_ordinal(
  nyha_grp ~ age + female + ef + diabetes,
  data = dta
)
print(obj)
head(obj$data[, c("id", "nyha_grp", "prob_I", "prob_II", "prob_III")])
} # }
```
