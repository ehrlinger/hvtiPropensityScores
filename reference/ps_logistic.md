# Estimate propensity scores via binary logistic regression

Fits `stats::glm(..., family = binomial())` and appends the estimated
propensity score, its logit, and a matching weight to the dataset.
Supports multiply-imputed data in a stacked "long" format.

## Usage

``` r
ps_logistic(
  formula,
  data,
  treatment_col = NULL,
  id_col = "id",
  imputation_col = NULL,
  score_col = "prob_t",
  logit_col = "logit_t",
  weight_col = "mt_wt",
  covariates = NULL
)
```

## Arguments

- formula:

  A formula with the binary treatment on the left-hand side, e.g.
  `tavr ~ age + female + ef + diabetes`.

- data:

  A data frame. If `imputation_col` is set, `data` must be a stacked
  "long" MI data frame as produced by
  `mice::complete(mids, action = "long")`.

- treatment_col:

  Name of the binary treatment column (0/1 or logical). If `NULL`
  (default), extracted from the LHS of `formula`.

- id_col:

  Patient identifier column. Required when `imputation_col` is set.
  Default `"id"`.

- imputation_col:

  Name of the imputation-index column in a stacked MI data frame. `NULL`
  (default) means a single complete dataset is supplied.

- score_col:

  Output column name for the propensity score (probability of
  treatment). Default `"prob_t"`.

- logit_col:

  Output column name for the logit of the propensity score,
  `log(p / (1 - p))`. Default `"logit_t"`.

- weight_col:

  Output column name for the matching weight. Default `"mt_wt"`.

- covariates:

  Character vector of covariate column names for SMD balance
  diagnostics. If `NULL` (default), all numeric columns other than
  `treatment_col`, `score_col`, `logit_col`, `weight_col`, `id_col`,
  `"quintile"`, and `"decile"` are used.

## Value

An object of class `c("ps_logistic", "ps_data")` with:

- `$data`:

  The base data frame with `score_col`, `logit_col`, `weight_col`,
  `quintile`, and `decile` columns appended.

- `$meta`:

  Named list: `formula`, `treatment_col`, `id_col`, `imputation_col`,
  `score_col`, `logit_col`, `weight_col`, `method`, `n_imputations`,
  `n_total`.

- `$tables`:

  Named list: `smd`, `group_counts`.

## Details

When `imputation_col` is supplied, the model is fit separately on each
imputed dataset; predicted probabilities are averaged per patient across
imputations (Rubin's combination rule for predictions), and the first
imputed dataset's covariate values form the base output data frame. This
mirrors the SAS template workflow of `PROC LOGISTIC ... BY _IMPUTATION_`
followed by `PROC SUMMARY mean=`.

The **matching weight** appended as `weight_col` is
`min(p, 1-p) / (p * trt + (1-p) * (1-trt))`, which equals the `mt_wt`
column produced in the SAS templates (Li & Greene, 2013). It is the ATM
(average treatment effect among the matched) estimand weight and is used
as input to
[`ps_match()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_match.md)
or
[`ps_weight()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_weight.md).

## See also

[`ps_match()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_match.md),
[`ps_weight()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_weight.md),
[`ps_ordinal()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_ordinal.md),
[`ps_nominal()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/ps_nominal.md),
[`sample_ps_data()`](https://ehrlinger.github.io/hvtiPropensityScores/reference/sample_ps_data.md)

## Examples

``` r
dta <- sample_ps_data(n = 200, seed = 42)

# Single complete dataset
obj <- ps_logistic(
  tavr ~ age + female + ef + diabetes + hypertension,
  data = dta
)
print(obj)
#> <ps_logistic>
#>   N total     : 400
#>   Treatment   : tavr
#>   PS column   : prob_t
#>   Weight col  : mt_wt
#>   Method      : logistic
#>   Tables      : smd, group_counts 
summary(obj)
#> Summary of <ps_logistic>
#> 
#> Balance (SMD):
#>                  variable     smd
#> age                   age -0.7994
#> female             female  0.3989
#> ef                     ef  0.6139
#> diabetes         diabetes -0.0554
#> hypertension hypertension -0.0310
#> match               match      NA
#> 
#> Group counts:
#>     group   n
#> 1 control 200
#> 2 treated 200
#> 

# Pass the scored data to ps_match()
matched <- ps_match(obj$data, score_col = obj$meta$score_col)
nrow(matched$data[matched$data$match == 1L, ])
#> [1] 400
```
