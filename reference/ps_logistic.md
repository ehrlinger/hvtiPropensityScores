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

# --- Single complete dataset (mirrors tp.lm.logistic_propensity_score.nomi.sas)
# Equivalent to: PROC LOGISTIC data=built descending; model tavr = ...;
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
#> Smd:
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

# The function appends:
#   prob_t   -- propensity score (p-hat; SAS _p_ / _propen_)
#   logit_t  -- log(p/(1-p))    (SAS _logit_)
#   mt_wt    -- matching weight  (SAS mt_wt = min(p,1-p)/(p*trt+(1-p)*(1-trt)))
#   quintile -- rank-based quintile (SAS int(_n_/(nobs/5))+1)
#   decile   -- rank-based decile
head(obj$data[, c("id", "tavr", "prob_t", "logit_t", "mt_wt",
                  "quintile", "decile")])
#>   id tavr    prob_t    logit_t     mt_wt quintile decile
#> 1  1    0 0.1112484 -2.0780524 0.1251738        1      1
#> 2  2    0 0.5270896  0.1084644 1.0000000        3      6
#> 3  3    0 0.3951929 -0.4255358 0.6534196        2      4
#> 4  4    0 0.2836606 -0.9263756 0.3959863        2      3
#> 5  5    0 0.3751620 -0.5101344 0.6004149        2      4
#> 6  6    0 0.3148051 -0.7777498 0.4594387        2      3

# Pass the scored data to ps_match() for downstream matching
matched <- ps_match(obj$data, score_col = obj$meta$score_col, seed = 42)
nrow(matched$data[matched$data$match == 1L, ])
#> [1] 400

# --- Multiply-imputed data (mirrors tp.lm.logistic_propensity_score.sas)
# Equivalent to:
#   PROC LOGISTIC data=built descending; BY _IMPUTATION_; model tavr = ...;
#   PROC SUMMARY data=decile; class ccfid; var _p_; output out=... mean=_propen;
# \donttest{
# Simulate a stacked MI dataset (2 imputations, column "_Imputation_")
dta_mi <- rbind(
  cbind(dta, `_Imputation_` = 1L),
  cbind(dta, `_Imputation_` = 2L)
)
names(dta_mi)[names(dta_mi) == "_Imputation_"] <- "imp"

obj_mi <- ps_logistic(
  tavr ~ age + female + ef + diabetes + hypertension,
  data           = dta_mi,
  imputation_col = "imp",
  id_col         = "id"
)
print(obj_mi)
#> <ps_logistic>
#>   N total     : 400
#>   Treatment   : tavr
#>   PS column   : prob_t
#>   Weight col  : mt_wt
#>   Method      : logistic-MI (2 imputations)
#>   Tables      : smd, group_counts 
# Per-patient PS is the average across the two imputed-dataset predictions,
# matching the PROC SUMMARY mean= step in the SAS template.
head(obj_mi$data[, c("id", "tavr", "prob_t")])
#>   id tavr    prob_t
#> 1  1    0 0.1112484
#> 2  2    0 0.5270896
#> 3  3    0 0.3951929
#> 4  4    0 0.2836606
#> 5  5    0 0.3751620
#> 6  6    0 0.3148051
# }
```
