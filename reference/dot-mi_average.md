# Run a fitting function over each imputation and return averaged predictions

Run a fitting function over each imputation and return averaged
predictions

## Usage

``` r
.mi_average(data, imputation_col, id_col, fit_fn, pred_fn)
```

## Arguments

- data:

  Stacked MI data frame.

- imputation_col:

  Column name holding the imputation index.

- id_col:

  Column name holding the patient ID.

- fit_fn:

  Function (sub_df) -\> fit object.

- pred_fn:

  Function (fit, sub_df) -\> numeric vector / matrix of predicted
  values. Must return results in the same row order as sub_df.

## Value

A list with `preds` (averaged predictions, same type as pred_fn output
for one dataset), `base_data` (first imputation, imputation col
removed), and `n_imp`.
