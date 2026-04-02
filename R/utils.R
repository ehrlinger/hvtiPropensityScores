###############################################################################
## utils.R
##
## Internal validation helpers shared across ps_*() constructors.
##
## None of these are exported; they are called as hvtiPropensityScores:::fn()
## in tests that need to exercise them directly.
###############################################################################


# ---------------------------------------------------------------------------
# Column-existence checks
# ---------------------------------------------------------------------------

#' Assert that required columns exist in a data frame
#'
#' @param data      A data frame.
#' @param cols      Character vector of required column names.
#' @param call_env  Environment for the error call (use `rlang::caller_env()`).
#' @return Invisible NULL.  Signals a classed error on failure.
#' @keywords internal
.check_cols <- function(data, cols, call_env = rlang::caller_env()) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0L) {
    rlang::abort(
      sprintf(
        "Required column(s) not found in `data`: %s",
        paste(missing, collapse = ", ")
      ),
      call = call_env
    )
  }
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Data-frame check
# ---------------------------------------------------------------------------

#' Assert that an object is a non-empty data frame
#'
#' @param data      Object to check.
#' @param call_env  Environment for the error call.
#' @return Invisible NULL.
#' @keywords internal
.check_df <- function(data, call_env = rlang::caller_env()) {
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.", call = call_env)
  }
  if (nrow(data) == 0L) {
    rlang::abort("`data` must contain at least one row.", call = call_env)
  }
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Binary-column check
# ---------------------------------------------------------------------------

#' Assert that a column contains only 0/1 (or logical) values
#'
#' @param data      A data frame.
#' @param col       Name of the column to check.
#' @param call_env  Environment for the error call.
#' @return Invisible NULL.
#' @keywords internal
.check_binary <- function(data, col, call_env = rlang::caller_env()) {
  x <- data[[col]]
  ok <- is.logical(x) ||
    (is.numeric(x) && all(x %in% c(0, 1, NA_real_)))
  if (!ok) {
    rlang::abort(
      sprintf("Column `%s` must be binary (0/1 or logical).", col),
      call = call_env
    )
  }
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Probability-range check
# ---------------------------------------------------------------------------

#' Assert that a column contains valid probabilities (0 <= x <= 1)
#'
#' @param data      A data frame.
#' @param col       Name of the column to check.
#' @param call_env  Environment for the error call.
#' @return Invisible NULL.
#' @keywords internal
.check_probability <- function(data, col, call_env = rlang::caller_env()) {
  x <- data[[col]]
  if (!is.numeric(x)) {
    rlang::abort(
      sprintf("Column `%s` must be numeric.", col),
      call = call_env
    )
  }
  if (any(!is.na(x) & (x < 0 | x > 1))) {
    rlang::abort(
      sprintf("Column `%s` must contain values in [0, 1].", col),
      call = call_env
    )
  }
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Standardised mean difference (SMD)
# ---------------------------------------------------------------------------

#' Compute standardised mean difference between two numeric vectors
#'
#' Uses the pooled standard deviation: `SMD = (mean1 - mean0) / sd_pooled`.
#' Returns `NA_real_` if either group is empty or has zero variance.
#'
#' @param x         Numeric vector of covariate values.
#' @param treatment Binary vector (0/1 or logical) of treatment indicators.
#' @return A single numeric SMD value.
#' @keywords internal
.calc_smd <- function(x, treatment) {
  x0 <- x[treatment == 0 & !is.na(x)]
  x1 <- x[treatment == 1 & !is.na(x)]

  if (length(x0) < 2L || length(x1) < 2L) return(NA_real_)

  m0 <- mean(x0);  s0 <- stats::sd(x0)
  m1 <- mean(x1);  s1 <- stats::sd(x1)

  sp <- sqrt((s0^2 + s1^2) / 2)
  if (sp == 0) return(NA_real_)

  (m1 - m0) / sp
}


# ---------------------------------------------------------------------------
# SMD table for a data frame
# ---------------------------------------------------------------------------

#' Build a tidy SMD table across all numeric covariates
#'
#' @param data       A data frame.
#' @param treatment  Name of the binary treatment column.
#' @param covariates Character vector of covariate column names.  If `NULL`,
#'   all numeric columns other than `treatment` are used.
#' @return A data frame with columns `variable` and `smd`.
#' @keywords internal
.smd_table <- function(data, treatment, covariates = NULL) {
  if (is.null(covariates)) {
    covariates <- setdiff(
      names(data)[vapply(data, is.numeric, logical(1))],
      treatment
    )
  }

  smds <- vapply(covariates, function(cn) {
    .calc_smd(data[[cn]], data[[treatment]])
  }, numeric(1))

  data.frame(
    variable = covariates,
    smd      = round(smds, 4L),
    stringsAsFactors = FALSE
  )
}
