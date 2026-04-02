###############################################################################
## sample-data.R
##
## Reproducible synthetic dataset for examples and tests.
##
## `sample_ps_data()` generates a data frame that mimics the structure of
## a typical cardiac-surgery SAS export used for propensity score analyses:
##   - binary treatment indicator (`tavr`)
##   - true propensity score (`prob_t`)  computed from a logistic model
##   - several continuous and binary covariates
##   - a match indicator column (`match`) populated with 0s (filled in by
##     `ps_match()`)
###############################################################################


#' Generate sample data for propensity score examples
#'
#' Produces a synthetic data frame that resembles a typical SAS export from a
#' cardiac-surgery propensity score analysis.  Propensity scores are generated
#' from a logistic model driven by the supplied covariates, so the resulting
#' `prob_t` column has realistic overlap and the SMD values are non-trivial.
#'
#' @param n         Number of patients per treatment arm.  Total rows =
#'   `2 * n`.  Must be a positive integer.  Default `500`.
#' @param seed      Random seed for reproducibility.  Default `42L`.
#' @param separation Numeric scalar controlling how well the covariates
#'   separate the two groups.  Higher values → less overlap.  Default `1.0`.
#'
#' @return A data frame with `2 * n` rows and the following columns:
#'   \describe{
#'     \item{`id`}{Integer patient identifier.}
#'     \item{`tavr`}{Binary treatment indicator: `0` = SAVR, `1` = TF-TAVR.}
#'     \item{`age`}{Continuous: age in years.}
#'     \item{`female`}{Binary: `1` = female.}
#'     \item{`ef`}{Continuous: ejection fraction (%).}
#'     \item{`diabetes`}{Binary: `1` = diabetic.}
#'     \item{`hypertension`}{Binary: `1` = hypertensive.}
#'     \item{`prob_t`}{Estimated propensity score (numeric, 0–1).}
#'     \item{`match`}{Match indicator, initialised to `0`.  Populated by
#'       [ps_match()].}
#'   }
#'
#' @seealso [ps_match()], [ps_weight()]
#'
#' @examples
#' dta <- sample_ps_data(n = 200, seed = 42)
#' head(dta)
#' table(dta$tavr)
#' hist(dta$prob_t, breaks = 20, main = "Propensity score distribution")
#'
#' @export
sample_ps_data <- function(n = 500L, seed = 42L, separation = 1.0) {
  if (!is.numeric(n) || length(n) != 1L || n < 1L)
    stop("`n` must be a positive integer.", call. = FALSE)
  if (!is.numeric(separation) || length(separation) != 1L || separation <= 0)
    stop("`separation` must be a positive number.", call. = FALSE)

  set.seed(seed)
  n <- as.integer(n)

  # Control arm (tavr = 0): older, lower EF
  age0      <- rnorm(n, mean = 75 + separation * 3, sd = 8)
  female0   <- stats::rbinom(n, 1, 0.40)
  ef0       <- rnorm(n, mean = 55 - separation * 3, sd = 10)
  diabetes0 <- stats::rbinom(n, 1, 0.30)
  hyp0      <- stats::rbinom(n, 1, 0.65)

  # Treatment arm (tavr = 1): younger, higher EF
  age1      <- rnorm(n, mean = 75 - separation * 3, sd = 8)
  female1   <- stats::rbinom(n, 1, 0.45)
  ef1       <- rnorm(n, mean = 55 + separation * 3, sd = 10)
  diabetes1 <- stats::rbinom(n, 1, 0.25)
  hyp1      <- stats::rbinom(n, 1, 0.60)

  dta <- data.frame(
    id           = seq_len(2L * n),
    tavr         = c(rep(0L, n), rep(1L, n)),
    age          = c(age0, age1),
    female       = c(female0, female1),
    ef           = c(ef0, ef1),
    diabetes     = c(diabetes0, diabetes1),
    hypertension = c(hyp0, hyp1)
  )

  # Fit a quick logistic model to produce realistic propensity scores
  fit <- stats::glm(
    tavr ~ age + female + ef + diabetes + hypertension,
    data   = dta,
    family = stats::binomial()
  )
  dta$prob_t <- as.numeric(stats::predict(fit, type = "response"))
  dta$match  <- 0L

  dta
}
