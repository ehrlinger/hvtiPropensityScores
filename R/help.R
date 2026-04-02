###############################################################################
## help.R
##
## Package-level documentation for hvtiPropensityScores.
###############################################################################

#' hvtiPropensityScores: Propensity Score Methods for Cardiac Surgery Studies
#'
#' @description
#' Provides functions for propensity score analysis used in cardiac surgery
#' comparative-effectiveness research.  Supports two balancing strategies:
#'
#' - **Nearest-neighbour matching** (`ps_match`): 1:1 greedy matching without
#'   replacement on the propensity score, with optional caliper.
#' - **Inverse-probability-of-treatment weighting** (`ps_weight`): IPTW weights
#'   for the ATE, ATT, or ATC estimand, with optional stabilisation and
#'   winsorisation.
#'
#' Both functions return a `ps_data` S3 object whose `$data` slot holds the
#' original dataset with the match indicator or weight column appended, and
#' whose `$tables` slot holds standardised mean difference (SMD) diagnostics
#' compatible with [hvtiPlotR::hv_mirror_hist()] and
#' [hvtiPlotR::hv_balance()].
#'
#' @section Typical workflow:
#' ```r
#' library(hvtiPropensityScores)
#'
#' # Generate or load your data
#' dta <- sample_ps_data(n = 500, seed = 42)
#'
#' # --- Matching ---
#' m   <- ps_match(dta)
#' print(m)
#' summary(m)
#' matched <- m$data[m$data$match == 1L, ]
#'
#' # --- Weighting ---
#' w   <- ps_weight(dta, estimand = "ATE")
#' print(w)
#' summary(w)
#' ```
#'
#' @seealso
#' - [ps_match()] — nearest-neighbour 1:1 matching
#' - [ps_weight()] — IPTW weighting (ATE / ATT / ATC)
#' - [sample_ps_data()] — reproducible synthetic dataset for examples / tests
#' - [is_ps_data()] — predicate for `ps_data` objects
#'
#' @name hvtiPropensityScores-package
#' @aliases hvtiPropensityScores
#' @keywords internal
"_PACKAGE"
