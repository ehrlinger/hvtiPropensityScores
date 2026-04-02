###############################################################################
## sa-sensitivity.R
##
## Sensitivity analysis functions for propensity score analyses.
##
## Four complementary methods:
##
##   sa_rosenbaum()  -- Rosenbaum gamma bounds for matched analyses
##                     (requires rbounds)
##   sa_evalue()     -- E-values (VanderWeele & Ding 2017); no extra dependency
##   sa_overlap()    -- Overlap / positivity diagnostics on the PS distribution
##   sa_trim_sweep() -- IPTW trim-threshold sensitivity sweep
##
###############################################################################


# ---------------------------------------------------------------------------
# sa_rosenbaum -- Rosenbaum sensitivity bounds
# ---------------------------------------------------------------------------

#' Rosenbaum Sensitivity Bounds for Matched Analyses
#'
#' Assesses how sensitive a matched-pair analysis is to unmeasured confounding
#' using the Rosenbaum bounds framework.  For each value of the sensitivity
#' parameter Gamma -- the maximum odds ratio for differential treatment
#' assignment within a matched pair due to hidden bias -- the function reports
#' the upper and lower bounds on the p-value of the Wilcoxon signed-rank test.
#'
#' A Gamma of 1 corresponds to no hidden bias (equivalent to perfect
#' randomisation within pairs).  The *sensitivity value* is the largest Gamma
#' at which the upper p-value bound remains below `alpha`; a larger sensitivity
#' value indicates a more robust result.
#'
#' @param x             A `ps_match` object (returned by [ps_match()]) **or** a
#'   plain data frame whose rows are matched patients.
#' @param outcome_col   Name of the numeric outcome column.  Required.
#' @param gamma_max     Maximum Gamma to test.  Default `3`.
#' @param gamma_inc     Increment between Gamma values.  Default `0.25`.
#' @param alpha         Significance threshold.  Default `0.05`.
#' @param treatment_col Name of the binary treatment column (0/1).  Taken from
#'   `x$meta$treatment_col` when `x` is a `ps_match` object.  Required when
#'   `x` is a plain data frame.
#' @param pair_id_col   Name of the pair-ID column created by [ps_match()].
#'   Default `"pair_id"`.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{`$bounds`}{Data frame with columns `gamma`, `p_upper`, `p_lower`,
#'     `reject_upper` (logical -- is the upper bound < `alpha`?).}
#'   \item{`$sensitivity_value`}{Largest Gamma for which the analysis remains
#'     significant (upper p-value < `alpha`), or `NA` if the initial analysis
#'     is not significant.}
#'   \item{`$n_pairs`}{Number of complete matched pairs used.}
#' }
#'
#' @seealso [ps_match()], [sa_evalue()], [sa_overlap()], [sa_trim_sweep()]
#'
#' @references
#' Rosenbaum, P. R. (2002). *Observational Studies* (2nd ed.). Springer.
#'
#' @examples
#' dta <- sample_ps_data(n = 300, seed = 1)
#' m   <- ps_match(dta)
#' \dontrun{
#'   # requires rbounds
#'   res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 2)
#'   res$bounds
#'   res$sensitivity_value
#' }
#'
#' @export
sa_rosenbaum <- function(x,
                         outcome_col,
                         gamma_max     = 3,
                         gamma_inc     = 0.25,
                         alpha         = 0.05,
                         treatment_col = NULL,
                         pair_id_col   = "pair_id") {

  if (!requireNamespace("rbounds", quietly = TRUE)) {
    rlang::abort(
      c("Package 'rbounds' is required for sa_rosenbaum().",
        i = "Install it with: install.packages('rbounds')"),
      call. = FALSE
    )
  }

  # ---- Extract data frame --------------------------------------------------
  if (is_ps_data(x)) {
    if (is.null(treatment_col)) treatment_col <- x$meta$treatment_col
    dta <- x$data
  } else if (is.data.frame(x)) {
    if (is.null(treatment_col)) {
      rlang::abort(
        "`treatment_col` must be specified when `x` is a plain data frame.",
        call. = FALSE
      )
    }
    dta <- x
  } else {
    rlang::abort("`x` must be a ps_match object or a data frame.", call. = FALSE)
  }

  # ---- Validate columns ----------------------------------------------------
  .check_cols(dta, c(outcome_col, treatment_col, pair_id_col))

  # ---- Extract matched pairs -----------------------------------------------
  matched <- dta[!is.na(dta[[pair_id_col]]), , drop = FALSE]

  if (nrow(matched) == 0L) {
    rlang::abort(
      "No matched pairs found -- all `pair_id` values are NA.",
      call. = FALSE
    )
  }

  trt_rows <- matched[as.integer(matched[[treatment_col]]) == 1L, , drop = FALSE]
  ctl_rows <- matched[as.integer(matched[[treatment_col]]) == 0L, , drop = FALSE]

  pair_ids <- sort(unique(trt_rows[[pair_id_col]]))
  trt_out  <- trt_rows[[outcome_col]][match(pair_ids, trt_rows[[pair_id_col]])]
  ctl_out  <- ctl_rows[[outcome_col]][match(pair_ids, ctl_rows[[pair_id_col]])]

  # Drop pairs with missing outcome values
  complete <- !is.na(trt_out) & !is.na(ctl_out)
  trt_out  <- trt_out[complete]
  ctl_out  <- ctl_out[complete]
  n_pairs  <- sum(complete)

  if (n_pairs < 2L) {
    rlang::abort(
      "At least 2 complete matched pairs are required.",
      call. = FALSE
    )
  }

  # ---- Wilcoxon signed-rank sensitivity analysis ---------------------------
  pair_mat <- cbind(trt_out, ctl_out)

  res <- rbounds::psens(
    x        = pair_mat,
    Gamma    = gamma_max,
    GammaInc = gamma_inc
  )

  bounds <- data.frame(
    gamma        = res[["Gamma"]],
    p_upper      = res[["pval+"]],
    p_lower      = res[["pval-"]],
    reject_upper = res[["pval+"]] < alpha,
    stringsAsFactors = FALSE
  )

  sig_gammas        <- bounds$gamma[bounds$reject_upper]
  sensitivity_value <- if (length(sig_gammas) > 0L) max(sig_gammas) else NA_real_

  list(
    bounds            = bounds,
    sensitivity_value = sensitivity_value,
    n_pairs           = n_pairs
  )
}


# ---------------------------------------------------------------------------
# sa_evalue -- E-value calculation
# ---------------------------------------------------------------------------

#' E-Value Sensitivity Analysis
#'
#' Computes the E-value for a point estimate and optionally for a confidence
#' interval bound, following VanderWeele & Ding (2017).  The E-value is the
#' minimum strength of association (on the risk-ratio scale) that an unmeasured
#' confounder would need with **both** the treatment and the outcome,
#' conditional on measured covariates, to fully explain away the observed
#' association.
#'
#' For odds ratios and hazard ratios the rare-outcome approximation
#' `RR_approx = sqrt(OR)` / `sqrt(HR)` is used to convert to the risk-ratio
#' scale before applying the E-value formula.  For risk differences a baseline
#' outcome probability `p0` must be supplied.
#'
#' @param estimate  Numeric.  The point estimate (RR, OR, HR, or RD on the
#'   original scale).  Must be > 0 for RR/OR/HR; any finite value for RD.
#' @param ci_lo     Optional numeric.  Lower confidence interval bound.
#' @param ci_hi     Optional numeric.  Upper confidence interval bound.
#'   Both `ci_lo` and `ci_hi` must be supplied together to compute the CI
#'   E-value.
#' @param type      One of `"RR"` (default), `"OR"`, `"HR"`, or `"RD"`.
#' @param p0        Baseline outcome probability in the unexposed group.
#'   Required when `type = "RD"`.  Must be in (0, 1).
#'
#' @return A list with:
#' \describe{
#'   \item{`$evalue_estimate`}{E-value for the point estimate.}
#'   \item{`$evalue_ci`}{E-value for the confidence interval bound closest to
#'     the null.  `NA` if `ci_lo`/`ci_hi` are not supplied.  Returns `1` if
#'     the CI already crosses the null.}
#'   \item{`$rr_estimate`}{The point estimate converted to the risk-ratio
#'     scale (same as `estimate` when `type = "RR"`).}
#'   \item{`$type`}{The measure type supplied.}
#' }
#'
#' @references
#' VanderWeele, T. J., & Ding, P. (2017).  Sensitivity analysis in
#' observational research: Introducing the E-value.  *Annals of Internal
#' Medicine*, 167(4), 268-274.  \doi{10.7326/M16-2607}
#'
#' @seealso [sa_rosenbaum()], [sa_overlap()], [sa_trim_sweep()]
#'
#' @examples
#' # Risk ratio of 2.5 -- how robust is the result?
#' sa_evalue(2.5, ci_lo = 1.8, ci_hi = 3.5)
#'
#' # Odds ratio of 0.6 (protective effect)
#' sa_evalue(0.6, ci_lo = 0.4, ci_hi = 0.9, type = "OR")
#'
#' # Risk difference of 0.10 with baseline probability 0.20
#' sa_evalue(0.10, ci_lo = 0.02, ci_hi = 0.18, type = "RD", p0 = 0.20)
#'
#' @export
sa_evalue <- function(estimate,
                      ci_lo    = NULL,
                      ci_hi    = NULL,
                      type     = c("RR", "OR", "HR", "RD"),
                      p0       = NULL) {

  type <- match.arg(type)

  if (!is.numeric(estimate) || length(estimate) != 1L || is.na(estimate)) {
    rlang::abort("`estimate` must be a single, non-missing numeric value.",
                 call. = FALSE)
  }

  if (type %in% c("RR", "OR", "HR") && estimate <= 0) {
    rlang::abort(
      sprintf("`estimate` must be positive for type = '%s'.", type),
      call. = FALSE
    )
  }

  if (type == "RD" && is.null(p0)) {
    rlang::abort(
      "`p0` (baseline outcome probability) must be supplied when type = 'RD'.",
      call. = FALSE
    )
  }

  if (!is.null(p0)) {
    if (!is.numeric(p0) || length(p0) != 1L || p0 <= 0 || p0 >= 1) {
      rlang::abort("`p0` must be a single number in (0, 1).", call. = FALSE)
    }
  }

  if (xor(is.null(ci_lo), is.null(ci_hi))) {
    rlang::abort(
      "Both `ci_lo` and `ci_hi` must be supplied together, or both omitted.",
      call. = FALSE
    )
  }

  # ---- Convert estimate to RR scale ----------------------------------------
  .to_rr <- function(val) {
    switch(type,
      RR  = val,
      OR  = sqrt(val),          # rare-outcome approximation
      HR  = sqrt(val),          # rare-outcome approximation
      RD  = (val + p0) / p0    # RR = p1/p0, p1 = RD + p0
    )
  }

  # ---- E-value formula for RR >= 1 ------------------------------------------
  # For a protective effect (rr < 1) flip to rr > 1 first.
  .evalue_rr <- function(rr) {
    if (is.na(rr) || rr == 1) return(1)
    if (rr < 1) rr <- 1 / rr
    rr + sqrt(rr * (rr - 1))
  }

  rr_est <- .to_rr(estimate)
  ev_est <- .evalue_rr(rr_est)

  # ---- E-value for CI bound ------------------------------------------------
  ev_ci <- NA_real_
  if (!is.null(ci_lo) && !is.null(ci_hi)) {
    null_val <- if (type == "RD") 0 else 1

    # Pick the CI bound closer to the null
    ci_bound <- if (estimate > null_val) ci_lo else ci_hi

    # If CI already crosses or touches null the E-value is 1 by convention
    null_crossed <- (estimate > null_val && ci_lo <= null_val) ||
                    (estimate < null_val && ci_hi >= null_val) ||
                    estimate == null_val

    if (null_crossed) {
      ev_ci <- 1
    } else {
      ev_ci <- .evalue_rr(.to_rr(ci_bound))
    }
  }

  list(
    evalue_estimate = ev_est,
    evalue_ci       = ev_ci,
    rr_estimate     = rr_est,
    type            = type
  )
}


# ---------------------------------------------------------------------------
# sa_overlap -- overlap / positivity diagnostics
# ---------------------------------------------------------------------------

#' Propensity Score Overlap and Positivity Diagnostics
#'
#' Examines the propensity score (PS) distribution for potential positivity
#' violations -- regions of covariate space where the probability of receiving
#' one treatment is near zero -- and quantifies the common support region.
#'
#' The function reports:
#' - PS summary statistics by treatment group.
#' - The empirical overlap region: \[max(min_treated, min_control),
#'   min(max_treated, max_control)\].
#' - The proportion of each group whose PS lies outside the overlap region.
#' - The number and proportion of patients with PS below `trim_threshold` or
#'   above `1 - trim_threshold` (near-positivity-violation flags).
#'
#' @param x             A `ps_data` object **or** a plain data frame.
#' @param score_col     Name of the propensity score column.  Taken from
#'   `x$meta$score_col` when `x` is a `ps_data` object.  Default `"prob_t"`.
#' @param treatment_col Name of the binary treatment column (0/1).  Taken from
#'   `x$meta$treatment_col` when `x` is a `ps_data` object`.  Default `"tavr"`.
#' @param trim_threshold Numeric in (0, 0.5).  PS values below this threshold
#'   or above `1 - trim_threshold` are flagged as near-positivity violations.
#'   Default `0.05`.
#'
#' @return A list with:
#' \describe{
#'   \item{`$summary`}{Data frame: one row per treatment group with `group`,
#'     `n`, `ps_min`, `ps_q25`, `ps_median`, `ps_mean`, `ps_q75`, `ps_max`.}
#'   \item{`$overlap_region`}{Named numeric vector: `lower` and `upper` bounds
#'     of the empirical common support region.}
#'   \item{`$outside_overlap`}{Data frame: `group`, `n_outside`,
#'     `pct_outside` -- patients whose PS is outside the overlap region.}
#'   \item{`$positivity_flags`}{Data frame: `group`, `n_near_zero`
#'     (PS < `trim_threshold`), `n_near_one` (PS > `1 - trim_threshold`),
#'     `pct_near_zero`, `pct_near_one`.}
#' }
#'
#' @seealso [ps_match()], [ps_weight()], [sa_rosenbaum()], [sa_evalue()],
#'   [sa_trim_sweep()]
#'
#' @examples
#' dta <- sample_ps_data(n = 400, seed = 7)
#' res <- sa_overlap(dta, score_col = "prob_t", treatment_col = "tavr")
#' res$overlap_region
#' res$positivity_flags
#'
#' # From a ps_match object
#' m   <- ps_match(dta)
#' res <- sa_overlap(m)
#' res$outside_overlap
#'
#' @export
sa_overlap <- function(x,
                       score_col     = "prob_t",
                       treatment_col = "tavr",
                       trim_threshold = 0.05) {

  # ---- Extract data frame --------------------------------------------------
  if (is_ps_data(x)) {
    score_col     <- rlang::`%||%`(x$meta$score_col,     score_col)
    treatment_col <- rlang::`%||%`(x$meta$treatment_col, treatment_col)
    dta <- x$data
  } else if (is.data.frame(x)) {
    dta <- x
  } else {
    rlang::abort("`x` must be a ps_data object or a data frame.", call. = FALSE)
  }

  # ---- Validate inputs -----------------------------------------------------
  .check_cols(dta, c(score_col, treatment_col))
  .check_probability(dta, score_col)

  if (!is.numeric(trim_threshold) || length(trim_threshold) != 1L ||
      trim_threshold <= 0 || trim_threshold >= 0.5) {
    rlang::abort("`trim_threshold` must be a single number in (0, 0.5).",
                 call. = FALSE)
  }

  ps  <- dta[[score_col]]
  trt <- as.integer(dta[[treatment_col]])

  ps0 <- ps[trt == 0L]
  ps1 <- ps[trt == 1L]
  n0  <- length(ps0)
  n1  <- length(ps1)

  # ---- Summary statistics --------------------------------------------------
  .ps_summary <- function(psv, grp, n) {
    data.frame(
      group     = grp,
      n         = n,
      ps_min    = round(min(psv,    na.rm = TRUE), 4),
      ps_q25    = round(stats::quantile(psv, 0.25, na.rm = TRUE), 4),
      ps_median = round(stats::median(psv,  na.rm = TRUE), 4),
      ps_mean   = round(mean(psv,   na.rm = TRUE), 4),
      ps_q75    = round(stats::quantile(psv, 0.75, na.rm = TRUE), 4),
      ps_max    = round(max(psv,    na.rm = TRUE), 4),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
  summary_tbl <- rbind(.ps_summary(ps0, "control", n0),
                       .ps_summary(ps1, "treated", n1))

  # ---- Overlap region ------------------------------------------------------
  overlap_lo <- max(min(ps0, na.rm = TRUE), min(ps1, na.rm = TRUE))
  overlap_hi <- min(max(ps0, na.rm = TRUE), max(ps1, na.rm = TRUE))
  overlap_region <- c(lower = round(overlap_lo, 4), upper = round(overlap_hi, 4))

  # ---- Outside overlap -----------------------------------------------------
  .n_outside <- function(psv) sum(psv < overlap_lo | psv > overlap_hi, na.rm = TRUE)
  n_out0 <- .n_outside(ps0)
  n_out1 <- .n_outside(ps1)

  outside_overlap <- data.frame(
    group      = c("control", "treated"),
    n_outside  = c(n_out0, n_out1),
    pct_outside = round(100 * c(n_out0 / n0, n_out1 / n1), 1),
    stringsAsFactors = FALSE
  )

  # ---- Positivity flags ----------------------------------------------------
  nz0 <- sum(ps0 < trim_threshold, na.rm = TRUE)
  no0 <- sum(ps0 > 1 - trim_threshold, na.rm = TRUE)
  nz1 <- sum(ps1 < trim_threshold, na.rm = TRUE)
  no1 <- sum(ps1 > 1 - trim_threshold, na.rm = TRUE)

  positivity_flags <- data.frame(
    group        = c("control", "treated"),
    n_near_zero  = c(nz0, nz1),
    n_near_one   = c(no0, no1),
    pct_near_zero = round(100 * c(nz0 / n0, nz1 / n1), 1),
    pct_near_one  = round(100 * c(no0 / n0, no1 / n1), 1),
    stringsAsFactors = FALSE
  )

  list(
    summary          = summary_tbl,
    overlap_region   = overlap_region,
    outside_overlap  = outside_overlap,
    positivity_flags = positivity_flags
  )
}


# ---------------------------------------------------------------------------
# sa_trim_sweep -- IPTW trim threshold sensitivity sweep
# ---------------------------------------------------------------------------

#' IPTW Weight-Trimming Sensitivity Sweep
#'
#' Re-computes IPTW weights across a range of winsorisation thresholds and
#' reports how the weight distribution and effective sample size (ESS) change
#' as trimming increases.  This helps identify a trimming threshold that
#' adequately controls extreme weights without an unacceptable loss of effective
#' sample size.
#'
#' @param x           A `ps_weight` object (returned by [ps_weight()]) **or**
#'   a `ps_match` object **or** a plain data frame that contains a propensity
#'   score column and a binary treatment column.
#' @param trim_seq    Numeric vector of trimming values to evaluate.  Each
#'   value is the lower tail quantile; weights outside the
#'   (`trim`, `1 - trim`) range are winsorised.  `0` means no trimming.
#'   Default `seq(0, 0.10, by = 0.01)`.
#' @param estimand    One of `"ATE"`, `"ATT"`, or `"ATC"` -- the causal
#'   estimand used to compute IPTW weights.  Taken from `x$meta$estimand`
#'   when `x` is a `ps_weight` object; defaults to `"ATE"` otherwise.
#' @param stabilise   Logical.  Whether to stabilise weights.  Taken from
#'   `x$meta$stabilised` when available; defaults to `TRUE`.
#' @param score_col   Name of the propensity score column.  Taken from
#'   `x$meta$score_col` when `x` is a `ps_data` object.  Default `"prob_t"`.
#' @param treatment_col Name of the binary treatment column.  Taken from
#'   `x$meta$treatment_col` when `x` is a `ps_data` object`.  Default
#'   `"tavr"`.
#'
#' @return A data frame with one row per trim value and columns:
#' \describe{
#'   \item{`trim`}{The winsorisation threshold.}
#'   \item{`ess_control`}{Effective sample size for the control group.}
#'   \item{`ess_treated`}{Effective sample size for the treated group.}
#'   \item{`max_weight`}{Maximum weight across all patients after trimming.}
#'   \item{`sd_weight`}{Standard deviation of weights after trimming.}
#'   \item{`pct_trimmed`}{Percentage of patients whose weight was winsorised.}
#' }
#'
#' @seealso [ps_weight()], [sa_rosenbaum()], [sa_evalue()], [sa_overlap()]
#'
#' @examples
#' dta    <- sample_ps_data(n = 400, seed = 3)
#' w_obj  <- ps_weight(dta, estimand = "ATE")
#' sweep  <- sa_trim_sweep(w_obj)
#' head(sweep)
#'
#' # Directly from a data frame
#' sweep2 <- sa_trim_sweep(dta, score_col = "prob_t", treatment_col = "tavr",
#'                         estimand = "ATE")
#'
#' @export
sa_trim_sweep <- function(x,
                          trim_seq      = seq(0, 0.10, by = 0.01),
                          estimand      = c("ATE", "ATT", "ATC"),
                          stabilise     = TRUE,
                          score_col     = "prob_t",
                          treatment_col = "tavr") {

  # ---- Extract data frame and metadata -------------------------------------
  if (is_ps_data(x)) {
    score_col     <- rlang::`%||%`(x$meta$score_col,     score_col)
    treatment_col <- rlang::`%||%`(x$meta$treatment_col, treatment_col)
    estimand_use  <- rlang::`%||%`(x$meta$estimand,      match.arg(estimand))
    stabilise     <- rlang::`%||%`(x$meta$stabilised,    stabilise)
    dta <- x$data
  } else if (is.data.frame(x)) {
    estimand_use <- match.arg(estimand)
    dta <- x
  } else {
    rlang::abort("`x` must be a ps_data object or a data frame.", call. = FALSE)
  }

  estimand_use <- match.arg(estimand_use, c("ATE", "ATT", "ATC"))

  # ---- Validate inputs -----------------------------------------------------
  .check_cols(dta, c(score_col, treatment_col))
  .check_probability(dta, score_col)

  if (!is.numeric(trim_seq) || any(trim_seq < 0) || any(trim_seq >= 0.5)) {
    rlang::abort(
      "`trim_seq` must be a numeric vector with all values in [0, 0.5).",
      call. = FALSE
    )
  }

  ps  <- dta[[score_col]]
  trt <- as.integer(dta[[treatment_col]])
  n   <- length(ps)

  # ---- Compute raw IPTW weights once ---------------------------------------
  raw_w <- switch(estimand_use,
    ATE = ifelse(trt == 1L, 1 / ps, 1 / (1 - ps)),
    ATT = ifelse(trt == 1L, 1,      ps / (1 - ps)),
    ATC = ifelse(trt == 1L, (1 - ps) / ps, 1)
  )

  if (stabilise) {
    p_trt  <- mean(trt == 1L)
    raw_w  <- raw_w * ifelse(trt == 1L, p_trt, 1 - p_trt)
  }

  # ---- Sweep ---------------------------------------------------------------
  rows <- lapply(sort(trim_seq), function(tr) {
    if (tr == 0) {
      w <- raw_w
      n_trimmed <- 0L
    } else {
      q_lo <- stats::quantile(raw_w, tr,       na.rm = TRUE)
      q_hi <- stats::quantile(raw_w, 1 - tr,   na.rm = TRUE)
      w    <- pmin(pmax(raw_w, q_lo), q_hi)
      n_trimmed <- sum(raw_w < q_lo | raw_w > q_hi, na.rm = TRUE)
    }

    data.frame(
      trim        = tr,
      ess_control = .effective_n(w[trt == 0L]),
      ess_treated = .effective_n(w[trt == 1L]),
      max_weight  = round(max(w, na.rm = TRUE), 4),
      sd_weight   = round(stats::sd(w, na.rm = TRUE), 4),
      pct_trimmed = round(100 * n_trimmed / n, 1),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

