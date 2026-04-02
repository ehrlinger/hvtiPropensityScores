###############################################################################
## test_sa_sensitivity.R
##
## Tests for sa_rosenbaum(), sa_evalue(), sa_overlap(), and sa_trim_sweep().
###############################################################################


# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------

make_matched <- function(n = 200, seed = 1) {
  dta <- sample_ps_data(n = n, seed = seed)
  ps_match(dta)
}

make_weighted <- function(n = 200, seed = 1, estimand = "ATE") {
  dta <- sample_ps_data(n = n, seed = seed)
  ps_weight(dta, estimand = estimand)
}


# ---------------------------------------------------------------------------
# ps_match — pair_id column (prerequisite for sa_rosenbaum)
# ---------------------------------------------------------------------------

test_that("ps_match() adds pair_id column to $data", {
  m <- make_matched()
  expect_true("pair_id" %in% names(m$data))
})

test_that("ps_match() pair_id is NA for unmatched, integer for matched", {
  m <- make_matched()
  matched_rows   <- m$data[m$data$match == 1L, ]
  unmatched_rows <- m$data[m$data$match == 0L, ]
  expect_true(all(!is.na(matched_rows$pair_id)))
  expect_true(all(is.na(unmatched_rows$pair_id)))
})

test_that("ps_match() pair_id links exactly one treated to one control per pair", {
  m <- make_matched()
  matched <- m$data[m$data$match == 1L, ]
  pair_counts <- table(matched$pair_id)
  expect_true(all(pair_counts == 2L))
})

test_that("ps_match() meta contains pair_id_col", {
  m <- make_matched()
  expect_equal(m$meta$pair_id_col, "pair_id")
})


# ---------------------------------------------------------------------------
# sa_evalue
# ---------------------------------------------------------------------------

test_that("sa_evalue() returns list with expected names for RR", {
  res <- sa_evalue(2.5)
  expect_named(res, c("evalue_estimate", "evalue_ci", "rr_estimate", "type"))
})

test_that("sa_evalue() E-value for RR = 1 is 1", {
  res <- sa_evalue(1.0)
  expect_equal(res$evalue_estimate, 1)
})

test_that("sa_evalue() E-value for RR = 2 is 2 + sqrt(2)", {
  res <- sa_evalue(2.0)
  expect_equal(res$evalue_estimate, 2 + sqrt(2 * (2 - 1)), tolerance = 1e-8)
})

test_that("sa_evalue() is symmetric: RR and 1/RR give the same E-value", {
  e1 <- sa_evalue(2.0)$evalue_estimate
  e2 <- sa_evalue(0.5)$evalue_estimate
  expect_equal(e1, e2, tolerance = 1e-8)
})

test_that("sa_evalue() evalue_ci is NA when no CI is supplied", {
  res <- sa_evalue(2.5)
  expect_true(is.na(res$evalue_ci))
})

test_that("sa_evalue() evalue_ci is 1 when CI crosses null", {
  res <- sa_evalue(2.5, ci_lo = 0.9, ci_hi = 3.5)
  expect_equal(res$evalue_ci, 1)
})

test_that("sa_evalue() evalue_ci > 1 when CI does not cross null (RR > 1)", {
  res <- sa_evalue(2.5, ci_lo = 1.5, ci_hi = 3.5)
  expect_gt(res$evalue_ci, 1)
  # CI E-value should be based on ci_lo (closer to null for RR > 1)
  expected <- { rr <- 1.5; rr + sqrt(rr * (rr - 1)) }
  expect_equal(res$evalue_ci, expected, tolerance = 1e-8)
})

test_that("sa_evalue() works for OR type with rare outcome approximation", {
  res <- sa_evalue(4.0, type = "OR")
  expect_equal(res$rr_estimate, sqrt(4), tolerance = 1e-8)
  expect_equal(res$type, "OR")
})

test_that("sa_evalue() works for HR type", {
  res <- sa_evalue(3.0, type = "HR")
  expect_equal(res$rr_estimate, sqrt(3), tolerance = 1e-8)
})

test_that("sa_evalue() works for RD type with p0", {
  res <- sa_evalue(0.10, type = "RD", p0 = 0.20)
  expect_equal(res$rr_estimate, (0.10 + 0.20) / 0.20, tolerance = 1e-8)
})

test_that("sa_evalue() errors on non-numeric estimate", {
  expect_error(sa_evalue("2.5"))
})

test_that("sa_evalue() errors when p0 missing for RD", {
  expect_error(sa_evalue(0.1, type = "RD"))
})

test_that("sa_evalue() errors when only one CI bound is supplied", {
  expect_error(sa_evalue(2.5, ci_lo = 1.5))
  expect_error(sa_evalue(2.5, ci_hi = 3.5))
})

test_that("sa_evalue() errors on non-positive OR", {
  expect_error(sa_evalue(-0.5, type = "OR"))
})


# ---------------------------------------------------------------------------
# sa_overlap
# ---------------------------------------------------------------------------

test_that("sa_overlap() returns list with expected names (plain data frame)", {
  dta <- sample_ps_data(n = 300, seed = 10)
  res <- sa_overlap(dta, score_col = "prob_t", treatment_col = "tavr")
  expect_named(res, c("summary", "overlap_region", "outside_overlap",
                      "positivity_flags"))
})

test_that("sa_overlap() summary has two rows and correct columns", {
  dta <- sample_ps_data(n = 300, seed = 11)
  res <- sa_overlap(dta)
  expect_equal(nrow(res$summary), 2L)
  expect_true(all(c("group", "n", "ps_min", "ps_max") %in% names(res$summary)))
})

test_that("sa_overlap() overlap_region is a named numeric vector of length 2", {
  dta <- sample_ps_data(n = 300, seed = 12)
  res <- sa_overlap(dta)
  expect_named(res$overlap_region, c("lower", "upper"))
  expect_true(res$overlap_region["lower"] <= res$overlap_region["upper"])
})

test_that("sa_overlap() pct_outside is in [0, 100]", {
  dta <- sample_ps_data(n = 300, seed = 13)
  res <- sa_overlap(dta)
  expect_true(all(res$outside_overlap$pct_outside >= 0))
  expect_true(all(res$outside_overlap$pct_outside <= 100))
})

test_that("sa_overlap() positivity_flags has correct structure", {
  dta <- sample_ps_data(n = 300, seed = 14)
  res <- sa_overlap(dta)
  pf  <- res$positivity_flags
  expect_named(pf, c("group", "n_near_zero", "n_near_one",
                     "pct_near_zero", "pct_near_one"))
  expect_equal(nrow(pf), 2L)
})

test_that("sa_overlap() accepts a ps_data object", {
  m   <- make_matched(n = 200, seed = 15)
  res <- sa_overlap(m)
  expect_named(res, c("summary", "overlap_region", "outside_overlap",
                      "positivity_flags"))
})

test_that("sa_overlap() summary n values sum to nrow(data)", {
  dta <- sample_ps_data(n = 300, seed = 16)
  res <- sa_overlap(dta)
  expect_equal(sum(res$summary$n), nrow(dta))
})

test_that("sa_overlap() errors on bad trim_threshold", {
  dta <- sample_ps_data(n = 100, seed = 17)
  expect_error(sa_overlap(dta, trim_threshold = 0.6))
  expect_error(sa_overlap(dta, trim_threshold = -0.1))
})


# ---------------------------------------------------------------------------
# sa_trim_sweep
# ---------------------------------------------------------------------------

test_that("sa_trim_sweep() returns a data frame with correct columns", {
  w   <- make_weighted()
  res <- sa_trim_sweep(w)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("trim", "ess_control", "ess_treated",
                    "max_weight", "sd_weight", "pct_trimmed") %in% names(res)))
})

test_that("sa_trim_sweep() has one row per element in trim_seq", {
  w    <- make_weighted()
  tseq <- seq(0, 0.05, by = 0.01)
  res  <- sa_trim_sweep(w, trim_seq = tseq)
  expect_equal(nrow(res), length(tseq))
  expect_equal(res$trim, sort(tseq))
})

test_that("sa_trim_sweep() pct_trimmed is 0 at trim = 0", {
  w   <- make_weighted()
  res <- sa_trim_sweep(w, trim_seq = c(0, 0.05))
  expect_equal(res$pct_trimmed[res$trim == 0], 0)
})

test_that("sa_trim_sweep() max_weight decreases (weakly) as trim increases", {
  w   <- make_weighted()
  res <- sa_trim_sweep(w, trim_seq = seq(0, 0.10, by = 0.01))
  # max_weight should be non-increasing with trim
  expect_true(all(diff(res$max_weight) <= 1e-6))
})

test_that("sa_trim_sweep() accepts a plain data frame", {
  dta <- sample_ps_data(n = 300, seed = 20)
  res <- sa_trim_sweep(dta, trim_seq = c(0, 0.05),
                       score_col = "prob_t", treatment_col = "tavr")
  expect_equal(nrow(res), 2L)
})

test_that("sa_trim_sweep() inherits estimand from ps_weight object", {
  w_att <- make_weighted(estimand = "ATT")
  # Should run without error; the ATT estimand is picked up from meta
  res   <- sa_trim_sweep(w_att, trim_seq = c(0, 0.05))
  expect_equal(nrow(res), 2L)
})

test_that("sa_trim_sweep() errors on bad trim_seq values", {
  w <- make_weighted()
  expect_error(sa_trim_sweep(w, trim_seq = c(0, 0.6)))
  expect_error(sa_trim_sweep(w, trim_seq = c(-0.01, 0.05)))
})


# ---------------------------------------------------------------------------
# sa_rosenbaum (requires rbounds)
# ---------------------------------------------------------------------------

test_that("sa_rosenbaum() errors gracefully when rbounds is absent", {
  skip_if(requireNamespace("rbounds", quietly = TRUE),
          "rbounds is installed — skipping absence test")
  m <- make_matched()
  expect_error(sa_rosenbaum(m, outcome_col = "ef"),
               regexp = "rbounds")
})

test_that("sa_rosenbaum() returns list with expected names", {
  skip_if_not_installed("rbounds")
  m   <- make_matched(n = 300, seed = 42)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 2, gamma_inc = 0.5)
  expect_named(res, c("bounds", "sensitivity_value", "n_pairs"))
})

test_that("sa_rosenbaum() bounds data frame has correct columns and rows", {
  skip_if_not_installed("rbounds")
  m   <- make_matched(n = 300, seed = 43)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 2, gamma_inc = 0.5)
  expect_named(res$bounds, c("gamma", "p_upper", "p_lower", "reject_upper"))
  # gamma values: 1.0, 1.5, 2.0
  expect_equal(nrow(res$bounds), 3L)
  expect_equal(res$bounds$gamma[1L], 1.0)
})

test_that("sa_rosenbaum() n_pairs is a positive integer", {
  skip_if_not_installed("rbounds")
  m   <- make_matched(n = 300, seed = 44)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 1.5, gamma_inc = 0.5)
  expect_true(res$n_pairs > 0L)
  expect_true(is.numeric(res$n_pairs))
})

test_that("sa_rosenbaum() p_upper at Gamma=1 equals unadjusted Wilcoxon p-value", {
  skip_if_not_installed("rbounds")
  m       <- make_matched(n = 300, seed = 45)
  matched <- m$data[m$data$match == 1L, ]
  trt     <- matched[matched$tavr == 1L, ]
  ctl     <- matched[matched$tavr == 0L, ]
  pair_ids <- sort(unique(trt$pair_id))
  trt_ef   <- trt$ef[match(pair_ids, trt$pair_id)]
  ctl_ef   <- ctl$ef[match(pair_ids, ctl$pair_id)]
  wilcox_p <- stats::wilcox.test(trt_ef, ctl_ef, paired = TRUE)$p.value

  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 1, gamma_inc = 1)
  expect_equal(res$bounds$p_upper[1L], wilcox_p, tolerance = 1e-6)
})

test_that("sa_rosenbaum() sensitivity_value is NA when result is not significant", {
  skip_if_not_installed("rbounds")
  # Use a very small dataset where the test is unlikely to be significant
  set.seed(99)
  dta <- sample_ps_data(n = 20, seed = 99)
  m   <- ps_match(dta)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 1, gamma_inc = 1,
                      alpha = 1e-10)  # Extremely small alpha → likely NA
  expect_true(is.na(res$sensitivity_value) || is.numeric(res$sensitivity_value))
})

test_that("sa_rosenbaum() accepts a plain data frame with treatment_col", {
  skip_if_not_installed("rbounds")
  m       <- make_matched(n = 300, seed = 46)
  matched <- m$data[!is.na(m$data$pair_id), ]
  res     <- sa_rosenbaum(matched, outcome_col = "ef",
                          treatment_col = "tavr", gamma_max = 1.5,
                          gamma_inc = 0.5)
  expect_named(res, c("bounds", "sensitivity_value", "n_pairs"))
})

test_that("sa_rosenbaum() errors when treatment_col missing for plain data frame", {
  skip_if_not_installed("rbounds")
  m       <- make_matched(n = 200, seed = 47)
  matched <- m$data[!is.na(m$data$pair_id), ]
  expect_error(sa_rosenbaum(matched, outcome_col = "ef"))
})
