###############################################################################
## test_sample_data.R
##
## Tests for sample_ps_data().
###############################################################################

test_that("sample_ps_data() returns a data frame with correct dimensions", {
  dta <- sample_ps_data(n = 100, seed = 42)
  expect_s3_class(dta, "data.frame")
  expect_equal(nrow(dta), 200L)
})

test_that("sample_ps_data() contains the expected columns", {
  dta <- sample_ps_data(n = 50, seed = 1)
  expected_cols <- c("id", "tavr", "age", "female", "ef",
                     "diabetes", "hypertension", "prob_t", "match")
  expect_true(all(expected_cols %in% names(dta)))
})

test_that("tavr column is binary 0/1 with equal group sizes", {
  dta <- sample_ps_data(n = 200, seed = 7)
  expect_true(all(dta$tavr %in% c(0L, 1L)))
  expect_equal(sum(dta$tavr == 0L), 200L)
  expect_equal(sum(dta$tavr == 1L), 200L)
})

test_that("prob_t is numeric and in [0, 1]", {
  dta <- sample_ps_data(n = 100, seed = 99)
  expect_true(is.numeric(dta$prob_t))
  expect_true(all(dta$prob_t >= 0 & dta$prob_t <= 1))
})

test_that("match column is initialised to 0", {
  dta <- sample_ps_data(n = 50, seed = 5)
  expect_true(all(dta$match == 0L))
})

test_that("sample_ps_data() is reproducible with the same seed", {
  d1 <- sample_ps_data(n = 100, seed = 42)
  d2 <- sample_ps_data(n = 100, seed = 42)
  expect_equal(d1, d2)
})

test_that("sample_ps_data() differs across seeds", {
  d1 <- sample_ps_data(n = 100, seed = 1)
  d2 <- sample_ps_data(n = 100, seed = 2)
  expect_false(identical(d1$prob_t, d2$prob_t))
})

test_that("sample_ps_data() errors on invalid n", {
  expect_error(sample_ps_data(n = 0))
  expect_error(sample_ps_data(n = -10))
  expect_error(sample_ps_data(n = "abc"))
})

test_that("sample_ps_data() errors on invalid separation", {
  expect_error(sample_ps_data(separation = -1))
  expect_error(sample_ps_data(separation = "x"))
})

test_that("separation parameter affects group overlap", {
  dta_lo <- sample_ps_data(n = 200, seed = 42, separation = 0.1)
  dta_hi <- sample_ps_data(n = 200, seed = 42, separation = 3.0)
  # Higher separation → more extreme propensity scores (closer to 0 or 1)
  ps_hi <- dta_hi$prob_t
  ps_lo <- dta_lo$prob_t
  # SD of extreme scores should be lower (pushed to tails)
  expect_true(stats::sd(ps_hi) < stats::sd(ps_lo) ||
                mean(abs(ps_hi - 0.5)) > mean(abs(ps_lo - 0.5)))
})
