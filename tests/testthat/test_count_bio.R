# tests/testthat/test_count_bio.R

test_that("count_bio basic grouped counts match dplyr::count_tally()", {
  set.seed(2000)

  df <- dplyr::tibble(
    arm = factor(sample(c("Drug", "Placebo"), 200, TRUE)),
    py  = rexp(200, rate = 1 / 0.5),
    ae  = rbinom(200, 1, 0.25)
  )

  dplyr_res <- df %>%
    dplyr::count(arm) %>%
    dplyr::arrange(arm)

  # mycount, no person-time
  bio_res <- df %>%
    count_bio(arm) %>%
    dplyr::arrange(arm)

  # group labels are same
  expect_equal(bio_res$arm, dplyr_res$arm)
  expect_equal(bio_res$n, dplyr_res$n)
})

test_that("count_bio incidence rates and CI match manual calculations", {
  set.seed(2000)

  df <- dplyr::tibble(
    arm = factor(sample(c("Drug", "Placebo"), 200, TRUE)),
    py  = rexp(200, rate = 1 / 0.5),
    ae  = rbinom(200, 1, 0.25)
  )

  # calculate by count_bio
  bio_rates <- df %>%
    count_bio(arm, events = ae, person_time = py, per = 100) %>%
    dplyr::arrange(arm)

  # calculate by hand
  manual <- df %>%
    dplyr::group_by(arm) %>%
    dplyr::summarise(
      n  = sum(ae, na.rm = TRUE),
      pt = sum(py, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(arm)

  alpha <- 1 - 0.95
  k     <- manual$n
  Tt    <- manual$pt

  rate      <- ifelse(Tt > 0, k / Tt, NA_real_)
  lower_cnt <- 0.5 * stats::qchisq(alpha / 2,     2 * k)
  upper_cnt <- 0.5 * stats::qchisq(1 - alpha / 2, 2 * (k + 1))

  rate_lcl <- ifelse(Tt > 0, lower_cnt / Tt, NA_real_)
  rate_ucl <- ifelse(Tt > 0, upper_cnt / Tt, NA_real_)

  manual$rate     <- rate     * 100
  manual$rate_lcl <- rate_lcl * 100
  manual$rate_ucl <- rate_ucl * 100

  # check equal or not
  expect_equal(bio_rates$n,  manual$n)
  expect_equal(bio_rates$pt, manual$pt)
  expect_equal(bio_rates$rate,     manual$rate,     tolerance = 1e-10)
  expect_equal(bio_rates$rate_lcl, manual$rate_lcl, tolerance = 1e-10)
  expect_equal(bio_rates$rate_ucl, manual$rate_ucl, tolerance = 1e-10)
})

test_that("count_bio errors when person_time is negative", {
  df_bad <- dplyr::tibble(
    arm = factor(c("Drug", "Drug")),
    py  = c(1, -2),
    ae  = c(0, 1)
  )

  expect_error(
    count_bio(df_bad, arm, events = ae, person_time = py),
    "negative values"
  )
})
