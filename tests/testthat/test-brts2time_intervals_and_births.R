context("brts2time_intervals_and_births")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

test_that("approximate_brts", {
  brts <- c(6, 5, 4, 3, 1e-14)
  brts_precision <- 6
  approx_brts <- mbd::approximate_brts(
    brts = brts,
    brts_precision = brts_precision
  )
  testthat::expect_true(all(approx_brts > 0))
})

test_that("basic use", {

  brts <- c(10, 9, 8, 7, 7, 6, 6, 6, 5, 5, 5, 5)
  out <- mbd::brts2time_intervals_and_births(brts)
  time_intervals <- out$time_intervals
  births <- out$births

  testthat::expect_equal(
    time_intervals,
    c(0, abs(diff(unique(c(brts, 0)))))
  )
  n_mbd_events <- length(unique(brts[duplicated(brts)]))
  for (i in 1:n_mbd_events) {
    testthat::expect_equal(
      births[births != 1][i + 1],
      sum(unique(brts[duplicated(brts)])[i] == brts)
    )
  }
})

test_that("advanced check", {

  pars <- c(0.2, 0.1, 2.5, 0.12)
  n_0 <- 2
  age <- 10
  cond <- 1
  max_sims <- 5 + (is_on_ci() * 15)

  for (seed in 1:max_sims) {
    brts <- mbd::mbd_sim(
      pars = pars,
      n_0 = n_0,
      age = age,
      cond = cond,
      seed = seed
    )$brts
    out <- mbd::brts2time_intervals_and_births(brts)
    time_intervals <- out$time_intervals
    births <- out$births

    testthat::expect_equal(
      time_intervals,
      c(0, abs(diff(unique(c(brts, 0)))))
    )
    n_mbd_events <- length(unique(brts[duplicated(brts)]))
    for (i in 1:n_mbd_events) {
      testthat::expect_equal(
        births[births != 1][i + 1],
        sum(unique(brts[duplicated(brts)])[i] == brts)
      )
    }
    testthat::expect_equal(
      length(time_intervals),
      length(births)
    )
  }
})
