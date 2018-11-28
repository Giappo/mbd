context("brts2time_intervals_and_births")

test_that("basic use", {

  brts <- c(10, 9, 8, 7, 7, 6, 6, 6, 5, 5, 5, 5)
  out <- brts2time_intervals_and_births(brts) # nolint internal function
  time_intervals <- out$time_intervals
  births <- out$births

  testthat::expect_equal(
    time_intervals,
    c(0, abs(diff(unique(c(brts, 0)))))
  )
  for (i in 1:length(unique(brts[duplicated(brts)]))) {
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
  max_sims <- 5 + (ribir::is_on_travis() * 15)

  for (s in 1:max_sims) {
    set.seed(s)
    brts <- mbd_sim(
      pars = pars,
      n_0 = n_0,
      age = age,
      cond = cond
    )$brts
    out <- brts2time_intervals_and_births(brts) # nolint internal function
    time_intervals <- out$time_intervals
    births <- out$births

    testthat::expect_equal(
      time_intervals,
      c(0, abs(diff(unique(c(brts, 0)))))
    )
    for (i in 1:length(unique(brts[duplicated(brts)]))) {
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