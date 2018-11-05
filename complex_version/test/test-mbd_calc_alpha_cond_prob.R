context("mbd_calc_alpha_cond_prob")

test_that("abuse", {

  testthat::expect_error(
    mbd_calc_alpha_cond_prob(
      brts = c(1, 2, 3),
      pars = c(0.1, 0.2, 0.3, 0.4),
      alpha = 10^20 # Big
    ),
    "Maximum number of species exceeds R's memory limits"
  )
})
