context("mbd_loglik_utils")

# check_sum_probs  ----
test_that("check_sum_probs ", {

  testthat::expect_error(
    utils::capture.output(mbd::check_sum_probs(
      sum_probs_1 = -1,
      sum_probs_2 = 1,
      debug_mode = FALSE
    ))
  )

  testthat::expect_error(
    utils::capture.output(mbd::check_sum_probs(
      sum_probs_1 = 1,
      sum_probs_2 = -1,
      debug_mode = FALSE
    ))
  )

  testthat::expect_error(
    utils::capture.output(mbd::check_sum_probs(
      sum_probs_1 = NA,
      sum_probs_2 = 1,
      debug_mode = FALSE
    ))
  )

  testthat::expect_error(
    utils::capture.output(mbd::check_sum_probs(
      sum_probs_1 = 1,
      sum_probs_2 = NA,
      debug_mode = FALSE
    ))
  )

})
