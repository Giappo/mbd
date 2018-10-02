context("alpha_conditional_probability")

test_that("abuse", {
  
  expect_error(
    alpha_conditional_probability(
      brts = c(1, 2, 3), 
      pars = c(0.1, 0.2, 0.3, 0.4), 
      alpha = 10^20 # Big
    ),
    "Maximum number of species exceeds R's memory limits"
  )
})
