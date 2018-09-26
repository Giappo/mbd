context("mbd_sim")

test_that("abuse", {
  expect_error(
    mbd_sim(pars = "nonsense"),
    "'pars' must have four parameters"
  )
  expect_error(
    mbd_sim(pars = c(-1.0, 1.0, 1.0, 1.0)),
    "The sympatric speciation rate 'pars\\[1\\]' must be positive"
  )
  expect_error(
    mbd_sim(pars = c(1.0, -1.0, 1.0, 1.0)),
    "The extinction rate 'pars\\[2\\]' must be positive"
  )
  expect_error(
    mbd_sim(pars = c(1.0, 1.0, -1.0, 1.0)),
    paste0("The multiple allopatric speciation trigger rate ",
      "'pars\\[3\\]' must be positive"
    )
  )
  expect_error(
    mbd_sim(pars = c(1.0, 1.0, 1.0, -1.0)),
    "The single-lineage speciation probability 'pars\\[4\\]' must be positive"
  )
})
