context("mbd_sim_checked")

test_that("same as classic interface", {

  lambda <- 0.2 # sympatric speciation rate
  mu <- 0.15 # extinction rate
  nu <- 2.0 # multiple allopatric speciation trigger rate
  q <- 0.1 # single-lineage speciation probability
  sim_pars <- c(lambda, mu, nu, q)
  crown_age <- 1

  set.seed(42)
  classic_sim <- mbd_sim(
    pars = sim_pars,
    soc = 2,
    age = crown_age,
    cond = 1
  )

  set.seed(42)
  new_sim <- mbd_sim_checked(
    mbd_params = create_mbd_params(lambda = lambda, mu = mu, nu = nu, q = q),
    crown_age = crown_age,
    conditioned_on = "non_extinction"
  )
  expect_equal(names(classic_sim), names(new_sim))
  expect_equal(classic_sim$btrs, new_sim$btrs)
})

test_that("abuse", {

  mbd_params <- create_mbd_params(0.1, 0.2, 0.3, 0.4)
  crown_age <- 15

  expect_error(
    mbd_sim_checked(
      mbd_params = "nonsense",
      crown_age = crown_age
    ),
    "'mbd_params' must be a valid MBD parameter set"
  )
  expect_error(
    mbd_sim_checked(
      mbd_params = mbd_params,
      crown_age = -12.34
    ),
    "'crown_age' must be positive"
  )
  expect_error(
    mbd_sim_checked(
      mbd_params = mbd_params,
      stem_age = -12.34
    ),
    "'stem_age' must be positive"
  )
  expect_error(
    mbd_sim_checked(
      mbd_params = mbd_params,
      crown_age = NA,
      stem_age = NA
    ),
    "'crown_age' or 'stem_age' must be set"
  )
  expect_error(
    mbd_sim_checked(
      mbd_params = mbd_params,
      crown_age = 1.2,
      stem_age = 2.3
    ),
    "'crown_age' or 'stem_age' must be set exclusively"
  )
  expect_error(
    mbd_sim_checked(
      mbd_params = mbd_params,
      crown_age = 1.2,
      conditioned_on = "nonsense"
    ),
    "'conditioned_on' must be either 'nothing' or 'non_extinction'"
  )

})
