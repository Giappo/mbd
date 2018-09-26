context("mbd_ml2")

test_that("use", {

  skip("WIP")
  set.seed(10)
  lambda <- 0.3
  mu <- 0.1
  nu <- 0.11
  q < 0.15
  test_pars <- c(lambda, mu, nu, q)
  idparsfix <- c(1, 2, 3)
  out_classic <- mbd_ml(
    brts = c(1, 2, 3),
    initparsopt = 0.11,
    idparsopt = 4,
    idparsfix = idparsfix,
    parsfix = test_pars[idparsfix],
    missnumspec = 0,
    cond = 1,
    soc = 2,
    verbose = FALSE
  )
  expect_equal(lambda, out_classic$lambda)
  expect_equal(mu, out_classic$mu)
  expect_equal(nu, out_classic$nu)
  expect_equal(0.0010006043479098515256, out_classic$q)
  
  out <- mbd_ml2(
    branching_times = c(1, 2, 3),
    init_param_values = create_mbd_params(
      lambda = lambda,
      mu = mu,
      nu = nu,
      q = q
    ),
    fixed_params = c("lambda", "mu", "nu"),
    estimated_params = c("q"),
    init_n_species = 2,
    n_missing_species = 0,
    conditioned_on = "non_extinction"
  )
})
  
test_that("abuse", {

  skip("WIP")
  
  mbd_params <- create_mbd_params(0.1, 0.2, 0.3, 0.4)
  fixed_params <- c("lambda", "mu")
  estimated_params <- c("nu", "q")

  expect_silent(
    mbd_ml2(
      branching_times = c(1, 2, 3),
      init_param_values = mbd_params,
      fixed_params = fixed_params,
      estimated_params = estimated_params
    )
  )
  expect_error(
    mbd_ml2(
      branching_times = "nonsense",
      init_param_values = mbd_params,
      fixed_params = fixed_params,
      estimated_params = estimated_params
    ),
    "'branching_times' must be numeric"
  )

  expect_error(
    mbd_ml2(
      branching_times = c(1, 2, -34.56),
      init_param_values = mbd_params,
      fixed_params = fixed_params,
      estimated_params = estimated_params
    ),
    "All 'branching_times' must be positive"
  )
  
  expect_error(
    mbd_ml2(
      branching_times = c(1, 2, 3),
      init_param_values = "nonsense",
      fixed_params = fixed_params,
      estimated_params = estimated_params
    ),
    paste0(
      "'init_param_values' must be an mbd_params, ",
      "as created by 'create_mbd_params'"
    )
  )
  expect_error(
    mbd_ml2(
      branching_times = c(1, 2, 3),
      init_param_values = mbd_params,
      fixed_params = "nonsense",
      estimated_params = estimated_params
    ),
    "'fixed_params' must be a set of MBD parameter names"
  )
  expect_error(
    mbd_ml2(
      branching_times = c(1, 2, 3),
      init_param_values = mbd_params,
      fixed_params = c(fixed_params, fixed_params),
      estimated_params = estimated_params
    ),
    "'fixed_params' must contain unique entries only"
  )
  expect_error(
    mbd_ml2(
      branching_times = c(1, 2, 3),
      init_param_values = mbd_params,
      fixed_params = fixed_params,
      estimated_params = "nonsense"
    ),
    "'estimated_params' must be a set of MBD parameter names"
  )
  expect_error(
    mbd_ml2(
      branching_times = c(1, 2, 3),
      init_param_values = mbd_params,
      fixed_params = fixed_params,
      estimated_params = c(estimated_params, estimated_params) 
    ),
    "'estimated_params' must contain unique entries only"
  )
})
  

