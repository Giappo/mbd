context("mbd_ml2")

test_that("compare style", {

  branching_times <- c(1, 2, 3)
  lambda <- 0.3
  mu <- 0.1
  nu <- 0.11
  q <- 0.15
  #############################
  # Classic interface
  #############################
  set.seed(10)
  out_classic <- mbd_ml(
    brts = branching_times,
    initparsopt = q,
    idparsopt = 4,
    idparsfix = c(1, 2, 3),
    parsfix = c(lambda, mu, nu),
    missnumspec = 0,
    cond = 1,
    soc = 2,
    verbose = FALSE
  )
  #############################
  # More intuitive interface
  #############################
  set.seed(10)
  out_new <- mbd_ml2(
    branching_times = branching_times,
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
  expect_equal(out_classic, out_new)
})
  
test_that("abuse", {

  mbd_params <- create_mbd_params(0.1, 0.2, 0.3, 0.4)
  fixed_params <- c("lambda", "mu")
  estimated_params <- c("nu", "q")

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
  expect_error(
    mbd_ml2(
      branching_times = c(1, 2, 3),
      init_param_values = mbd_params,
      fixed_params = c("lambda", "mu", "q"), # q twice
      estimated_params = c("nu", "q")
    ),
    paste0(
      "'fixed_params' and 'estimated_params' together must contain each ",
      "of the four MBD parameter names"
    )
  )
  expect_error(
    mbd_ml2(
      branching_times = c(1, 2, 3),
      init_param_values = mbd_params,
      fixed_params = fixed_params,
      estimated_params = estimated_params,
      init_n_species = 0
    ),
    "'init_n_species' must be 1 or 2"
  )
  expect_error(
    mbd_ml2(
      branching_times = c(1, 2, 3),
      init_param_values = mbd_params,
      fixed_params = fixed_params,
      estimated_params = estimated_params,
      n_missing_species = -12345
    ),
    "'n_missing_species' must be positive"
  )
  expect_error(
    mbd_ml2(
      branching_times = c(1, 2, 3),
      init_param_values = mbd_params,
      fixed_params = fixed_params,
      estimated_params = estimated_params,
      conditioned_on = "nonsense"
    ),
    "'conditioned_on' must be either 'nothing' or 'non_extinction'"
  )

})
