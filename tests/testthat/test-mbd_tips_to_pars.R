context("mbd_tips_to_pars")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

test_that("calc_diff_percentile", {

  pars <- c(0.2, 0.10, 2, 0.1)
  n_0 <- 2
  tips <- 50
  age <- 10
  percentile <- 0.9
  out <- calc_diff_percentile(
    pars = pars,
    n_0 = n_0,
    tips = tips,
    age = age,
    percentile = percentile
  )
  testthat::expect_true(
    out >= 0 && out <= 1
  )
})

test_that("mbd_tips_to_pars_single", {

  if (!is_on_ci()) {
    skip("This must be performed only on ci")
  }

  tips <- 60
  age <- 7
  optim_ids <- c(FALSE, FALSE, TRUE, TRUE)
  out <- mbd_tips_to_pars_single(
    tips = tips,
    age = age,
    start_pars = c(0.2, 0.15, 1, 0.15),
    optim_ids = optim_ids,
    verbose = FALSE
  )
  testthat::expect_true(
    length(out) == sum(optim_ids)
  )
  testthat::expect_true(is.numeric(out))
  testthat::expect_true(
    all(names(out) == get_param_names()[optim_ids])
  )
})

test_that("mbd_tips_to_pars", {

  if (!is_on_ci()) {
    skip("This must be performed only on ci")
  }

  tips <- 60
  age <- 7
  optim_ids <- c(FALSE, FALSE, TRUE, TRUE)
  chain_length <- 2
  out <- mbd_tips_to_pars(
    tips = tips,
    age = age,
    start_pars = c(0.2, 0.15, 1, 0.15),
    optim_ids = optim_ids,
    chain_length = chain_length,
    verbose = FALSE
  )
  testthat::expect_true(
    length(out) == chain_length
  )
  testthat::expect_true(
    all(unlist(
      lapply(out, FUN = is.numeric)
    ))
  )
  testthat::expect_equal(
    unlist(
      lapply(out, FUN = names)
    ),
    rep(get_param_names()[optim_ids], chain_length)
  )
})

test_that("find_nu_q_relation", {

  if (!is_on_ci()) {
    skip("This must be performed only on ci")
  }

  tips <- 60
  age <- 6
  chain_length <- 2
  out <- find_nu_q_relation(
    tips = tips,
    age = age,
    chain_length = chain_length,
    verbose = FALSE
  )
  testthat::expect_true(is.numeric(out))
  testthat::expect_true(out > 0)
})
