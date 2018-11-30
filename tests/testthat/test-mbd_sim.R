context("mbd_sim")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

# abuse ----
test_that("abuse", {
  expect_error(
    mbd_sim(pars = "nonsense"),
    "'pars' must have four parameters"
  )
  expect_error(
    mbd_sim(pars = c(-12.34, 1.0, 1.0, 1.0)),
    "The sympatric speciation rate 'pars\\[1\\]' must be positive"
  )
  expect_error(
    mbd_sim(pars = c(1.0, -12.34, 1.0, 1.0)),
    "The extinction rate 'pars\\[2\\]' must be positive"
  )
  expect_error(
    mbd_sim(pars = c(1.0, 1.0, -12.34, 1.0)),
    paste0(
      "The multiple allopatric speciation trigger rate ",
      "'pars\\[3\\]' must be positive"
    )
  )
  expect_error(
    mbd_sim(pars = c(1.0, 1.0, 1.0, -12.34)),
    "The single-lineage speciation probability 'pars\\[4\\]' must be positive"
  )
  expect_error(
    mbd_sim(
      pars = c(1.0, 1.0, 1.0, 1.0),
      tips_interval = c(2, 1)
    ),
    "'tips_interval' must contain two values, of which the second is larger"
  )
  expect_error(
    mbd_sim(
      pars = c(1.0, 1.0, 1.0, 1.0),
      tips_interval = c(-12, -11)
    ),
    "'tips_interval' must contain two positive values"
  )
  expect_error(
    mbd_sim(
      pars = c(1.0, 1.0, 1.0, 1.0),
      n_0 = "nonsense"
    ),
    "'n_0' must be numeric"
  )
  expect_error(
    mbd_sim(
      pars = c(1.0, 1.0, 1.0, 1.0),
      n_0 = 42
    ),
    "'n_0' must be '1' of '2'"
  )
})

# mbd_sim ----
test_that("mbd_sim", {

  pars <- c(0.2, 0.15, 2, 0.1)
  n_0 <- 2
  age <- 5
  max_sims <- 15 + (is_on_ci() * 35)

  # test with cond == 0
  cond <- 0
  for (s in 1:max_sims) {
    set.seed(s)
    out <- mbd_sim(
      pars = pars,
      n_0 = n_0,
      age = age,
      cond = cond
    )
    expect_true(
      is.numeric(out$brts) && all(out$brts >= 0) &&
        all(out$brts == sort(out$brts, decreasing = TRUE))
    )
    expect_true(
      is.numeric(out$l_matrix),
      all(out$l_matrix[-1, 2] %in% out$l_matrix[, 3])
    )
    expect_true(
      all(
        (floor(out$brts * 1e2) * 1e-2) %in%
          (floor(out$l_matrix[, 1] * 1e2) * 1e-2)
      )
    )
    expect_true(
      ncol(out$l_matrix) == 4
    )
    expect_true(
      class(out$reconstructed_tree) == "phylo"
    )
    expect_true(
      class(out$full_tree) == "phylo"
    )
    expect_true(
      sum(out$full_tree$edge.length) >=
        sum(out$reconstructed_tree$edge.length)
    )
    expect_true(
      length(out$reconstructed_tree$tip.label) == sum(out$l_matrix[, 4] == -1)
    )
  }

  # test with cond == 1
  cond <- 1
  for (s in (max_sims + 1):(2 * max_sims)) {
    set.seed(s)
    out <- mbd_sim(
      pars = pars,
      n_0 = n_0,
      age = age,
      cond = cond
    )
    expect_true(
      is.numeric(out$brts) && all(out$brts >= 0) &&
        all(out$brts == sort(out$brts, decreasing = TRUE))
    )
    expect_true(
      is.numeric(out$l_matrix),
      all(out$l_matrix[-1, 2] %in% out$l_matrix[, 3])
    )
    expect_true(
      all(
        (floor(out$brts * 1e2) * 1e-2) %in%
          (floor(out$l_matrix[, 1] * 1e2) * 1e-2)
      )
    )
    expect_true(
      ncol(out$l_matrix) == 4
    )
    expect_true(
      class(out$reconstructed_tree) == "phylo"
    )
    expect_true(
      class(out$full_tree) == "phylo"
    )
    expect_true(
      sum(out$full_tree$edge.length) >=
        sum(out$reconstructed_tree$edge.length)
    )
    expect_true(
      length(out$reconstructed_tree$tip.label) == sum(out$l_matrix[, 4] == -1)
    )
    expect_true(
      # check conditioning
      any(out$l_matrix[, 4] == -1)
    )
  }
})
