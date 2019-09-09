context("package style")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

test_that("cyclomatic complexity", {

  if (is_on_ci()) {
    skip("Already performed by lintrbot.")
  }

  max_complexity <- 15
  fun_list <- ls(paste0("package:", get_pkg_name())) # nolint internal function
  for (i in seq_along(fun_list)) {
    expect_true(
      cyclocomp::cyclocomp(get(fun_list[[i]])) < max_complexity
    )
  }
})

test_that("package style", {
  lintr::expect_lint_free()
})
