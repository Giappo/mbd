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
  out <- cyclocomp::cyclocomp_package("mbd")
  too_complex <- !all(out$cyclocomp < max_complexity)
  if (too_complex) {
    print(
      "These functions are too complex:"
    )
    print(
      out$name[out$cyclocomp >= max_complexity]
    )
  }
  expect_false(too_complex)
})

test_that("package style", {
  lintr::expect_lint_free()
})
