context("package style")

lint_folder <- function(foldername, ...) {
  if (dir.exists(foldername)) {
    lintr::lint_dir(foldername, ...)
  }
}

test_that("package style", {
  testthat::expect_true(
    length(lint_folder("R")) == 0
  )
  testthat::expect_true(
    length(lint_folder("vignettes")) == 0
  )
})
