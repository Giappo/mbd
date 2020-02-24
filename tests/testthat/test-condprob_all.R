context("condprob_all")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

# test all condprobs ----
testthat::test_that("test all condprobs", {

  if (is_on_ci()) {
    absorbs <- c(TRUE, FALSE)
    brts <- c(4)
  } else {
    absorbs <- c(TRUE)
    brts <- c(3)
  }
  fortrans <- c(TRUE, FALSE)
  eqs <- c("p_eq", "q_eq")
  pars <- c(0.3, 0.15, 1.8, 0.11)

  nas <- rep(NA, length(absorbs) * length(fortrans) * length(eqs))
  df <- data.frame(
    pc = nas,
    eq = nas,
    fortran = nas,
    absorb = nas
  )
  i <- 1
  for (absorb in absorbs) {
    for (eq in eqs) {
      if (absorb) {
        lx <- 30
      } else {
        lx <- 65
      }
      for (fortran in fortrans) {
        print(i)
        df$pc[i] <- mbd::calculate_condprob(
          pars = pars,
          brts = brts,
          lx = lx,
          eq = eq,
          fortran = fortran,
          absorb = absorb
        )
        df$absorb[i] <- absorb
        df$fortran[i] <- fortran
        df$eq[i] <- eq
        i <- i + 1
      }
    }
  }

  pc_best <- df$pc[df$eq == "p_eq" & df$absorb == TRUE & df$fortran == FALSE]
  for (i in seq_len(nrow(df))) {
    testthat::expect_equal(
      df$pc[i],
      pc_best,
      tolerance = 1e-4
    )
  }

  # test nee's approx
  pc_nee <- mbd::calculate_condprob(
    pars = pars,
    brts = brts,
    lx = lx,
    eq = "nee",
    fortran = fortran,
    absorb = absorb
  )
  testthat::expect_equal(
    pc_nee,
    pc_best,
    tolerance = 1e-3
  )
})
