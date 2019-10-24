#' @noRd
my_try_catch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(
      expr, error = function(e) {
      err <<- e
      NULL
    }
    ), warning = function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(
    value = value,
    warning = warn,
    error = err
  )
}

#' Function to build a matrix, used in creating the A and B operators.
#' It produces the structure
#'  q ^ (m - n) * (1 - q) ^ (k + 2 * n-m) *
#'  sum_j 2 ^ j choose(k, j) * choose(n, m - n - j)
#' @inheritParams default_params_doc
#' @examples
#'   m <- hyper_a_hanno(n_species = 2, k = 2, q = 0.1)
#'   testthat::expect_equal(m[1, 1], 0.81)
#'   testthat::expect_equal(m[1, 2], 0.00)
#'   testthat::expect_equal(m[1, 3], 0.00)
#'   testthat::expect_equal(m[2, 1], 0.36)
#'   testthat::expect_equal(m[2, 2], 0.729)
#'   testthat::expect_equal(m[2, 3], 0.00)
#'   testthat::expect_equal(m[3, 1], 0.04)
#'   testthat::expect_equal(m[3, 2], 0.405)
#'   testthat::expect_equal(m[3, 3], 0.6561)
#' @noRd
#' @author Hanno Hildenbrandt, adapted by Richel J.C. Bilderbeek
hyper_a_hanno <- function(
  n_species,
  k,
  q
) {
  if (n_species > 46340) {
    stop("'n_species' must be below 46340. ",
         "Cannot allocate matrix with 2^31 elements")
  }
  # HG function: fast O(N), updated after Moulis meeting
  j <- 0:k
  a_1 <- (1 - q) ^ (k) * choose(k, j) * (2) ^ j
  n_species <- n_species + 1
  matrix_a <- diag(a_1[1], nrow = n_species + 2, ncol = n_species + 2)
  matrix_a[1:(k + 1), 1] <- a_1
  for (dst in 2:n_species) {
    src <- dst - 1
    s <- src:min(n_species, 2 * src + k - 1)
    matrix_a[s + 2, dst] <- matrix_a[s, src] + matrix_a[s + 1, src]
    m <- s - 1
    n <- src - 1
    matrix_a[s, src] <- matrix_a[s, src] * q ^ (m - n) * (1 - q) ^ (2 * n - m)
  }
  matrix_a[n_species, n_species] <- matrix_a[n_species, n_species] *
    (1 - q) ^ (n_species - 1)
  matrix_a[1:n_species, 1:n_species]
}

#' @title mbd ODE system integrator
#' @description Integrates "func" in the time interval
# *if* this function returns, the result doesn't contains
# any negative number
#' @inheritParams default_params_doc
#' @param func function for the right hand side of the ODE
#' @export
#' @author Hanno Hildenbrandt, adapted by Giovanni Laudanno
mbd_solve <- function(
  vector,
  time_interval,
  func = mbd_loglik_rhs,
  parms
) {

  y <- vector
  t1 <- time_interval

  g <- 10 # granularity
  t0 <- 0
  start_rtol <- 1e-8
  atol <- 1e-100 # realistically zero
  rtol <- start_rtol # something reasonable, hopefully
  while (TRUE) {
    tseq <- seq(t0, t1, length.out = g)
    out <- deSolve::ode(
      y = y,
      times = tseq,
      func = func,
      parms = parms,
      atol = atol,
      rtol = rtol,
      tcrit = t1
    )
    # it might be useful for debug istate = attributes(out)$istate
    # it might be useful for debug rstate = attributes(out)$rstate
    lkg <- 0 # last known good
    for (ff in 1:g) {
      a <- any(out[ff, -1] < 0)
      if (!is.na(a) && a) {
        break;
      }
      lkg <- lkg + 1
    }
    if (lkg == g) {
      break; # done and dusted
    }
    if (lkg > 1) {
      # trace back to last known good and try from there
      t0 <- as.numeric(out[lkg, 1])
      y <- as.numeric(out[lkg, -1])
      # relax tol to default
      rtol <- start_rtol
    }
    else {
      # no progress, make tol more strict
      rtol <- rtol / 100
    }
  }
  out[g, -1]
}

#' @title N matrix
#' @description Creates the N matrix,
#' used to express the probabilities for multiple speciations.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @author Giovanni Laudanno
#' @export
create_n <- function(
  pars,
  k,
  b,
  lx,
  matrix_builder = hyper_a_hanno
) {
  if (b > k) {
    stop("you can't have more births than species present in the phylogeny") # nolint
  }

  q <- pars[4]

  matrix <- choose(k, b) * (q ^ b) *
    matrix_builder(n_species = lx, k = k - b, q = q)
  matrix
}

#' @title The A matrix
#' @description Creates the A matrix,
#' used for likelihood integration between branching times.
#' @inheritParams default_params_doc
#' @param no_species_out_of_the_matrix If true prevents interactions from states
#'  describing number of species greater than lx
#' @details This is not to be called by the user.
#' @author Giovanni Laudanno
#' @export
create_a <- function(
  pars,
  k,
  lx,
  matrix_builder = mbd:::hyper_a_hanno,
  no_species_out_of_the_matrix = FALSE
) {
  if (k > lx) {
    stop("The matrix is too small. Increase lx.")
  }
  lambda <- pars[1]
  mu <- pars[2]
  nu <- pars[3]
  q <- pars[4]

  testit::assert(lx < 2 ^ 31)
  nvec <- 0:lx

  matrix <- nu * create_n(
    pars = pars,
    k = k,
    b = 0,
    lx = lx,
    matrix_builder = matrix_builder
  )

  # mu terms
  matrix[row(matrix) == col(matrix) - 1] <- mu * nvec[2:(lx + 1)]

  # lambda terms
  matrix[row(matrix) == col(matrix) + 1] <-
    matrix[row(matrix) == col(matrix) + 1] + lambda * (nvec[1:(lx)] + 2 * k)

  # diagonal
  # (it is forbidden to speciate outside of the matrix)
  if (no_species_out_of_the_matrix == TRUE) {
    for (n in 0:(lx - 1)) {
      limit <- min(n + k, lx - n)
      avec <- 1:limit
      matrix[n + 1, n + 1] <- -nu *
        sum(
          choose(n + k, avec) * (q ^ avec) * (1 - q) ^ (n + k - avec)
        )
    }
    diag(matrix) <- diag(matrix) - (lambda + mu) * (k + nvec)
    matrix[length(nvec), length(nvec)] <- (-mu) * (k + nvec[lx + 1])
    # check if probabilities are imported into the system
    testit::assert(
      colSums(matrix)[-c(lx + 1)] >= -1e-13
    )
  } else {
    diag(matrix) <- -nu * (1 - (1 - q) ^ (k + nvec)) -
      (lambda + mu) * (k + nvec)
  }

  matrix_a <- matrix
  matrix_a
}

#' @title B matrix
#' @description Creates the B matrix,
#' used to modify the likelihood on branching times.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @author Giovanni Laudanno
#' @export
create_b <- function(
  pars,
  k,
  b,
  lx,
  matrix_builder = hyper_a_hanno
) {
  if (b > k) {
    stop("you can't have more births than species present in the phylogeny") # nolint
  }

  lambda <- pars[1]
  nu <- pars[3]

  n_matrix <- create_n(
    pars = pars,
    k = k,
    b = b,
    lx = lx,
    matrix_builder = matrix_builder
  )
  matrix <- lambda * k * diag(lx + 1) * (b == 1) + nu * n_matrix
  matrix
}

#' @title Builds the right hand side of the ODE set for multiple birth model
#' @description Builds the right hand side of the ODE set
#' for multiple birth model
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @author Hanno Hildebrandt
#' @export
mbd_loglik_rhs <- function(t, x, params) {
  list(params %*% x)
}

#' Converts branching times to 'time intervals between branching times'
#'   and 'birth at nodes' vectors
#' @inheritParams default_params_doc
#' @noRd
brts2time_intervals_and_births <- function(brts, brts_precision = 8) {

  brts <- DDD::roundn(brts, digits = brts_precision)
  branching_times <- unlist(unname(sort(abs(brts), decreasing = TRUE)))
  unique_branching_times <- unique(branching_times)
  time_intervals <- c(0, -diff(c(unique_branching_times, 0)))
  multiple_births_coords <- duplicated(branching_times)
  multiple_branching_times <- unique(branching_times[multiple_births_coords])
  births <- c(
    0,
    rev(c(
      unname(table(branching_times))
    ))[-1],
    0
  )

  testit::assert(length(births) == length(time_intervals))
  testit::assert(
    abs(sum(time_intervals) - max(abs(brts))) <=
      .Machine$double.eps * max(abs(brts))
  )
  testit::assert(
    all.equal(
      rev(cumsum(rev(time_intervals)))[-1], unique_branching_times
    )
  )
  testit::assert(
    all(time_intervals[-1] > .Machine$double.neg.eps) # nolint
  )
  testit::assert(
    all(time_intervals[-1] > 10 ^ (-brts_precision))
  )
  testit::assert(
    all(unique_branching_times[births > 1] == multiple_branching_times)
  )
  testit::assert(all(multiple_branching_times %in% branching_times))
  for (i in seq_along(multiple_branching_times)) {
    testit::assert(
      sum(multiple_branching_times[i] == branching_times) > 1
    )
  }

  list(
    time_intervals = time_intervals,
    births = births
  )
}

#' @noRd
max_lx <- function() {
  maximum_lx <- 1400
  maximum_lx
}

#' @noRd
check_sum_probs <- function(
  sum_probs_1,
  sum_probs_2,
  debug_mode
) {
  # Removing sum_probs_1 and sum_probs_2 effects from the LL
  if (!(all(sum_probs_1 > 0))) {
    cat("The value of sum_probs_1 is: ", sum_probs_1, "\n")
    if (debug_mode == FALSE) {
      stop("problems: sum_probs_1 is non positive!")
    }
  }
  if (!(all(sum_probs_2 > 0))) {
    cat("The value of sum_probs_2 is: ", sum_probs_2, "\n")
    if (debug_mode == FALSE) {
      stop("problems: sum_probs_2 is non positive!")
    }
  }
  if (any(is.na(sum_probs_1) | is.nan(sum_probs_1))) {
    cat("The value of sum_probs_1 is: ", sum_probs_1, "\n")
    if (debug_mode == FALSE) {
      stop("problems: sum_probs_1 is Na or NaN!")
    }
  }
  if (any(is.na(sum_probs_2) | is.nan(sum_probs_2))) {
    cat("The value of sum_probs_2 is: ", sum_probs_2, "\n")
    if (debug_mode == FALSE) {
      stop("problems: sum_probs_2 is Na or NaN!")
    }
  }
}

#' @noRd
deliver_loglik <- function(
  likelihood,
  sum_probs_1,
  sum_probs_2,
  cond,
  pc,
  debug_mode
) {
  if (debug_mode == TRUE) {
    cat("The value of the likelihood is: ", likelihood, "\n")
  }

  # Removing sum_probs_1 and sum_probs_2 effects from the LL
  check_sum_probs(
    sum_probs_1 = sum_probs_1,
    sum_probs_2 = sum_probs_2,
    debug_mode = debug_mode
  )
  loglik <- log(likelihood) + sum(log(sum_probs_1)) + sum(log(sum_probs_2))

  # Various checks
  loglik <- as.numeric(loglik)
  if (is.nan(loglik) | is.na(loglik)) {
    loglik <- -Inf
  } else {
    # conditioned likelihood
    loglik <- loglik - log(pc) * (cond > 0)
  }
  loglik
}
