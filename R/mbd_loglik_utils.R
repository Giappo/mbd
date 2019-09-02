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

# integrates func from 0 to t1
# *if* this function returns, the result doesn't contains
# any negative number
#' @export
mbd_solve2 <- function(
  q_vector,
  time_interval,
  func = mbd_loglik_rhs,
  matrix_a,
  debug_mode = FALSE
) {

  atol_min <- 1e-6
  rtol_min <- 1e-10
  atol_max <- 1e-24
  rtol_max <- 1e-100
  atol_step <- -2
  rtol_step <- -10

  atol_vec <-
    10 ^ seq(from = log10(atol_min), to = log10(atol_max), by = atol_step)
  rtol_vec <-
    10 ^ seq(from = log10(rtol_min), to = log10(rtol_max), by = rtol_step)
  
  atol_len <- length(atol_vec)
  rtol_len <- length(rtol_vec)

  y <- q_vector
  parms <- matrix_a
  t1 <- time_interval
  if (debug_mode == TRUE) {
    cat("----------------------------------\n")
  }
  g <- 10         # granularity
  t0 <- 0
  # atol <- 1e-4;   # something reasonable, hopefully
  # rtol <- 1e-5;   # something reasonable, hopefully
  aa <- rr <- 1
  while (TRUE) {
    atol <- atol_vec[aa]
    rtol <- rtol_vec[rr]
    if (debug_mode == TRUE) {
      cat(atol, rtol, t0, t1, "\n")
    }
    tseq <- seq(t0, t1, length.out = g)
    x <- capture.output(my_try_catch(
      out <- deSolve::ode(
        y = y,
        times = tseq,
        func = func,
        parms = parms,
        atol = atol,
        rtol = rtol,
        tcrit = t1
      )
    ))
    # istate = attributes(out)$istate
    # rstate = attributes(out)$rstate
    lkg <- 0    # last known good
    for (ff in 1:g) {
      a <- any(out[ff, -1] < 0)
      if (!is.na(a) && a) { 
        break;
      }
      lkg <- lkg + 1
    }
    if (lkg == g) {
      break;  # done and dusted
    }
    if (lkg > 1) {
      # trace back to last known good and try from there
      t0 <- as.numeric(out[lkg, 1])
      y <- as.numeric(out[lkg, -1])
      # relax tol to default
      # atol <- 1e-4
      # rtol <- 1e-5
      aa <- rr <- 1
    } else {
      # no progress, make tol more strict
      # atol <- atol / 100
      # rtol <- rtol / 100
      if (aa >= atol_len) {
        if (debug_mode == TRUE) {
          plot(y, main = "Vector obtained with smallest tolerance")
        }
        stop("Integration failed")
      }
      if (any(grepl(x = x, pattern = "too much accuracy requested"))) {
        rr <- 1
        aa <- aa + 1
      } else {
        if (rr < rtol_len) {
          rr <- rr + 1
        } else {
          rr <- 1
          aa <- aa + 1
        }
      }
    }
  }
  out[g, -1]
}

# integrates func from 0 to t1 (old version)
# *if* this function returns, the result doesn't contains
# any negative number
#' @export
mbd_solve <- function(
  q_vector,
  time_interval,
  func = mbd_loglik_rhs,
  matrix_a,
  debug_mode = FALSE
) {

  y <- q_vector
  parms <- matrix_a
  t1 <- time_interval

  if (debug_mode == TRUE) {
    cat("----------------------------------\n")
  }
  g = 10         # granularity
  t0 = 0
  start_rtol <- 1e-8
  atol = 1e-100;  # realistically zero
  rtol = start_rtol;   # something reasonable, hopefully
  while (TRUE) {
    tseq = seq(t0, t1, length.out = g)
    if (debug_mode == TRUE) {
      cat(atol, rtol, t0, t1, "\n")
    }
    out <- deSolve::ode(
      y = y,
      times = tseq,
      func = func,
      parms = parms,
      atol = atol,
      rtol = rtol,
      tcrit = t1
    )
    istate = attributes(out)$istate
    rstate = attributes(out)$rstate
    lkg = 0    # last known good
    for (ff in 1:g) {
      a = any(out[ff, -1] < 0)
      if (!is.na(a) && a) { 
        break;
      }
      lkg = lkg + 1
    }
    if (lkg == g) {
      break;  # done and dusted
    }
    if (lkg > 1) {
      # trace back to last known good and try from there
      t0 = as.numeric(out[lkg, 1])
      y = as.numeric(out[lkg, -1])
      # relax tol to default
      rtol = start_rtol
    }
    else {
      # no progress, make tol more strict
      rtol = rtol / 100
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
    testit::assert(colSums(matrix)[-(lx + 1)] >= 0)
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
  testit::assert(sum(time_intervals) == max(abs(brts)))
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

#' @author Giovanni Laudanno
#' @title The q-vector
#' @description Yields all the values of the q-vector obtained from a likelihood
#'  computation.
#' @inheritParams default_params_doc
#' @param correct_negatives Do you want to use the negative correction? See
#'  \link{negatives_correction}.
#' @return The q-vector in time
#' @noRd
mbd_calculate_q_vector <- function(
  pars,
  brts,
  n_0 = 2,
  missnumspec = 0,
  lx = 1 + 2 * (length(brts) + length(missnumspec)),
  methode = "lsodes",
  q_threshold = 1e-3,
  abstol = 1e-16,
  reltol = 1e-10
) {
  # BASIC SETTINGS AND CHECKS
  check_brts(brts = brts, n_0 = n_0)
  if (
    check_pars(
      pars = pars,
      safety_checks = TRUE,
      q_threshold = q_threshold
    ) == "wrong"
  ) {
    return(-Inf)
  }

  # Adjusting data
  data <- brts2time_intervals_and_births(brts) # nolint internal function
  time_intervals <- data$time_intervals
  births <- data$births
  lt <- length(time_intervals)
  testit::assert(n_0 - 1 + length(brts) == n_0 + sum(births)) #every tip is born

  # LIKELIHOOD INTEGRATION

  # Setting initial conditions (there's always a +1 because of Q0)
  q_i <- c(1, rep(0, lx))
  q_t <- matrix(0, ncol = (lx + 1), nrow = lt)
  q_t[1, ] <- q_i
  dimnames(q_t)[[2]] <- paste0("Q", 0:lx)
  k <- n_0 # n_0 is the number of species at t = 1
  # t is starting from 2 so all is ok with births[t] and time_intervals[t]
  t <- 2
  D <- C <- rep(1, lt)

  # Evolving the initial state to the present
  while (t <= lt) {

    # Creating A matrix
    matrix_a <- create_a(
      pars = pars,
      lx = lx,
      k = k
    )

    # Applying A operator
    q_t[t, ] <- mbd_solve(
      q_vector = q_t[(t - 1), ],
      matrix_a = matrix_a,
      time_interval = time_intervals[t],
      debug_mode = TRUE
    )

    # Applying C operator (this is a trick to avoid precision issues)
    C[t] <- 1 / sum(sort(q_t[t, ]))
    q_t[t, ] <- q_t[t, ] * C[t]

    # Loop has to end after integrating to t_p
    if (!(t < lt)) {
      break
    }

    # Creating B matrix
    matrix_b <- create_b(
      pars = pars,
      lx = lx,
      k = k,
      b = births[t]
    )

    # Applying B operator
    q_t[t, ] <- (matrix_b %*% q_t[t, ])

    # Applying D operator (this works exactly like C)
    D[t] <- 1 / sum(sort(q_t[t, ]))
    q_t[t, ] <- q_t[t, ] * D[t]

    # Updating running parameters
    k <- k + births[t]
    t <- t + 1
  }
  q_f <- q_t[t, ]
  m <- 0:(ncol(q_t) - 1)
  one_over_cm <- (3 * (m + 1)) / (m + 3)
  one_over_qm_binom <- 1 / choose(m + n_0, n_0)
  vm <- one_over_qm_binom * one_over_cm
  q_f2 <- vm * q_f / (prod(C) * prod(D))
  q_f2
}
