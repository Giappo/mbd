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

#' The A operator is given by the integration of a set of differential equations
#' between two consecutive nodes. So, defined the set in the time interval
#' [t_{i-1}, t_i], where k species are present in the phylogeny, as:
#'
#' d
#' --Q^k_m(t) = SUM_n(M^k_m,n * Q^k_n(t)
#' dt
#'
#' where m, n, label the amount of unseen species in the phylogeny,
#' A is thus defined as:
#'
#' A(t_i - t_{i-1}) = exp(M(t_k - t_{k-1})
#' @inheritParams default_params_doc
#' @param atol_min min for atol
#' @param atol_max max for atol
#' @param atol_step step for atol
#' @param rtol_min min for rtol
#' @param rtol_max max for rtol
#' @param rtol_step step for rtol
#' @author Giovanni Laudanno
#' @export
a_operator <- function(
  q_vector,
  transition_matrix,
  time_interval,
  methode = "lsoda",
  abstol = 1e-10,
  reltol = 1e-12,
  precision = 50L,
  atol_min = 1e-8,
  rtol_min = 1e-10,
  atol_max = 1e-20,
  rtol_max = 1e-88,
  atol_step = -2,
  rtol_step = -6
) {
  matrix_a <- transition_matrix
  pippobaudo <- abstol * reltol * precision; rm(pippobaudo)

  atol_vec <-
    10 ^ seq(from = log10(atol_min), to = log10(atol_max), by = atol_step)
  rtol_vec <-
    10 ^ seq(from = log10(rtol_min), to = log10(rtol_max), by = rtol_step)
  rtol_vec <- c(0, rtol_vec)

  atol_len <- length(atol_vec)
  rtol_len <- length(rtol_vec)

  good_result <- FALSE
  a <- r <- 1
  while (
    !good_result &&
    a <= atol_len &&
    r <= rtol_len
  ) {
    atol <- atol_vec[a]
    rtol <- rtol_vec[r]

    x <- capture.output(
      my_try_catch(
        result <- deSolve::ode(
          y = q_vector,
          times = c(0, time_interval),
          func = mbd_loglik_rhs,
          parms = matrix_a,
          atol = atol,
          rtol = rtol,
          method = methode
        )[2, -1]
      )
    )
    result <- correct_negatives(result)
    good_result <-
      all(!is.na(result)) &&
      all(!is.nan(result)) &&
      any(result != q_vector) &&
      all(result >= 0)
    if (!good_result) {
      if (any(grepl(x = x, pattern = "too much accuracy requested"))) {
        r <- 1
        a <- a + 1
      } else {
        if (r < rtol_len) {
          r <- r + 1
        } else {
          r <- 1
          a <- a + 1
        }
      }
    } else {
      return(result)
    }
  }

  # If errors are present, are they negligible?
  tolerance <- 1e-10
  if (!good_result) {
    negative_coords <- which(result < 0)
    for (coord in negative_coords) {
      if (abs(result[coord]) / sum(abs(result)) <= tolerance) {
        result[coord] <- 0
      }
    }
  }
  good_result <-
    all(!is.na(result)) &&
    all(!is.nan(result)) &&
    any(result != q_vector) &&
    all(result >= 0)
  if (good_result) {
    return(result)
  } else {
    stop("Integration failed.")
  }
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
mbd_solve <- function(
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
#' @noRd
mbd_solve2 <- function(y, t1, func, parms) {
  cat("----------------------------------\n")
  g = 10         # granularity
  t0 = 0
  atol = 1e-4;   # something reasonable, hopefully
  rtol = 1e-5;   # something reasonable, hopefully
  while (TRUE) {
    cat(atol, rtol, t0, t1, "\n")
    tseq = seq(t0, t1, length.out = g)
    out <- deSolve::ode(
      y = y,
      times = tseq,
      func = func,
      parms = parms,
      atol = atol,
      rtol = rtol,
      tcrit = t1
    )
    # istate = attributes(out)$istate
    # rstate = attributes(out)$rstate
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
      atol = 1e-4
      rtol = 1e-5
    }
    else {
      # no progress, make tol more strict
      atol = atol / 100
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
#' @details This is not to be called by the user.
#' @author Giovanni Laudanno
#' @export
create_a2 <- function(
  pars,
  k,
  lx,
  matrix_builder = mbd:::hyper_a_hanno,
  adjust_last_entry = TRUE
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
  diag(matrix) <- -nu * (1 - (1 - q) ^ (k + nvec)) - (lambda + mu) * (k + nvec)
  matrix[row(matrix) == col(matrix) - 1] <- mu * nvec[2:(lx + 1)]
  matrix[row(matrix) == col(matrix) + 1] <-
    matrix[row(matrix) == col(matrix) + 1] + lambda * (nvec[1:(lx)] + 2 * k)
  if (adjust_last_entry == TRUE) {
    matrix[length(nvec), length(nvec)] <- (-mu) * (k + nvec[lx + 1])
  }

  transition_matrix <- matrix
  transition_matrix
}

#' @title The A matrix
#' @description Creates the A matrix,
#' used for likelihood integration between branching times.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @author Giovanni Laudanno
#' @export
create_a <- function(
  pars,
  k,
  lx,
  matrix_builder = mbd:::hyper_a_hanno,
  adjust_last_entry = TRUE
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
  
  # diagonal
  for (n in 0:(lx - 1)) {
    limit <- min(n + k, lx - n)
    avec <- 1:limit
    matrix[n + 1, n + 1] <- -nu *
      sum(
        choose(n + k, avec) * (q ^ avec) * (1 - q) ^ (n + k - avec)
      )
  }
  diag(matrix) <- diag(matrix) - (lambda + mu) * (k + nvec)
  
  # mu terms
  matrix[row(matrix) == col(matrix) - 1] <- mu * nvec[2:(lx + 1)]
  
  # lambda terms
  matrix[row(matrix) == col(matrix) + 1] <-
    matrix[row(matrix) == col(matrix) + 1] + lambda * (nvec[1:(lx)] + 2 * k)
  
  # it is forbidden to speciate outside of the matrix
  if (adjust_last_entry == TRUE) {
    matrix[length(nvec), length(nvec)] <- (-mu) * (k + nvec[lx + 1])
  }

  # check if probabilities are imported into the system
  testit::assert(colSums(matrix)[-(lx + 1)] >= 0)

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

#' @title Test consistency of branching times
#' @description Test if the given branching times
#' can actually be generated by a MBD process.
#' In any moment you CANNOT have more births than number of species.
#' @inheritParams default_params_doc
#' @return TRUE or FALSE.
#' @author Giovanni Laudanno
#' @export
check_brts_consistency <- function(brts, n_0) {
  if (sum(brts == max(brts)) > 1) {
    stop("Crown/stem age has to be reported only once in the branching times.")
  }
  births <- brts2time_intervals_and_births(brts)$births # nolint internal function
  kvec <- n_0 + cumsum(c(0, births))
  kvec
  all(births <= kvec[-length(kvec)])
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

#' Checks for NA, NaN or negative components in a vector (usually used for q_t)
#' @inheritParams default_params_doc
#' @param v a vector
#' @param display_output If TRUE it prints the flags
#' @noRd
negatives_correction <- function(v, pars, display_output = FALSE) {
  problems <- 0
  if (any(is.na(v))) {
    problems <- 1
    na_components <- which(is.na(v) & !is.nan(v))
    nan_components <- which(is.nan(v))
    if (display_output == TRUE && !missing(pars)) {
      cat("There are non-numeric components for par values:", pars, "\n")
      if (length(na_components) > 0) {
        cat("NA component are:", na_components)
      }
      if (length(nan_components) > 0) {
        cat("NaN component are:", nan_components)
      }
    }
  }

  if (any(v < 0) && problems == 0) {
    v[v < 0 & (abs(v) / abs(max(v))) < 1e-10] <- 0
  }
  v
}

#' Puts to zero extremely small negative values
#' @inheritParams default_params_doc
#' @param v a vector
#' @param tolerance tolerance
#' @noRd
correct_negatives <- function(v, tolerance = 1e-10) {
  if (any(is.na(v)) | any(is.nan(v))) {
    return(v)
  }
  if (any(v < 0)) {
    coords <- (v < 0 & (abs(v) / sum(abs(v))) < tolerance)
    v[coords] <- abs(v[coords])
  }
  return(v)
}

#' Count the number of multiple speciation events
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
mbd_count_n_spec_events <- function(brts) {
  births <- brts2time_intervals_and_births(brts)$births # nolint internal function
  sum(births > 1)
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
    q_t[t, ] <- a_operator(
      q_vector = q_t[(t - 1), ],
      transition_matrix = matrix_a,
      time_interval = time_intervals[t],
      precision = 50L,
      methode = methode,
      abstol = abstol,
      reltol = reltol
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

#' @title Check q_t
#' @description Check q_t
#' @inheritParams default_params_doc
#' @return Nothing
#' @author Giovanni Laudanno
#' @export
check_q_vector <- function(
  q_t,
  t,
  pars,
  brts,
  debug_mode = FALSE
) {
  q_vector <- q_t[t, ]
  if (any(q_vector < 0)) {
    if (debug_mode == TRUE) {
      brts2 <- c(brts, 0)
      w <- data.frame(matrix(NA, nrow = length(q_vector), ncol = 0))
      w$values <- q_vector
      w$x <- 1:length(q_vector)
      w$cols <- ifelse(sign(q_vector) > 0, "blue", "red")
      print(w$values)
      plot(
        values ~ x,
        w,
        pch = 15,
        xlab = "m",
        ylab = "Q_m^k(t)",
        main = paste0(
          "pars = ", pars[1], ", ", pars[2], ", ", pars[3], ", ", pars[4]
        ),
        sub = paste0(
          "Problems in time interval: (", brts2[t - 1], ", ", brts2[t], ")"
        )
      )
    } else {
      stop("problems: q_t is negative!")
    }
  }
}
