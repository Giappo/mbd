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
#' @author Giovanni Laudanno
#' @noRd
a_operator <- function(
  q_vector,
  transition_matrix,
  time_interval,
  precision = 50L,
  abstol = 1e-16,
  reltol = 1e-12,
  methode = "lsodes"
) {
  times <- c(0, time_interval)
  ode_matrix <- transition_matrix
  out <- list(
    value = q_vector,
    warning = 1,
    error = 1
  )
  methodes <- mbd_methodes()
  other_methodes <- methodes[!methodes == methode]
  methodes <- c(methode, other_methodes)
  max_rtol <- 1e-32
  for (methode in methodes) {
    rtol <- reltol
    while (
      (
        !is.null(out$warning) ||
        !is.null(out$error) ||
        any(out$value < 0)
      ) &&
      rtol >= max_rtol
    ) {
      x <- utils::capture.output(R.utils::withTimeout(
        out <- my_try_catch(deSolve::ode( # nolint internal function
          y = q_vector,
          times = times,
          func = mbd_loglik_rhs,
          parms = ode_matrix,
          atol = abstol,
          rtol = rtol,
          method = methode
        )[2, -1]),
        timeout = 1001
      ))
      rtol <- rtol * 1e-4
    }
    if (is.null(out$warning) && is.null(out$error)) {
      break
    } else {
      out$value <- rep(NA, length(q_vector))
    }
  }

  if (!is.null(out$warning) || !is.null(out$error)) {
    result <- rep(NA, length(q_vector))
  } else {
    result <- out$value
  }
  rm(x)
  result
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
#' @author Hanno Hildenbrand, adapted by Richel J.C. Bilderbeek
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
create_a <- function(
  pars,
  k,
  lx,
  matrix_builder = hyper_a_hanno
) {
  lambda <- pars[1]
  mu <- pars[2]
  nu <- pars[3]
  q <- pars[4]

  testit::assert(lx < 2 ^ 31)
  nvec <- 0:lx

  # matrix <- nu * matrix_builder(n_species = lx, k = k, q = q)
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
  # matrix[length(nvec), length(nvec)] <- (-mu) * (k + nvec[length(nvec)])
  matrix
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
  # q <- pars[4]

  # k2 <- k - b
  # matrix <- matrix_builder(n_species = lx, k = k2, q = q)
  # lambda * k * diag(lx + 1) * (b == 1) + nu * choose(k, b) * (q ^ b) * matrix
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
#' @author Giovanni Laudanno
#' @export
mbd_loglik_rhs <- function(t, x, params) {
  with(as.list(x), {
    starting_vector <- x
    transition_matrix <- params
    dx <- rep(0, length(starting_vector))
    dx <- drop(transition_matrix %*% starting_vector)
    out <- (dx)
    names(out) <- names(x)
    return(list(out))
  })
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
  if (any(q_vector < 0)) { # debug
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
        col = cols,
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
