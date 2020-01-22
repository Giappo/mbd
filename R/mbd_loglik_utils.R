# Tools to calculate the likelihood -----

#' Function to build a matrix, used in creating the A and B operators.
#' It produces the structure
#'  q ^ (m - n) * (1 - q) ^ (k + 2 * n-m) *
#'  sum_j 2 ^ j choose(k, j) * choose(n, m - n - j)
#' @inheritParams default_params_doc
#' @examples
#'   m <- mbd::hyper_matrix0(n_species = 2, k = 2, q = 0.1)
#'   testthat::expect_equal(m[1, 1], 0.81)
#'   testthat::expect_equal(m[1, 2], 0.00)
#'   testthat::expect_equal(m[1, 3], 0.00)
#'   testthat::expect_equal(m[2, 1], 0.36)
#'   testthat::expect_equal(m[2, 2], 0.729)
#'   testthat::expect_equal(m[2, 3], 0.00)
#'   testthat::expect_equal(m[3, 1], 0.04)
#'   testthat::expect_equal(m[3, 2], 0.405)
#'   testthat::expect_equal(m[3, 3], 0.6561)
#' @author Hanno Hildenbrandt, adapted by Richel J.C. Bilderbeek
#' @export
hyper_matrix0 <- function(
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
    matrix_a[s + 2, dst] <- matrix_a[s, src] + matrix_a[s + 1, src] # big number
    m <- s - 1
    n <- src - 1
    matrix_a[s, src] <- matrix_a[s, src] * q ^ (m - n) * (1 - q) ^ (2 * n - m)
  }
  matrix_a[n_species, n_species] <- matrix_a[n_species, n_species] *
    (1 - q) ^ (n_species - 1)
  matrix_a[1:n_species, 1:n_species]
}

#' Function to build a matrix, used in creating the A and B operators.
#' It produces the structure
#'  q ^ (m - n) * (1 - q) ^ (k + 2 * n-m) *
#'  sum_j 2 ^ j choose(k, j) * choose(n, m - n - j)
#' @inheritParams default_params_doc
#' @author Hanno Hildenbrandt
#' @export
hyper_matrix1 <- function(
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
  mat <- diag(a_1[1], nrow = n_species + 2, ncol = n_species)
  mat[1:(k + 1), 1] <- a_1
  for (dst in 2:n_species) {
    src <- dst - 1
    s <- src:min(n_species, 2 * src + k - 1)
    m <- s - 1
    n <- src - 1
    mat[s + 2, dst] <- mat[s, src] + mat[s + 1, src]
    log_col <- log2(mat[s, src])
    a <- log_col / log2(q)
    mat[s, src] <- (q ^ (m - n + a) * (1 - q) ^ (2 * n - m))
  }
  mat[n_species, n_species] <-
    mat[n_species, n_species] * (1 - q) ^ (n_species - 1)
  mat[1:n_species, 1:n_species]
}

#' Function to build a matrix, used in creating the A and B operators.
#' It produces the structure
#'  q ^ (m - n) * (1 - q) ^ (k + 2 * n-m) *
#'  sum_j 2 ^ j choose(k, j) * choose(n, m - n - j)
#' @inheritParams default_params_doc
#' @author Hanno Hildenbrandt, adapted by Giovanni Laudanno
#' @export
hyper_matrix2 <- function(
  n_species,
  k,
  q
) {
  if (n_species > 46340) {
    stop("'n_species' must be below 46340. ",
         "Cannot allocate matrix with 2^31 elements")
  }
  # HG function: fast O(N), updated after Moulis meeting
  func <- function(m, n) (q ^ (m - n) * (1 - q) ^ (2 * n - m))
  m_min <- k
  n_min <- ceiling(k / 2)
  mins <- func(m_min, n_min)

  j <- 0:k
  a_1 <- (1 - q) ^ (k) * choose(k, j) * (2) ^ j
  a_1 <- a_1 * mins
  n_species <- n_species + 1
  mat <- diag(a_1[1], nrow = n_species + 2, ncol = n_species)
  mat[1:(k + 1), 1] <- a_1
  for (dst in 2:n_species) {
    src <- dst - 1
    s <- src:min(n_species, 2 * src + k - 1)
    m <- s - 1
    n <- src - 1
    mat[s + 2, dst] <- mat[s, src] + mat[s + 1, src]
    log_col <- log2(mat[s, src])
    a <- log_col / log2(q)
    mat[s, src] <- (q ^ (a + (m - m_min) - (n - n_min)) *
                    (1 - q) ^ (2 * (n - n_min) - (m - m_min)))
  }
  mat[n_species, n_species] <-
    mat[n_species, n_species] * (1 - q) ^ (n_species - 1)
  mat[!is.finite(mat)] <- 0

  mat[1:n_species, 1:n_species]
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
  matrix_builder = mbd::hyper_matrix2
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
  matrix_builder = mbd::hyper_matrix2,
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

  matrix <- nu * mbd::create_n(
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

#' @title The A matrix
#' @description Creates the A matrix but in a slower, more understandable way.
#' @inheritParams default_params_doc
#' @param no_species_out_of_the_matrix If true prevents interactions from states
#'  describing number of species greater than lx
#' @details This is not to be called by the user.
#' @author Giovanni Laudanno
#' @export
create_a_slow <- function(
  pars,
  k,
  lx,
  no_species_out_of_the_matrix = FALSE
) {

  lambda <- pars[1]
  mu <- pars[2]
  nu <- pars[3]
  q <- pars[4]
  mvec <- 1:lx - 1

  mat <- matrix(0, lx + 1, lx + 1)
  mat[col(mat) == row(mat) + 1] <- mu * (mvec + 1)

  # lower triangular matrix: m > n
  griglia <- expand.grid(m = 1:lx, n = 0:(lx - 1))
  griglia <- unname(as.matrix(griglia[griglia[, 1] > griglia[, 2], ]))
  for (i in seq_len(nrow(griglia))) {
    m <- griglia[i, 1]
    n <- griglia[i, 2]
    j <- 0:min(m - n, k)
    j <- j[n >= (m - n - j)]
    if (length(j) > 0) {
      mn1 <- log(nu) +
        log(1 - q) * k +
        log(q) * (m - n) +
        log(1 - q) * (2 * n - m)
      mn2s <- log(2) * j + lchoose(k, j) + lchoose(n, m - n - j)

      min_mn2 <- min(mn2s)
      a_matrix_mn <- sum(exp(mn2s - min_mn2)) * exp(min_mn2 + mn1)
      if (
        is.nan(a_matrix_mn) ||
        a_matrix_mn == 0 ||
        is.infinite(a_matrix_mn)
      ) {
        max_mn2 <- max(mn2s)
        a_matrix_mn <- sum(exp(mn2s - max_mn2)) * exp(max_mn2 + mn1)
      }
    } else {
      a_matrix_mn <- 0
    }
    mat[m + 1, n + 1] <- a_matrix_mn
  }

  mat[col(mat) == row(mat) - 1] <- mat[col(mat) == row(mat) - 1] +
    lambda * (mvec + 2 * k)

  # main diagonal
  m <- 0:lx
  nu_terms <- rep(0, lx + 1)
  if (no_species_out_of_the_matrix == TRUE) {
    for (n in 0:(lx - 1)) {
      limit <- min(n + k, lx - n)
      avec <- 1:limit
      nu_terms[n + 1] <-
        sum(
          choose(n + k, avec) * (q ^ avec) * (1 - q) ^ (n + k - avec)
        )
    }
  } else {
    nu_terms <- (1 - (1 - q) ^ (m + k))
  }
  diag(mat) <-
    -nu * nu_terms +
    -mu * (m + k) +
    -lambda * (m + k) +
    no_species_out_of_the_matrix * lambda * c(rep(0, lx), 1) * (m + k)
  mat
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
  matrix_builder = mbd::hyper_matrix2
) {
  if (b > k) {
    stop("you can't have more births than species present in the phylogeny") # nolint
  }

  lambda <- pars[1]
  nu <- pars[3]

  n_matrix <- mbd::create_n(
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

# Utility tools -----

#' Maximum allowed value for lx
#' @export
max_lx <- function() {
  maximum_lx <- 1400
  maximum_lx
}

#' Default value for lx
#' @export
default_lx <- function(
  brts,
  missnumspec = 0
) {
  def_lx <- min(1 + 2 * (length(brts) + max(missnumspec)), mbd::max_lx())
  def_lx
}

#' Test for pure birth in \link{mbd_loglik}
#' @inheritParams default_params_doc
#' @export
is_pbd <- function(
  pars,
  tips_interval,
  cond,
  missnumspec,
  n_0
) {
  is_it_pure_birth <-
    pars[2] == 0 &&
    all(tips_interval == c(n_0 * (cond > 0), Inf)) &&
    missnumspec == 0
  is_it_pure_birth
}

#' Approximate the branching times to a fixed resolution. Events separated by
#'  a time smaller than 10^-brts_precision can be considered simultaneous
#' @inheritParams default_params_doc
#' @export
approximate_brts <- function(
  brts,
  brts_precision = 8
) {
  brts <- DDD::roundn(brts, digits = brts_precision)
  brts[brts <= 10 ^ (-brts_precision)] <- 10 ^ (-brts_precision)
  brts
}

#' Converts branching times to 'time intervals between branching times'
#'   and 'birth at nodes' vectors
#' @inheritParams default_params_doc
#' @export
brts2time_intervals_and_births <- function(
  brts,
  brts_precision = 8
) {

  brts <- mbd::approximate_brts(brts = brts, brts_precision = brts_precision)

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
    all(time_intervals[-1] >= .Machine$double.neg.eps) # nolint
  )
  testit::assert(
    all(time_intervals[-1] >= 10 ^ (-brts_precision))
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

#' Initialize the vector q_t for \link{mbd_loglik}
#' @inheritParams default_params_doc
#' @param lt length of time intervals vector
#' @export
initialize_q_t <- function(
  lx,
  lt
) {
  q_i <- c(1, rep(0, lx))
  q_t <- matrix(0, ncol = (lx + 1), nrow = lt)
  q_t[1, ] <- q_i
  dimnames(q_t)[[2]] <- paste0("Q", 0:lx)
  q_t
}

#' Check the vectors of probability sums computed during likelihood integration
#' @inheritParams default_params_doc
#' @export
check_sum_probs <- function(
  sum_probs_1,
  sum_probs_2
) {
  # Removing sum_probs_1 and sum_probs_2 effects from the LL
  if (!(all(sum_probs_1 > 0))) {
    cat("The value of sum_probs_1 is: ", sum_probs_1, "\n")
    stop("problems: sum_probs_1 is non positive!")
  }
  if (!(all(sum_probs_2 > 0))) {
    cat("The value of sum_probs_2 is: ", sum_probs_2, "\n")
    stop("problems: sum_probs_2 is non positive!")
  }
  if (any(is.na(sum_probs_1) | is.nan(sum_probs_1))) {
    cat("The value of sum_probs_1 is: ", sum_probs_1, "\n")
    stop("problems: sum_probs_1 is Na or NaN!")
  }
  if (any(is.na(sum_probs_2) | is.nan(sum_probs_2))) {
    cat("The value of sum_probs_2 is: ", sum_probs_2, "\n")
    stop("problems: sum_probs_2 is Na or NaN!")
  }
}

#' Delivers the likelihood at the end of \link{mbd_loglik}
#' @inheritParams default_params_doc
#' @export
deliver_loglik <- function(
  likelihood,
  sum_probs_1,
  sum_probs_2,
  cond,
  pc
) {

  # Removing sum_probs_1 and sum_probs_2 effects from the LL
  mbd::check_sum_probs(
    sum_probs_1 = sum_probs_1,
    sum_probs_2 = sum_probs_2
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
