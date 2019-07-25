#' @author Giovanni Laudanno
#' @title Calculates the likelihood for a multiple birth-death process
#' @description mbd_loglik provides the likelihood for a process
#'   in which multiple births (from different parents) at the same time
#'   are possible, along with usual sympatric speciation and extinction events.
#' @inheritParams default_params_doc
#' @return The function returns the natural logarithm
#'   of the likelihood for the process.
#' @export
mbd_loglik <- function(
  pars,
  brts,
  n_0 = 2,
  cond = 1,
  missnumspec = 0,
  lx = 1 + 2 * (length(brts) + length(missnumspec)),
  tips_interval = c(n_0 * (cond > 0), Inf),
  methode = "lsodes",
  q_threshold = 1e-3,
  abstol = 1e-16,
  reltol = 1e-10
) {
  # BASIC SETTINGS AND CHECKS
  check_brts(brts = brts, n_0 = n_0)
  check_cond(cond = cond, tips_interval = tips_interval, n_0 = n_0)
  if (
    check_pars(
      pars = pars,
      safety_checks = TRUE,
      q_threshold = q_threshold
    ) == "wrong"
  ) {
    return(-Inf)
  }

  # Calculate conditional probability
  pc <- calculate_conditional_prob(
    brts = brts,
    pars = pars,
    cond = cond,
    n_0 = n_0,
    lx = lx + 100,
    tips_interval = tips_interval,
    methode = methode,
    abstol = abstol,
    reltol = reltol
  )

  # Use Pure Multiple Birth when there is no extinction
  if (
    pars[2] == 0 &&
    all(tips_interval == c(n_0 * (cond > 0), Inf))
    && missnumspec == 0
  ) {
    return(
      mbd::pmb_loglik(pars = pars, brts = brts, n_0 = n_0) - log(pc)
    )
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
  matrix_a <- vector("list", lt)

  # Evolving the initial state to the present
  while (t <= lt) {

    # Creating A matrix
    matrix_a[[t]] <- create_a(
      pars = pars,
      lx = lx,
      k = k
    )

    # Applying A operator
    q_t[t, ] <- a_operator(
      q_vector = q_t[(t - 1), ],
      transition_matrix = matrix_a[[t]],
      time_interval = time_intervals[t],
      precision = 50L,
      methode = methode,
      abstol = abstol,
      reltol = reltol
    )
    # it removes some small negative values that can occur as bugs from the
    # integration process
    q_t[t, ] <- negatives_correction(q_t[t, ], pars)  # nolint internal function
    if (any(q_t[t, ] < 0)) { # debug
      print(q_t[t, ]) # debug
      plot(q_t[t, ]) # debug
      stop("problems: q_t is negative!") # debug
    } # debug

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
    q_t[t, ] <- negatives_correction(q_t[t, ], pars)  # nolint internal function
    if (any(q_t[t, ] < 0)) { # debug
      print(q_t[t, ]) # debug
      plot(q_t[t, ], xlab = "m", ylab = "Q_m^k(t)") # debug
      stop("problems: q_t is negative!") # debug
    } # debug

    # Applying D operator (this works exactly like C)
    D[t] <- 1 / sum(sort(q_t[t, ]))
    q_t[t, ] <- q_t[t, ] * D[t]

    # Updating running parameters
    k <- k + births[t]
    t <- t + 1
  }
  # testit::assert(all(q_t >= 0)) # q_t has been correctly integrated
  if (any(q_t[t, ] < 0)) { # debug
    print(q_t[t, ]) # debug
    plot(q_t[t, ]) # debug
    stop("problems: q_t is negative!") # debug
  } # debug
  testit::assert(k == n_0 + sum(births)) # k is the number of tips
  testit::assert(t == lt) # t refers to the last time interval

  # Selecting the state I am interested in
  vm <- 1 / choose(k + missnumspec, k)
  likelihood <- vm * q_t[t, (missnumspec + 1)]

  # Removing C and D effects from the LL
  # testit::assert(all(C >= 0))
  if (!(all(C >= 0))) {
    print(pars)
    stop("problems: C is negative!") 
  }
  if (!(all(D >= 0))) {
    print(pars)
    stop("problems: D is negative!") 
  }
  loglik <- log(likelihood) - sum(log(C)) - sum(log(D))

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
