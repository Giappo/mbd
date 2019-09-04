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
  lx = min(1 + 2 * (length(brts) + max(missnumspec)), max_lx()),
  tips_interval = c(n_0 * (cond > 0), Inf),
  q_threshold = 1e-3,
  debug_mode = FALSE
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
  pc <- cond_prob(
    pars = pars,
    brts = brts,
    cond = cond,
    n_0 = n_0,
    tips_interval = tips_interval,
    lx = min(lx, 31),
    debug_mode = debug_mode
  )

  # Use Pure Multiple Birth when there is no extinction
  if (
    pars[2] == 0 &&
    all(tips_interval == c(n_0 * (cond > 0), Inf)) &&
    missnumspec == 0
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
      debug_mode = debug_mode
    )
    check_q_vector(
      q_t = q_t,
      t = t,
      pars = pars,
      brts = brts,
      debug_mode = debug_mode
    )

    # Applying C operator (this is a trick to avoid precision issues)
    C[t] <- sum(sort(q_t[t, ]))
    q_t[t, ] <- q_t[t, ] / C[t]

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
    check_q_vector(
      q_t = q_t,
      t = t,
      pars = pars,
      brts = brts,
      debug_mode = debug_mode
    )

    # Applying D operator (this works exactly like C)
    D[t] <- sum(sort(q_t[t, ]))
    q_t[t, ] <- q_t[t, ] / D[t]

    # Updating running parameters
    k <- k + births[t]
    t <- t + 1
  }

  testit::assert(k == n_0 + sum(births)) # k is the number of tips
  testit::assert(t == lt) # t refers to the last time interval

  # Selecting the state I am interested in
  vm <- 1 / choose(k + missnumspec, k)
  likelihood <- vm * q_t[t, (missnumspec + 1)]

  # Removing C and D effects from the LL
  if (!(all(C > 0))) {
    cat("The value of C is: ", C, "\n")
    if (debug_mode == FALSE) {
      stop("problems: C is non positive!")
    }
  }
  if (!(all(D > 0))) {
    cat("The value of D is: ", D, "\n")
    if (debug_mode == FALSE) {
      stop("problems: D is non positive!")
    }
  }
  if (any(is.na(C) || is.nan(C))) {
    cat("The value of C is: ", C, "\n")
    if (debug_mode == FALSE) {
      stop("problems: C is Na or NaN!")
    }
  }
  if (any(is.na(D) || is.nan(D))) {
    cat("The value of D is: ", D, "\n")
    if (debug_mode == FALSE) {
      stop("problems: D is Na or NaN!")
    }
  }
  if (debug_mode == TRUE) {
   cat("The value of the likelihood is: ", likelihood, "\n")
  }
  loglik <- log(likelihood) + sum(log(C)) + sum(log(D))

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
