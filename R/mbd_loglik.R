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
  lx = min(1 + 3 * (length(brts) + max(missnumspec)), mbd::max_lx()),
  tips_interval = c(n_0 * (cond > 0), Inf),
  q_threshold = 1e-3,
  debug_mode = FALSE,
  fortran = TRUE
) {

  # BASIC SETTINGS AND CHECKS
  mbd::check_brts(brts = brts, n_0 = n_0)
  mbd::check_cond(cond = cond, tips_interval = tips_interval, n_0 = n_0)
  if (
    mbd::check_pars(
      pars = pars,
      safety_checks = TRUE,
      q_threshold = q_threshold
    ) == "wrong"
  ) {
    return(-Inf)
  }

  # FORTRAN?
  if (fortran == TRUE) {
    rhs_function <- "mbd_runmod"
  } else {
    rhs_function <- mbd::mbd_loglik_rhs
  }

  # Calculate conditional probability
  if (cond == 1) {
    pc <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = mbd::get_lx_condprob(pars = pars, n_0 = n_0, age = max(abs(brts))),
      eq = mbd::condprob_select_eq(
        pars = pars,
        brts = brts,
        fortran = fortran
      ),
      fortran = fortran
    )
  } else {
    pc <- 1
  }

  # Use Pure Multiple Birth when there is no extinction
  if (mbd::is_pbd(
    pars = pars,
    cond = cond,
    n_0 = n_0,
    tips_interval = tips_interval,
    missnumspec = missnumspec
  )) {
    return(mbd::pmb_loglik(pars = pars, brts = brts, n_0 = n_0))
  }

  # Adjusting data
  data <- mbd::brts2time_intervals_and_births(brts) # nolint internal function
  time_intervals <- data$time_intervals
  births <- data$births
  lt <- length(time_intervals)
  testit::assert(n_0 - 1 + length(brts) == n_0 + sum(births)) #every tip is born

  # Setting initial conditions (there's always a +1 because of Q0)
  q_t <- mbd::initialize_q_t(lx = lx, lt = lt)
  k <- n_0 # n_0 is the number of species at t = 1
  t <- 2 # t is starts from 2 so all is ok with births[t] and time_intervals[t]
  sum_probs_2 <- sum_probs_1 <- rep(1, lt)

  # Evolving the initial state to the present
  while (t <= lt) {

    # Creating A matrix
    matrix_a <- mbd::create_a(
      pars = pars,
      lx = lx,
      k = k
    )

    # Applying A operator
    q_t[t, ] <- mbd::mbd_solve(
      vector = q_t[(t - 1), ],
      parms = matrix_a,
      func = rhs_function,
      time_interval = time_intervals[t]
    )
    mbd::check_q_vector(
      q_t = q_t,
      t = t,
      pars = pars,
      brts = brts,
      debug_mode = debug_mode
    )

    # Normalizing the q_vector
    sum_probs_1[t] <- sum(sort(q_t[t, ]))
    q_t[t, ] <- q_t[t, ] / sum_probs_1[t]

    # Loop has to end after integrating to t_p
    if (!(t < lt)) {
      break
    }

    # Creating B matrix
    matrix_b <- mbd::create_b(
      pars = pars,
      lx = lx,
      k = k,
      b = births[t]
    )

    # Applying B operator
    q_t[t, ] <- (matrix_b %*% q_t[t, ])
    mbd::check_q_vector(
      q_t = q_t,
      t = t,
      pars = pars,
      brts = brts,
      debug_mode = debug_mode
    )

    # Normalizing the q_vector
    sum_probs_2[t] <- sum(sort(q_t[t, ]))
    q_t[t, ] <- q_t[t, ] / sum_probs_2[t]

    # Updating running parameters
    k <- k + births[t]
    t <- t + 1
  }

  testit::assert(k == n_0 + sum(births)) # k is the number of tips
  testit::assert(t == lt) # t refers to the last time interval

  # Selecting the state I am interested in
  vm <- 1 / choose(k + missnumspec, k)
  likelihood <- vm * q_t[t, (missnumspec + 1)]

  loglik <- mbd::deliver_loglik(
    likelihood = likelihood,
    sum_probs_1 = sum_probs_1,
    sum_probs_2 = sum_probs_2,
    cond = cond,
    pc = pc,
    debug_mode = debug_mode
  )
  loglik
}
