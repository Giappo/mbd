#' @author Giovanni Laudanno
#' @title Calculates the likelihood for a multiple birth-death process
#' @description mbd_loglik provides the likelihood for a process
#'   in which multiple births (from different parents) at the same time
#'   are possible, along with usual sympatric speciation and extinction events.
#' @inheritParams default_params_doc
#' @return The function returns the natural logarithm
#'   of the likelihood for the process.
#' @export
mbd_loglik_debug <- function(
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
  pc <- calculate_conditional_prob_debug(
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
    if (any(q_t[t, ] < 0)) { # debug
      print(q_t[t, ]) # debug
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
    if (any(q_t[t, ] < 0)) { # debug
      print(q_t[t, ]) # debug
    } # debug

    # Applying D operator (this works exactly like C)
    D[t] <- 1 / sum(sort(q_t[t, ]))
    q_t[t, ] <- q_t[t, ] * D[t]

    # Updating running parameters
    k <- k + births[t]
    t <- t + 1
  }
  # testit::assert(all(q_t >= 0)) # q_t has been correctly integrated
  if (!all(q_t >= 0)) {
    print(pars)
    stop("problems!") #!!!
  }
  testit::assert(k == n_0 + sum(births)) # k is the number of tips
  testit::assert(t == lt) # t refers to the last time interval

  # Selecting the state I am interested in
  vm <- 1 / choose(k + missnumspec, k)
  likelihood <- vm * q_t[t, (missnumspec + 1)]

  # Removing C and D effects from the LL
  # testit::assert(all(C >= 0))
  if (!(all(C >= 0))) {
    print(pars)
    stop("problems!") 
  }
  testit::assert(all(D >= 0))
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

#' @noRd
a_operator_debug <- function(
  q_vector,
  time_interval,
  transition_matrix,
  abstol = 1e-16,
  reltol = 1e-12,
  methode = "lsodes"
) {
  out <- deSolve::ode( # nolint internal function
    y = q_vector,
    times = c(0, time_interval),
    func = mbd_loglik_rhs,
    parms = transition_matrix,
    atol = abstol,
    rtol = reltol,
    method = methode
  )[2, -1]
  out
}

#' Called by \link{mbd_loglik} if there is a conditioning != 0
#' @inheritParams default_params_doc
#' @return the pc. If \code{is.nan(log(pc))} the log-likelihood estimation
#'   by \link{mbd_loglik} is -Inf
#' @author Giovanni Laudanno
#' @noRd
calculate_conditional_prob_debug <- function(
  brts,
  pars,
  cond,
  lx = 1000,
  n_0 = 2,
  tips_interval = c(n_0 * (cond > 0), Inf),
  methode = "lsodes",
  abstol = 1e-16,
  reltol = 1e-10
) {
  check_cond(cond = cond, tips_interval = tips_interval, n_0 = n_0)
  if (cond == 0) {
    return(1)
  }
  pc0 <- pc1
  total_time <- max(abs(brts))
  m <- 0:lx
  one_over_cm <- (3 * (m + 1)) / (m + 3)
  one_over_qm_binom <- 1 / choose(m + n_0, n_0)
  q_i <- c(1, rep(0, lx))
  testit::assert(length(one_over_cm) == length(m))
  testit::assert(length(one_over_qm_binom) == length(m))
  testit::assert(length(q_i) == length(m))
  # creating a_matrix
  matrix_a <- create_a(pars = pars, k = n_0, lx = lx) # nolint internal function
  # integrating the starting q_vector to t_p
  a2_v1 <- a_operator(
    q_vector = q_i,
    transition_matrix = matrix_a,
    time_interval = total_time,
    precision = 250L,
    methode = methode,
    abstol = abstol,
    reltol = reltol
  )
  names(a2_v1) <- paste0("Q", 0:lx)

  total_product <- a2_v1 * one_over_cm * one_over_qm_binom
  missingspecies_min <- max(tips_interval[1] - n_0, 0)
  missingspecies_max <- min(tips_interval[2] - n_0, lx)
  # +1 is because of the zero-th component
  tips_components <- 1 + c(missingspecies_min, missingspecies_max)
  pc1 <- sum(total_product[tips_components[1]:tips_components[2]])

  pc <- pc0 * (cond == 0) + pc1 * (cond == 1)
  if (!((pc >= 0 && pc <= 1))) {
    print(paste0("problem, pc is ", pc)) 
  }
  pc
}

#' @title Maximization of the loglikelihood
#'   under a multiple birth-death diversification model
#' @description mbd_ml computes the maximum likelihood estimates of
#'   the parameters of a multiple birth-death diversification model
#'   for a given set of phylogenetic branching times.
#'   It also outputs the corresponding loglikelihood
#'   that can be used in model comparisons.
#'   Differently from mbd_ml it can account for three kind of events:
#'   sympatric (single) speciation, multiple (allopatric) speciation
#'   and extinction.
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return The output is a dataframe containing
#'   estimated parameters and maximum loglikelihood.
#'   The computed loglikelihood contains the factor q! m! / (q + m)!
#'   where q is the number of species in the phylogeny and m is the number of
#'   missing species, as explained in the supplementary
#'   material to Etienne et al. 2012.
#' @examples
#' set.seed(2)
#' lambda <- 0.2 # sympatric speciation rate
#' mu <- 0.15 # extinction rate;
#' nu <- 2.0 # multiple allopatric speciation trigger rate
#' q <- 0.1 # single-lineage speciation probability
#' sim_pars <- c(lambda, mu, nu, q)
#' crown_age <- 1
#' cond <- 1
#' n_0 <- 2
#' sim <- mbd_sim(
#'  pars = sim_pars,
#'  n_0 = n_0, # Use a crown age
#'  age = crown_age,
#'  cond = cond # Condition on non-extinction
#' )
#' start_pars <- c(0.2, 0.15, 2, 0.15)
#' optim_ids <- c(FALSE, FALSE, FALSE, TRUE)
#' graphics::plot(sim$reconstructed_tree)
#' out <- mbd::mbd_ml(
#'   start_pars = start_pars,
#'   optim_ids = optim_ids,
#'   brts = sim$brts,
#'   cond = cond,
#'   n_0 = n_0,
#'   verbose = FALSE
#' )
#' @export
mbd_ml_debug <- function(
  loglik_function = mbd_loglik_debug,
  brts,
  start_pars = c(0.5, 0.3, 0.5, 0.3),
  n_0 = 2,
  cond = 1,
  optim_ids = rep(TRUE, length(start_pars)),
  true_pars = start_pars,
  tips_interval = c(0, Inf),
  q_threshold = 1e-3,
  verbose = TRUE,
  lx = 1 + 2 * (length(brts)),
  methode = "lsodes"
) {
  # setup and checks
  if (true_pars[3] == 0 | true_pars[4] == 0) {
    q_threshold <- 0
  }
  par_names <- get_param_names() # nolint internal function
  testit::assert(length(optim_ids) == length(start_pars))
  testit::assert(length(true_pars) == length(start_pars))
  start_pars[!optim_ids] <- true_pars[!optim_ids]
  if (any(start_pars < 0)) {
    stop("You cannot start from negative parameters!")
  }
  out_names <- c(par_names, "loglik", "df", "conv")
  failpars <- rep(-1, length(start_pars))
  failout  <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
  colnames(failout) <- out_names

  # define function to optimize
  optim_fun <- function(tr_optim_pars) {
    pars2 <- rep(0, length(start_pars))
    optim_pars <- pars_transform_back(tr_optim_pars) # nolint internal function
    pars2[optim_ids] <- optim_pars
    pars2[!optim_ids] <- true_pars[!optim_ids]

    out <- -loglik_function(
      pars = pars2,
      brts = brts,
      cond = cond,
      n_0 = n_0,
      q_threshold = q_threshold,
      lx = lx,
      methode = methode
    )
    if (verbose == TRUE) {
      printed_values <- paste0(
        c(par_names, "loglik"),
        " = ",
        signif(c(pars2, -out), digits = 5)
      )
      print_this <- paste(printed_values, sep = ",")
      cat(print_this, "\n")
    }
    out
  }

  # initial likelihood
  tr_start_pars <- rep(0, length(start_pars))
  tr_start_pars <- pars_transform_forward(start_pars[optim_ids]) # nolint internal function
  initloglik <- -optim_fun(tr_start_pars)
  utils::flush.console()
  if (initloglik == -Inf) {
    cat(
      message = "The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n" # nolint
    )
    out2 <- failout
    return(invisible(out2))
  }

  # maximum likelihood
  out <- subplex::subplex(
    par = tr_start_pars,
    fn = function(x) optim_fun(x)
  )

  # report missed convergence
  if (out$conv > 0) {
    cat2(
      "Optimization has not converged. Try again with different initial values.\n", # nolint
      verbose = verbose
    )
    out2 <- data.frame(
      t(failpars),
      loglik = -1,
      df = -1,
      conv = unlist(out$conv)
    )
    names(out2) <- out_names
    return(invisible(out2))
  }

  # return mle results
  outpars <- rep(0, length(start_pars))
  outpars[optim_ids] <- pars_transform_back( # nolint internal function
    as.numeric(unlist(out$par))
  )
  outpars[!optim_ids] <- true_pars[!optim_ids]
  names(outpars) <- par_names

  out2 <- data.frame(
    row.names = NULL,
    outpars[1],
    outpars[2],
    outpars[3],
    outpars[4],
    -out$value,
    sum(optim_ids),
    unlist(out$conv)
  )
  names(out2) <- out_names
  return(out2)
}
