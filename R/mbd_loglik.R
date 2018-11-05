#' @author Giovanni Laudanno
#' @title Calculates the likelihood for a multiple birth-death process
#' @description mbd_loglik_basic provides the likelihood for a process
#'   in which multiple births (from different parents) at the same time
#'   are possible, along with usual sympatric speciation and extinction events.
#' @inheritParams default_params_doc
#' @param pars vector of parameters:
#' \itemize{
#'   \item pars[1] is lambda, the sympatric speciation rate;
#'   \item pars[2] is mu, the extinction rate;
#'   \item pars[3] is nu, the multiple allopatric speciation trigger rate;
#'   \item pars[4] is q, the single-lineage speciation probability;
#' }
#' @param safety_threshold It determines the precision on the parameters.
#' @return The function returns the natural logarithm
#'   of the likelihood for the process.
#' @export
mbd_loglik <- function(
  pars,
  brts,
  n_0 = 2,
  cond = 1,
  lx = 1 + 2 * (length(brts) + length(missnumspec)),
  tips_interval = c(0, Inf),
  missnumspec = 0,
  safety_threshold = 1e-3,
  methode = "expo",
  abstol = 1e-16,
  reltol = 1e-10
) {
  #BASIC SETTINGS AND CHECKS
  lambda <- pars[1]
  mu     <- pars[2]
  nu     <- pars[3]
  q      <- pars[4]
  
  if (cond == 0) {
    tips_interval <- c(0, Inf)
  }
  
  if (are_these_parameters_wrong(
    brts = brts,
    pars = pars,
    safety_threshold = safety_threshold,
    n_0 = n_0
  )) {return(-Inf)}
  
  # Use Pure Multiple Birth when there is no extinction
  if (mu == 0 &&
      all(tips_interval == c(0, Inf)) &&
      missnumspec == 0
  ) {
    return(mbd::pmb_loglik(pars = pars, brts = brts, n_0 = n_0))
  }
  # Use Multiple Birth Death model
  #ADJUSTING DATA
  data <- mbd:::brts2time_intervals_and_births(brts) # nolint internal function
  time_intervals <- c(0, data$time_intervals)
  births <- c(0, data$births)
  init_n_lineages <- n_0 #number of starting species
  k_interval <- init_n_lineages + cumsum(births)
  
  pc <- 1
  if (cond == 1) {
    pc <- mbd:::calculate_conditional_probability(
      brts = brts,
      pars = pars,
      n_0 = n_0,
      lx = lx,
      tips_interval = tips_interval,
      methode = methode,
      abstol = abstol,
      reltol = reltol
    )
  }
  
  if (is.nan(log(pc))) {
    # Whatever happens in the future, if pc is NaN,
    # the result will be NaN
    return(NaN)
  }
  
  #LIKELIHOOD INTEGRATION
  
  #SETTING INITIAL CONDITIONS (there's always a +1 because of Q0)
  q_i <- c(1, rep(0, lx))
  q_t <- matrix(0, ncol = (lx + 1), nrow = length(time_intervals))
  q_t[1, ] <- q_i
  dimnames(q_t)[[2]] <- paste0("Q", 0:lx)
  # init_n_lineages is the number of species at t = 1
  k <- init_n_lineages
  # t is starting from 2 so everything is ok with birth[t] and
  # time_intervals[t] vectors
  t <- 2
  D <- C <- rep(1, (length(time_intervals)))
  
  #EVOLVING THE INITIAL STATE TO THE LAST BRANCHING POINT
  while (t <= length(time_intervals)) {
    
    #Applying A operator
    matrix_a <- create_a(
      pars = pars,
      k = k,
      lx = lx
    )
    q_t[t, ] <- mbd:::a_operator(
      q_vector = q_t[(t - 1), ],
      transition_matrix = matrix_a,
      time_interval = time_intervals[t],
      precision = 50L,
      methode = methode,
      abstol = abstol,
      reltol = reltol
    )
    if (methode != "sexpm") {
      # it removes some small negative values that can occur
      # as bugs from the integration process
      q_t[t, ] <- mbd:::negatives_correction(q_t[t, ], pars) # nolint internal function
    }
    if (any(is.nan(q_t[t, ]))) {
      nan_values <- 1
      break
    }
    if (any(q_t[t, ] < 0)) {
      negative_values <- 1
      break
    }
    
    #Applying C operator (this is a trick to avoid precision issues)
    C[t] <- 1 / (sum(q_t[t, ])); q_t[t, ] <- q_t[t, ] * C[t]
    
    if (!(t < length(time_intervals))) {break}
    
    #Applying B operator
    matrix_b <- create_b(
      pars, 
      k = k, 
      b = births[t],
      lx = lx
    )
    q_t[t, ] <- (matrix_b %*% q_t[t, ])
    if (methode != "sexpm") {
      q_t[t, ] <- mbd:::negatives_correction(q_t[t, ], pars) # nolint internal function
    }
    if (any(is.nan(q_t[t, ]))) {
      if (Sys.info()[["sysname"]] == "Windows") {
        print(pars); print(q_t[t, ])
      }
      nan_values <- 1
      break
    }
    if (any(q_t[t, ] < 0)) {
      negative_values <- 1
      break
    }
    
    #Applying D operator (this works exactly like C)
    D[t] <- 1 / (sum(q_t[t, ])); q_t[t, ] <- q_t[t, ] * D[t]
    
    #Updating running parameters
    k <- k + births[t]
    t <- t + 1
  }
  
  #Selecting the state I am interested in
  vm <- 1 / choose((k + missnumspec), k)
  #I have to include +1 because of Q0
  likelihood  <- vm * q_t[t, (missnumspec + 1)]
  
  #Removing C and D effects from the LL
  loglik <- log(likelihood) - sum(log(C)) - sum(log(D))
  
  #Various checks
  loglik <- as.numeric(loglik)
  if (is.nan(loglik) | is.na(loglik)) {
    loglik <- -Inf
  } else {
    loglik <- loglik - log(pc) * (cond == 1) #conditioned likelihood
  }
  loglik
}
