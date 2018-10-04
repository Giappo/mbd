#' @author Giovanni Laudanno
#' @title Calculates the likelihood for a multiple birth-death process
#' @description mbd_loglik provides the likelihood for a process
#'   in which multiple births (from different parents) at the same time
#'   are possible, along with usual sympatric speciation and extinction events.
#' @inheritParams default_params_doc
#' @param pars vector of parameters:
#' \itemize{
#'   \item pars[1] is lambda, the sympatric speciation rate;
#'   \item pars[2] is mu, the extinction rate;
#'   \item pars[3] is nu, the multiple allopatric speciation trigger rate;
#'   \item pars[4] is q, the single-lineage speciation probability_
#' }
#' @param safety_threshold It determines the precision on the parameters.
#' @return The function returns the natural logarithm
#'   of the likelihood for the process.
#' @examples
#' set.seed(11)
#' simulated_data = mbd_sim(
#'   pars = c(0.6, 0.1, 2.2, 0.1), soc = 2, age = 10, cond = 1
#' )
#' graphics::plot(simulated_data$tas)
#' # @Giappo: too big too run
#' # mbd::mbd_loglik(
#' #   pars = c(0.8, 0.05, 2.2, 0.1),
#' #   brts = simulated_data$brts,
#' #   soc = 2, cond = 1, missnumspec = 0
#' # )
#'
#' @seealso use \link{pmb_loglik} to calculate the log-likelihood for
#'   a Pure Multiple Birth model.
#' @export
mbd_loglik <- function(
    pars,
    brts,
    soc = 2,
    cond = 1,
    lx0 = 900,
    alpha = 10,
    tips_interval = c(0, Inf),
    missnumspec = 0,
    safety_threshold = 1e-3,
    methode = "expo",
    minimum_multiple_births = 0,
    print_errors = TRUE
) {
  #BASIC SETTINGS AND CHECKS
  lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4]
  abstol <- 1e-16; reltol <- 1e-10
  if (cond == 0) {
    tips_interval <- c(0, Inf)
  }
  if (length(pars) != 4) {
    stop("'pars' must have a length of four")
  }
  if (any(is.nan(pars))) {
    stop("'pars' cannot contain NaNs")
  }
  if (any(is.infinite(pars))) {
    stop("'pars' cannot contain Infs")
  }
  if (lambda < 0.0) {
    stop("'lambda' must be positive")
  }
  if (mu < 0.0) {
    stop("'mu' must be positive")
  }
  if (nu < 0.0) {
    stop("'nu' must be positive")
  }
  if (q < 0.0) {
    stop("'q' must be positive")
  }
  if (q > 1.0) {
    stop("'q' must be less or equal to one")
  }
  if (minimum_multiple_births < 0) {
    stop("'minimum_multiple_births' must be positive")
  }

  testit::assert(q >= 0 - safety_threshold)
  testit::assert(q <= 1.0 + safety_threshold)
  # Use Pure Multiple Birth when there is no extinction
  if (mu == 0 &&
      all(tips_interval == c(0, Inf)) &&
      missnumspec == 0 &&
      minimum_multiple_births == 0
  ) {
    return(mbd::pmb_loglik(pars = pars, brts = brts, soc = soc))
  }
  # Use Multiple Birth Death model
  #ADJUSTING DATA
  data <- brts2time_intervals_and_births(brts) # nolint internal function
  time_intervals <- c(0, data$time_intervals)
  births <- c(0, data$births)
  init_n_lineages <- soc #number of starting species
  k_interval <- init_n_lineages + cumsum(births)
  max_k <- max(k_interval)

  #DETERMINE PC AND ALPHA (OR LX)
  pc_and_alpha <- mbd::alpha_analysis(
    brts = brts,
    pars = pars,
    tips_interval = tips_interval,
    cond = cond,
    soc = soc,
    alpha0 = alpha,
    max_k = max_k,
    methode = methode,
    abstol = abstol,
    reltol = reltol,
    minimum_multiple_births = minimum_multiple_births
  )

  pc    <- pc_and_alpha$pc
  # alpha is the proportionality factor between max_k and
  # the edge of the matrix
  alpha <- pc_and_alpha$alpha
  lx <- alpha * max_k;
  pc <- 1
  if (cond == 1) {
    pc <- mbd::calc_cond_prob(
      brts = brts,
      pars = pars,
      soc = soc,
      lx = lx,
      tips_interval = tips_interval,
      methode = methode,
      abstol = abstol,
      reltol = reltol
    )
  }

  #LIKELIHOOD INTEGRATION
  start_over_again <- 1; iterations <- 0; max_iterations <- 100
  negative_values <- nan_values <- 0;
  while (start_over_again == 1 & iterations < max_iterations) {
    #MATRIX DIMENSION SETUP
    # alpha is the proportionality factor between max_k
    # and the edge of the matrix
    lx <- alpha * max_k;

    #SETTING INITIAL CONDITIONS (there's always a +1 because of Q0)
    q_i <- c(1, rep(0, lx))
    #do I need a +1 in nrow?
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
      transition_matrix <- create_a(
        lambda = lambda,
        mu = mu,
        nu = nu,
        q = q,
        k = k,
        max_number_of_species = lx
      )
      q_t[t, ] <- a_operator(
        q_matrix = q_t[(t - 1), ],
        transition_matrix = transition_matrix,
        time_interval = time_intervals[t],
        precision = 50L,
        methode = methode,
        a_abstol = abstol,
        a_reltol = reltol
      )
      if (methode != "sexpm") {
        # it removes some small negative values that can occur
        # as bugs from the integration process
        q_t[t, ] <- negatives_correction(q_t[t, ], pars) # nolint internal function
      }
      if (any(is.nan(q_t[t, ]))) {
        if (Sys.info()[["sysname"]] == "Windows") {
          print(pars); print(q_t[t, ])
        }
        nan_values <- 1; break
      }
      if (any(q_t[t, ] < 0)) {
        negative_values <- 1
        break
      }

      #Applying C operator (this is a trick to avoid precision issues)
      C[t] <- 1 / (sum(q_t[t, ])); q_t[t, ] <- q_t[t, ] * C[t]

      if (t < length(time_intervals)) {
        #Applying B operator
        matrix_b <- create_b(
          lambda = lambda, nu = nu, q = q, k = k, b = births[t],
          max_number_of_species = lx
        )
        q_t[t, ] <- (matrix_b %*% q_t[t, ])
        if (methode != "sexpm") {
          q_t[t, ] <- negatives_correction(q_t[t, ], pars) # nolint internal function
        }
        if (any(is.nan(q_t[t, ]))) {
          if (Sys.info()[["sysname"]] == "Windows") {
            print(pars); print(q_t[t, ])
          }
          nan_values <- 1; break
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
      } else {
        break
      }
    }
    if (negative_values == 1) {
      alpha <- alpha + 2
    }
    if (nan_values == 1) {
      alpha <- alpha - 5;
      if (alpha <= 0) {
        alpha <- 6
      }
    }
    iterations <- iterations + 1
    start_over_again <- (nan_values) || (negative_values)
  }

  #Selecting the state I am interested in
  vm <- 1 / choose((k + missnumspec), k)
  #I have to include +1 because of Q0
  p_matrix  <- vm * q_t[t, (missnumspec + 1)]

  #Removing C and D effects from the LL
  loglik <- log(p_matrix) - sum(log(C)) - sum(log(D))

  #Various checks
  loglik <- as.numeric(loglik)
  if (is.nan(loglik) | is.na(loglik)) {
    loglik <- -Inf
  } else {
    loglik <- loglik - log(pc) * (cond == 1) #conditioned likelihood
  }
  loglik
}
