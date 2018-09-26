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

  condition1 <- (any(is.nan(pars)) != 0 | any(is.infinite(pars)) != 0)
  condition2 <- (lambda < 0 | mu < 0 | nu < 0 |
                 q <= 0 + safety_threshold |
                 q >= 1 - safety_threshold |
                 minimum_multiple_births < 0)
  condition3 <- (length(pars) != 4)
  if (condition1) {
    if (print_errors == TRUE) { 
      print("input parameters are either infinite or NaN")
    }
    loglik <- -Inf
  } else if (condition2) {
    if (print_errors == TRUE) { 
      print("input parameters have wrong values")
    }
    loglik <- -Inf
  } else if (condition3) {
    if (print_errors == TRUE) { 
      print("wrong number of input parameters")
    }
    loglik <- -Inf
  } else if (mu == 0 &&
      all(tips_interval == c(0, Inf)) &&
      missnumspec == 0 &&
      minimum_multiple_births == 0
  ) {
    #using pure birth analytical formula
    loglik <- pmb_loglik(pars = pars, brts = brts, soc = soc) 
  } else {
    #MAIN

    #ADJUSTING DATA
    data <- brts2time_intervals_and_births(brts)
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
    lx <- max_number_of_species <- alpha * max_k; 
    pc <- 1
    if (cond == 1) {
      pc <- mbd::calculate_conditional_probability(
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
      lx <- max_number_of_species <- alpha * max_k; 
      # lx <- 10 * alpha; #bonus line to simplify
      # nvec <- 0:lx

      #SETTING INITIAL CONDITIONS (there's always a +1 because of Q0)
      Qi <- c(1, rep(0, lx))
      #do I need a +1 in nrow?
      Qt <- matrix(0, ncol = (lx + 1), nrow = length(time_intervals)) 
      Qt[1,] <- Qi
      dimnames(Qt)[[2]] <- paste0("Q", 0:lx)
      # init_n_lineages is the number of species at t = 1
      k <- init_n_lineages 
      # t is starting from 2 so everything is ok with birth[t] and 
      # time_intervals[t] vectors
      t <- 2  
      D <- C <- rep(1, (length(time_intervals)))

      #EVOLVING THE INITIAL STATE TO THE LAST BRANCHING POINT
      while (t <= length(time_intervals)) {
        #Applying A operator
        transition_matrix <- create_A(
          lambda = lambda, 
          mu = mu, 
          nu = nu, 
          q = q, 
          k = k, 
          max_number_of_species = lx
        )
        Qt[t,] <- A_operator(
          Q = Qt[(t - 1),], 
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
          Qt[t,] <- negatives_correction(Qt[t,], pars)
        } 
        if (any(is.nan(Qt[t,])))
        {
          if (Sys.info()[['sysname']] == "Windows")
          {
            print(pars); print(Qt[t,])
          }
          nan_values <- 1; break
        }
        if (any(Qt[t,] < 0)) {
          negative_values <- 1
          break
        }

        #Applying C operator (this is a trick to avoid precision issues)
        C[t] <- 1 / (sum(Qt[t,])); Qt[t,] <- Qt[t,] * C[t]

        if (t < length(time_intervals))
        {
          #Applying B operator
          B <- create_B(lambda = lambda, nu = nu, q = q, k = k, b = births[t],
                              max_number_of_species = lx)
          Qt[t,] <- (B %*% Qt[t,])
          if (methode != "sexpm") {
            Qt[t,] <- negatives_correction(Qt[t,], pars)
          }
          if (any(is.nan(Qt[t,])))
          {
            if (Sys.info()[['sysname']] == "Windows")
            {
              print(pars); print(Qt[t,])
            }
            nan_values <- 1; break
          }
          if (any(Qt[t,] < 0)) {
            negative_values <- 1 
            break
          }

          #Applying D operator (this works exactly like C)
          D[t] <- 1 / (sum(Qt[t,])); Qt[t,] <- Qt[t,] * D[t]

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
    P  <- vm * Qt[t, (missnumspec + 1)] #I have to include +1 because of Q0

    #Removing C and D effects from the LL
    loglik <- log(P) - sum(log(C)) - sum(log(D))

    #Various checks
    loglik <- as.numeric(loglik)
    if (is.nan(loglik) | is.na(loglik)) {
      loglik <- -Inf
    } else {
      loglik <- loglik - log(pc) * (cond == 1) #conditioned likelihood
    }
  }
  loglik
}
