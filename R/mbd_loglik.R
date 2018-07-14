#' @author Giovanni Laudanno
#' @title Calculates the likelihood for a multiple birth-death process
#' @description mbd_loglik provides the likelihood for a process in which multiple births (from different parents) at the same time are possible, along with usual sympatric speciation and extinction events.
#' @param pars vector of parameters:
#' \itemize{
#'   \item pars[1] is lambda, the sympatric speciation rate;
#'   \item pars[2] is mu, the extinction rate;
#'   \item pars[3] is nu, the multiple allopatric speciation trigger rate;
#'   \item pars[4] is q, the single-lineage speciation probability.
#' }
#' @param brts A set of branching times of a phylogeny.
#' @param soc Sets whether stem or crown age should be used (1 or 2)
#' @param cond Set 1 if you want to condition on stem or crown age and non-extinction of the phylogeny. Set 0 otherwise.
#' @param lx It is the number of ODEs considered for the computation.
#' @param tips_interval It takes into account tips boundaries constrain on simulated dataset.
#' @param missnumspec The number of species that are in the clade but missing in the phylogeny.
#' @param methode Specifies how the integration must be performed: set "sexpm" if you want to use sexpm; set "expo" if you want to use expoRkit; set "lsoda" if you want to use the "lsoda" method with the "deSolve::ode" function.
#' @param safety_threshold It determines the precision on the parameters.
#' @return The function returns the natural logarithm of the likelihood for the process.
#'
#' @examples
#' set.seed(11)
#' simulated_data = MBD:::mbd_sim(pars = c(0.6, 0.1, 2.2, 0.1), soc = 2, age = 10, cond = 1)
#' plot(simulated_data$tas)
#' mbd_loglik(pars = c(0.8, 0.05, 2.2, 0.1), brts = simulated_data$brts, soc = 2, cond = 1, missnumspec = 0)
#'
#' @export
mbd_loglik <- function(pars, 
                       brts, 
                       soc = 2, 
                       cond = 1, 
                       lx0 = 600,
                       tips_interval = c(0, Inf),
                       missnumspec = 0, 
                       safety_threshold = 1e-3,
                       methode = "expo", 
                       minimum_multiple_births = 0, 
                       print_errors = TRUE){
  
  #Optional stuff that I might need to run the program one line at the time:
  #brts = sim_data[[1]]; missnumspec = 0;pars = sim_pars; missing_interval = c(1, Inf); methode = "expo"
  
  #BASIC SETTINGS AND CHECKS
  lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4]
  abstol <- 1e-16; reltol <- 1e-10
  if (cond == 0) {tips_interval <- c(0, Inf)}
  # min_tips <- tips_interval[1]; max_tips <- tips_interval[2]
  
  condition1 <- (any(is.nan(pars)) != 0 | any(is.infinite(pars)) != 0)
  condition2 <- (lambda < 0 | mu < 0 | nu < 0 |
                 q <= 0 + safety_threshold | 
                 q >= 1 - safety_threshold |
                 minimum_multiple_births < 0)
  condition3 <- (length(pars) != 4)
  if       (condition1)
  {
    if (print_errors == TRUE){print("input parameters are either infinite or NaN")}
    loglik <- -Inf
  }else if (condition2)
  {
    if (print_errors == TRUE){print("input parameters have wrong values")}
    loglik <- -Inf
  }else if (condition3)
  {
    if (print_errors == TRUE){print("wrong number of input parameters")}
    loglik <- -Inf
  }else if (mu == 0 && all(tips_interval == c(0, Inf)) && missnumspec == 0 && minimum_multiple_births == 0)
  {
    loglik <- MBD:::pmb_loglik(pars = pars, brts = brts, soc = soc) #using pure birth analytical formula
  }else
  {#MAIN
    
    #ADJUSTING DATA
    data <- MBD:::brts2time_intervals_and_births(brts)
    time_intervals <- c(0, data$time_intervals)
    births <- c(0, data$births)
    N0 <- soc #number of starting species
    k_interval <- N0 + cumsum(births)
    max_k <- max(k_interval)
    lx <- MBD:::determine_k_limit(pars = pars, brts = brts, lx = lx0, soc = soc, 
                                  methode = methode, abstol = abstol, reltol = reltol)
    alpha <- lx/10
    Pc <- 1
    if (cond == 1){
      Pc <- MBD::calculate_conditional_probability(brts = brts, 
                                                   pars = pars,
                                                   soc = soc,
                                                   lx = lx,
                                                   tips_interval = tips_interval,
                                                   methode = methode,
                                                   abstol = abstol,
                                                   reltol = reltol)
    }
    
    #LIKELIHOOD INTEGRATION
    start_over_again <- 1; iterations <- 0; max_iterations <- 100
    negative_values <- nan_values <- 0;
    while (start_over_again == 1 & iterations < max_iterations)
    {
      #MATRIX DIMENSION SETUP
      # max_number_of_species <- alpha * max_k; #alpha is the proportionality factor between max_k and the edge of the matrix
      lx <- 10 * alpha; #bonus line to simplify
      nvec <- 0:lx

      #SETTING INITIAL CONDITIONS (there's always a +1 because of Q0)
      Qi <- c(1, rep(0, lx))
      Qt <- matrix(0, ncol = (lx + 1), nrow = length(time_intervals))
      Qt[1,] <- Qi
      dimnames(Qt)[[2]] <- paste0("Q", 0:lx)
      k <- N0 #N0 is the number of species at t=1
      t <- 2  #t is starting from 2 so everything is ok with birth[t] and time_intervals[t] vectors
      D <- C <- rep(1, (length(time_intervals)))

      #EVOLVING THE INITIAL STATE TO THE LAST BRANCHING POINT
      while (t <= length(time_intervals))
      {
        #Applying A operator
        transition_matrix <- MBD:::create_A(lambda = lambda, mu = mu, nu = nu, q = q, k = k,max_number_of_species = lx)
        Qt[t,] <- MBD:::A_operator(Q = Qt[(t-1),], transition_matrix = transition_matrix, time_interval = time_intervals[t], precision = 50L, methode = methode, A_abstol = abstol, A_reltol = reltol)
        if (methode != "sexpm"){Qt[t,] <- MBD:::negatives_correction(Qt[t,], pars)} #it removes some small negative values that can occurr as bugs from the integration process
        if (any(is.nan(Qt[t,])))
        {
          if (Sys.info()[['sysname']] == "Windows")
          {
            print(pars); print(Qt[t,])
          }
          nan_values <- 1; break
        }
        if (any(Qt[t,] < 0)){negative_values <- 1; break}

        #Applying C operator (this is a trick to avoid precision issues)
        C[t] <- 1/(sum(Qt[t,])); Qt[t,] <- Qt[t,] * C[t]

        if (t < length(time_intervals))
        {
          #Applying B operator
          B <- MBD:::create_B(lambda = lambda, nu = nu, q = q, k = k, b = births[t],
                              max_number_of_species = lx)
          Qt[t,] <- (B %*% Qt[t,])
          if (methode != "sexpm"){Qt[t,] <- MBD:::negatives_correction(Qt[t,], pars)}
          if (any(is.nan(Qt[t,])))
          {
            if (Sys.info()[['sysname']] == "Windows")
            {
              print(pars); print(Qt[t,])
            }
            nan_values <- 1; break
          }
          if (any(Qt[t,] < 0)){negative_values <- 1; break}

          #Applying D operator (this works exactly like C)
          D[t] <- 1/(sum(Qt[t,])); Qt[t,] <- Qt[t,] * D[t]

          #Updating running parameters
          k <- k + births[t]
          t <- t + 1
        }else{break}
      }
      if (negative_values == 1){alpha <- alpha + 2}
      if (nan_values == 1)     {alpha <- alpha - 5; if(alpha <= 0){alpha <- 6}}
      iterations <- iterations + 1
      start_over_again <- (nan_values) || (negative_values)
    }

    #Selecting the state I am interested in
    vm <- 1/choose((k + missnumspec), k)
    P  <- vm * Qt[t, (missnumspec + 1)] #I have to include +1 because of Q0

    #Removing C and D effects from the LL
    loglik <- log(P) - sum(log(C)) - sum(log(D))

    #Various checks
    loglik <- as.numeric(loglik)
    if (is.nan(loglik) | is.na(loglik))
    {
      loglik <- -Inf
    }else
    {
      loglik <- loglik - log(Pc) * (cond == 1) #conditioned likelihood
    }
  }
  return(loglik)
}


