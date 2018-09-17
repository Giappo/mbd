#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
hyperA_HannoX <- function(N, k, q) {# HG function: fast O(N), updated after Moulis meeting
  #this is the matrix builder: helps to create A and B operators
  #it produces the structure q^(m-n)*(1-q)^(k+2*n-m)*sum_j 2^j choose(k,j)*choose(n,m-n-j)
  j <- 0:k
  A1 <- (1-q)^(k) * choose(k, j) * (2)^j
  N <- N + 1
  A <- diag(A1[1], nrow=N+2, ncol=N+2)
  A[1:(k+1),1] <- A1
  for (dst in 2:N) {
    src <- dst - 1
    s <- src:min(N,2*src+k-1)
    A[s+2,dst] <- A[s,src] + A[s+1,src]
    m <- s-1; n <- src-1;
    A[s,src] <- A[s,src] * q^(m-n)*(1-q)^(2*n-m)
  }
  A[N,N] = A[N,N] * (1-q)^(N-1);
  A[1:N,1:N]
}

# if (.Platform$OS.type=="windows"){matrix_builder = hyperA:::hyperA}
# if (.Platform$OS.type=="unix")   {matrix_builder = hyperA_HannoX} #this is set to work with home directory of the cluster

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
create_A0 <- function(max_number_of_species, lambda, mu, q, k, matrix_builder = hyperA_HannoX){
  nvec = 0:max_number_of_species
  M = lambda * matrix_builder(N = max_number_of_species, k = k, q = q)

  # diag(M) = (-lambda)*c( (1-(1-q)^(k+( nvec[1:max_number_of_species] ) )),0 ) - mu*(nvec+k)
  diag(M) = (-lambda) * ( 1 - (1 - q)^(k + nvec) ) - mu * (nvec + k) #new version to avoid the dumpster problem at the end of the matrix

  M[row(M) == col(M) - 1] = mu * nvec[2:(max_number_of_species + 1)]
  return(M)
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
create_B0 <- function(max_number_of_species, q, k, b, matrix_builder = hyperA_HannoX){#lambda * choose(k,b) * q^b  is going to be added in logB in the main script
  k2 <- k - b
  B <- matrix_builder(N = max_number_of_species, k = k2, q = q)
  return(B)
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
create_A <- function(lambda, mu, nu, q, k, max_number_of_species){
  nvec <- 0:max_number_of_species
  M <- create_A0(max_number_of_species = max_number_of_species,
                       lambda = nu, mu = mu, q = q, k = k)
  M[row(M) == col(M) + 1] <- M[row(M) == col(M) + 1] + lambda * (nvec[1:(max_number_of_species)]+2*k)

  # M[row(M) == col(M)] = M[row(M) == col(M)] - c(lambda*(nvec[-length(nvec)]+k),0)
  M[row(M) == col(M)] <- M[row(M) == col(M)] - lambda * (nvec + k) #new version to avoid the dumpster problem at the end of the matrix

  return(M)
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
create_B <- function(lambda, nu, q, k, b, max_number_of_species){
  M <- create_B0(max_number_of_species = max_number_of_species, q = q, k = k, b = b)
  B <- lambda * k * diag(max_number_of_species + 1) * (b == 1) + nu * choose(k, b) * (q^b) * M
  return(B)
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
create_A.no_mbd = function(lambda,mu,nu,q,k,max_number_of_species,minimum_multiple_births){
  nvec <- 0:max_number_of_species
  M <- create_A0(max_number_of_species = max_number_of_species,lambda = nu,mu = mu,q = q,k = k)
  M[row(M) == col(M) + 1] <- M[row(M) == col(M) + 1] + lambda * (nvec[1:(max_number_of_species)]+2*k)
  M[row(M) == col(M)]     <- M[row(M) == col(M)] - c(lambda * (nvec[-length(nvec)] + k), 0)
  M[row(M) > col(M) + minimum_multiple_births] <- 0
  return(M)
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
create_B.no_mbd = function(lambda,nu,q,k,b,max_number_of_species,minimum_multiple_births){
  M <- create_B0(max_number_of_species = max_number_of_species, q = q, k = k, b = b)
  B <- lambda * k * diag(max_number_of_species + 1) * (b == 1) + nu * choose(k,b) * (q^b) * M
  B[row(B) > col(B) + minimum_multiple_births] = 0
  return(B)
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
A_operator <- function(Q, transition_matrix, time_interval, precision = 50L,
                       A_abstol = 1e-16, A_reltol = 1e-10, methode = "expo"){
  
  precision_limit <- 2000
  precision_step1 <- 40
  precision_step2 <- 50
  max_repetitions <- 10
  result <- rep(-1, length(Q))
  bad_result <- 0
  
  if (methode == "sexpm")
  {
    exp_matrix <- rsexpm:::sexpm(transition_matrix * time_interval)
    result     <- exp_matrix %*% Q
  }
  
  if (methode == "expo")
  {
    result.nan <- result.negative <- 1
    repetition <- 1
    while ((result.nan == 1| result.negative == 1) & repetition < max_repetitions)
    {
      result <- try(expoRkit::expv(v = Q, x = transition_matrix, t = time_interval, m = precision), silent = TRUE)
      
      result.nan <- (any(!is.numeric(result)) || any(is.nan(result)))
      if (result.nan) 
      {
        precision <- precision - precision_step1
      }else
      {
        result.negative <- (any(result < 0))
        if (result.negative) 
        {
          precision <- precision + precision_step2
          if (precision > precision_limit) 
          {
            break
          }
        }
      }
      repetition <- repetition + 1
    }
  }
  
  bad_result <- (any(!is.numeric(result)) || any(is.nan(result)))
  if (!bad_result) {bad_result <- (any(result < 0))}
  
  if (methode == "lsoda" | bad_result)
  {
    times <- c(0, time_interval)
    ode_matrix <- transition_matrix
    R.utils::withTimeout(result <- deSolve::ode(y = Q, times = times, func = mbd_loglik_rhs, parms = ode_matrix, atol = A_abstol, rtol = A_reltol)[2,-1], timeout = 1001)
  }
  
  return(result)
}

#' #' @title Internal mbd function
#' #' @description Internal mbd function.
#' #' @details This is not to be called by the user.
#' #' @export
#' A_operator_old <- function(Q, transition_matrix, time_interval, precision = 50L,
#'                        A_abstol = 1e-16, A_reltol = 1e-10, methode = "expo"){
#' 
#'   precision_limit <- 3000
#' 
#'   if (methode == "sexpm")
#'   {
#'     exp_matrix <- rsexpm:::sexpm(transition_matrix * time_interval)
#'     result     <- exp_matrix %*% Q
#'   }else if (methode == "expo")
#'   {
#'     result <- try(expoRkit::expv(v = Q, x = transition_matrix, t = time_interval, m = precision), silent = T)
#'     while ( ( any(!is.numeric(result)) || any(is.nan(result)) ) && precision < precision_limit )
#'     {
#'       precision <- precision + 200
#'       result <- try(expoRkit::expv(v = Q, x = transition_matrix, t = time_interval, m = precision), silent = T)}
#'     if (!any(is.nan(result)))
#'     {
#'       while (any(result < 0) && precision < precision_limit)
#'       {
#'         precision <- precision + 200
#'         result <- try(expoRkit::expv(v = Q, x = transition_matrix, t = time_interval, m = precision), silent = T)
#'       }
#'     }
#'   }else if (methode=="lsoda")
#'   {
#'     times <- c(0, time_interval)
#'     ode_matrix <- transition_matrix
#'     # result<-deSolve::ode(y = Q, times = times, func = mbd_loglik_rhs, parms = ode_matrix,atol=A_abstol,rtol=A_reltol)[2,-1]
#'     R.utils::withTimeout(result <- deSolve::ode(y = Q, times = times, func = mbd_loglik_rhs, parms = ode_matrix,atol=A_abstol,rtol=A_reltol)[2,-1], timeout = 1000)
#'   }
#' 
#'   if (any(!is.numeric(result)) || any(is.nan(result))) #sometimes expoRkit gives weird negative values. In this case perform standard lsoda integration.
#'   {
#'     times <- c(0, time_interval)
#'     ode_matrix <- transition_matrix
#'     # result<-deSolve::ode(y = Q, times = times, func = mbd_loglik_rhs, parms = ode_matrix,atol=A_abstol,rtol=A_reltol)[2,-1]
#'     R.utils::withTimeout(result <- deSolve::ode(y = Q, times = times, func = mbd_loglik_rhs, parms = ode_matrix,atol=A_abstol,rtol=A_reltol)[2,-1], timeout = 1000)
#'   }
#'   else if (any(result < 0)) #sometimes expoRkit gives weird negative values. In this case perform standard lsoda integration.
#'   {
#'     times <- c(0, time_interval)
#'     ode_matrix <- transition_matrix
#'     # result=deSolve::ode(y = Q, times = times, func = mbd_loglik_rhs, parms = ode_matrix,atol=A_abstol,rtol=A_reltol)[2,-1]
#'     R.utils::withTimeout(result <- deSolve::ode(y = Q, times = times, func = mbd_loglik_rhs, parms = ode_matrix,atol=A_abstol,rtol=A_reltol)[2,-1], timeout = 1000)
#'   }
#' 
#'   # if ( ( any(!is.numeric(result)) || any(is.nan(result)) ) && methode!="sexpm"){
#'   #     # to run this you actually need sexpm to be installed which, on the cluster, might not be the case.
#'   #     exp_matrix = rsexpm:::sexpm(transition_matrix*time_interval)
#'   #     result = exp_matrix %*% Q
#'   # }
#' 
#'   return(result)
#' }

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
mbd_loglik_rhs <- function (t, x, pars){
  #builds right hand side of the ODE set for multiple birth model
  with(as.list(x), {
    starting_vector = x
    transition_matrix = pars
    dx  = rep(0, length(starting_vector))
    dx  = drop(transition_matrix %*% starting_vector)
    out = (dx)
    names(out) = names(x)
    return(list(out))
  })
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
determine_k_limit <- function(pars, brts, lx, soc, methode, abstol = 1e-16, reltol = 1e-10) {
  lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4]
  mvec <- 0:lx
  Qi <- c(1, rep(0, lx))
  total_time <- max(abs(brts));
  T0 <- create_A(lambda = lambda, mu = 0, nu = nu, q = q, k = soc,
                       max_number_of_species = lx); #dim(TM); max(is.na(TM)); max(is.infinite(TM))
  Pm <- A_operator(Q = Qi, transition_matrix = T0, time_interval = total_time,
                         precision = 250L, methode = methode, A_abstol = abstol, A_reltol = reltol)
  # graphics::plot((Pm/sum(Pm)))
  k_limit <- soc + max(mvec[(mvec %in% which((cumsum(Pm/sum(Pm))) <= 0.95))]); k_limit
  return(k_limit)
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
calculate_conditional_probability <- function (brts,
                                               pars,
                                               lx = 1000,
                                               soc = 2,
                                               tips_interval = c(0, Inf),
                                               methode = 'expo',
                                               abstol = 1e-16,
                                               reltol = 1e-10){
  
  lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4];
  total_time <- max(abs(brts));
  
  m <- 0:lx; length(m)
  one_over_Cm <- (3 * (m + 1))/(m + 3); length(one_over_Cm)
  one_over_qm_binom <- 1/choose((m + soc), soc); length(one_over_qm_binom)
  # Qi <- c(1, rep(0, lx)); length(Qi)
  Qi <- rep(0, lx + 1);  Qi[3] <- 1 #starting with k = 0 and m = 2 missing species
  k <- 0 #assuming 0 species
  
  TM <- create_A(lambda = lambda, mu = mu, nu = nu, q = q, k = 0,
                       max_number_of_species = lx); #dim(TM); max(is.na(TM)); max(is.infinite(TM))
  
  A2_v1 <- A_operator(Q = Qi, transition_matrix = TM, time_interval = total_time,
                            precision = 250L, methode = methode, A_abstol = abstol, A_reltol = reltol); A2_v1
  
  # A2_v1 <- try(expoRkit::expv(v = Qi, x = TM, t = total_time, m = 50L), silent = T)
  
  total_product <- A2_v1 * one_over_Cm * one_over_qm_binom
  tips_components <- 1 + 0:1 #these are the components I want to exclude (the one corresponding to 0 and 1 tips)
  Pc <- 1 - sum(total_product[tips_components]); Pc
  # Pc <- sum(total_product)
  
  return(Pc)
}


#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
calculate_conditional_probability0 <- function (brts,
                                                pars,
                                                lx = 1000,
                                                soc = 2,
                                                tips_interval = c(0, Inf),
                                                methode = 'expo',
                                                abstol = 1e-16,
                                                reltol = 1e-10){

  lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4];
  total_time <- max(abs(brts));

  m <- 0:lx; length(m)
  one_over_Cm <- (3 * (m + 1))/(m + 3); length(one_over_Cm)
  one_over_qm_binom <- 1/choose((m + soc), soc); length(one_over_qm_binom)
  Qi <- c(1, rep(0, lx)); length(Qi)

  TM <- create_A(lambda = lambda, mu = mu, nu = nu, q = q, k = soc,
                       max_number_of_species = lx); #dim(TM); max(is.na(TM)); max(is.infinite(TM))

  A2_v1 <- A_operator(Q = Qi, transition_matrix = TM, time_interval = total_time,
  precision = 250L, methode = methode, A_abstol = abstol, A_reltol = reltol)

  # A2_v1 <- try(expoRkit::expv(v = Qi, x = TM, t = total_time, m = 50L), silent = T)

  total_product <- A2_v1 * one_over_Cm * one_over_qm_binom
  missingspecies_min <- max((tips_interval[1] - 2), 0 )
  missingspecies_max <- min((tips_interval[2] - 2), lx)
  tips_components <- 1 + c(missingspecies_min, missingspecies_max) # +1 is because of the zero-th component
  Pc <- sum(total_product[tips_components[1]:tips_components[2]])
  # Pc <- sum(total_product)

  return(Pc)
}


#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
calculate_conditional_probability0PB <- function (brts,
                                                  pars,
                                                  lx = 200,
                                                  soc = 2,
                                                  tips_interval = c(0, Inf),
                                                  methode = 'expo',
                                                  abstol = 1e-16,
                                                  reltol = 1e-10){
  
  lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4];
  total_time <- max(abs(brts))
  if (mu != 0){cat('mu is supposed to be equal zero to use this function'); return(Pc <- NA)}
  
  m <- 0:lx; length(m)
  one_over_Cm <- (3 * (m + 1))/(m + 3); length(one_over_Cm)
  one_over_qm_binom <- 1/choose((m + soc), soc); length(one_over_qm_binom)
  Qi <- c(1, rep(0, lx)); length(Qi)
  
  TM <- create_A(lambda = lambda, mu = mu, nu = nu, q = q, k = soc,
                       max_number_of_species = lx); #dim(TM); max(is.na(TM)); max(is.infinite(TM))
  
  A2_v1 <- A_operator(Q = Qi, transition_matrix = TM, time_interval = total_time,
                            precision = 250L, methode = methode, A_abstol = abstol, A_reltol = reltol)
  
  # A2_v1 <- try(expoRkit::expv(v = Qi, x = TM, t = total_time, m = 50L), silent = T)
  
  total_product <- A2_v1 * one_over_Cm * one_over_qm_binom
  missingspecies_min <- max((tips_interval[1] - 2), 0 )
  missingspecies_max <- min((tips_interval[2] - 2), lx)
  tips_components <- 1 + c(missingspecies_min, missingspecies_max) # +1 is because of the zero-th component
  opposite_tips_components <- (m + 1)[!((m + 1) %in% (tips_components[1]:tips_components[2]))]
  # Pc <- sum(total_product[tips_components[1]:tips_components[2]])
  Pc <- 1 - sum(total_product[opposite_tips_components])
  # Pc <- sum(total_product)
  
  return(Pc)
}



#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
find_best_lx_for_Pc <- function(brts, 
                                pars,
                                soc = 2,
                                methode = 'expo',
                                abstol = 1e-16,
                                reltol = 1e-10,
                                iterations = 20,
                                interval.min = 500,
                                interval.max = 1400) {
  
  a <- iterations/2
  interval.width <- interval.max - interval.min
  step1 <- floor(interval.width/a)
  
  lx.test <- rep(NA, length(lxvec <- seq(interval.min + step1, interval.max - step1, step1))); i <- 1; right.lx.coord <- 0
  for (lx2 in lxvec)
  {
    lx.test[i] <- mbd::calculate_conditional_probability0(brts = brts,
                                                          pars = c(pars[1], 0, pars[3], pars[4]),
                                                          lx = lx2, 
                                                          soc = soc,
                                                          tips_interval = c(0, Inf),
                                                          methode = methode,
                                                          abstol = abstol,
                                                          reltol = reltol)
    if (!is.na(abs(lx.test[i]))) if (abs(lx.test[i] - 1) < 0.01) {right.lx.coord <- i; lx <- lxvec[right.lx.coord]; break}
    i <- i + 1
  }; lx.test
  if (right.lx.coord == 0)
  {
    right.lx.coord <- which(abs(lx.test - 1) == min(abs(lx.test - 1), na.rm = TRUE))
    lx <- lxvec[right.lx.coord]
  }
  
  lx.test2 <- rep(NA, length(lxvec2 <- floor(seq(lx - step1, lx + step1, 2 * step1/a)))); j <- 1; right.lx.coord2 <- 0
  for (lx2 in lxvec)
  {
    lx.test2[j] <- mbd::calculate_conditional_probability0(brts = brts,
                                                           pars = c(pars[1], 0, pars[3], pars[4]),
                                                           lx = lx2, 
                                                           soc = soc,
                                                           tips_interval = c(0, Inf),
                                                           methode = methode,
                                                           abstol = abstol,
                                                           reltol = reltol)
    if (!is.na(abs(lx.test2[i]))) if (abs(lx.test2[j] - 1) < 0.01) {right.lx.coord2 <- j; lx <- lxvec2[right.lx.coord2]; break}
    j <- j + 1
  }; lx.test2
  if (right.lx.coord2 == 0)
  {
    right.lx.coord2 <- which(abs(lx.test2 - 1) == min(abs(lx.test2 - 1), na.rm = T))
    lx <- lxvec2[right.lx.coord2]
  }
  
  return(lx)
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
calculate_conditional_probability1 <- function (brts,
                                               pars,
                                               soc = 2,
                                               tips_interval = c(0, Inf),
                                               methode = 'expo',
                                               abstol = 1e-16,
                                               reltol = 1e-10){
  
  lx <- find_best_lx_for_Pc(brts = brts, pars = pars, soc = soc)
  if (pars[2] == 0)
  {
    Pc <- mbd::calculate_conditional_probability0PB(brts = brts,
                                                    pars = pars,
                                                    lx = lx,
                                                    soc = soc,
                                                    tips_interval = tips_interval,
                                                    methode = methode,
                                                    abstol = abstol,
                                                    reltol = reltol)  
  }else
  {
    # @Giappo: this function does not exist any more
    # Pc <- mbd::calculate_conditional_probability02(brts = brts,
    #                                                pars = pars,
    #                                                lx = lx,
    #                                                soc = soc,
    #                                                tips_interval = tips_interval,
    #                                                methode = methode,
    #                                                abstol = abstol,
    #                                                reltol = reltol)
    stop("Function 'mbd::calculate_conditional_probability02' is absent")
  }
  return(list(Pc = Pc, lx = lx))
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
alpha_conditional_probability <- function (brts, pars, alpha, tips_interval = c(0, Inf),
                                           cond = 1, soc = 2, methode = "expo",
                                           abstol = 1e-16, reltol = 1e-10,
                                           minimum_multiple_births = 0){
  
  lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4];
  min_tips <- tips_interval[1]; max_tips <- tips_interval[2];
  min_tips <- max(min_tips, soc * cond) #check this
  N0 <- soc
  total_time <- max(abs(brts));
  births <- c(0, brts2time_intervals_and_births(brts)$births)
  k_interval <- N0 + cumsum(births)
  max_k <- max(k_interval)
  max_number_of_species <- alpha * max_k; #alpha is the proportionality factor between max_k and the edge of the matrix
  
  if (!(cond == 1 | tips_interval[1] > 0 | tips_interval[2] < Inf))
  {
    Pc <- 1; A2_v1 <- c(1, rep(0, max_number_of_species))
  }else
  {
    m <- 0:max_number_of_species;
    one_over_Cm <- (3 * (m + 1))/(m + 3)
    one_over_qm_binom <- 1/choose((m + N0), N0)
    tips_components <- (1 + min_tips):(1 + min(max_tips, max_number_of_species)) #applying tips constrain
    if (cond == 1){tips_components <- tips_components - N0} #I am already considering the starting species to survive. I must not double count them!
    
    Qi <- c(1, rep(0, max_number_of_species))
    Mk_N0 <- create_A(lambda = lambda, mu = mu, nu = nu, q = q, k = soc,
                            max_number_of_species = max_number_of_species)
    A2_v1 <- A_operator(Q = Qi, transition_matrix = Mk_N0, time_interval = total_time,
                              precision = 50L, methode = methode, A_abstol = abstol, A_reltol = reltol)
    if (methode != "sexpm"){A2_v1 <- negatives_correction(A2_v1, pars)} #it removes some small negative values that can occurr as bugs from the integration process
    
    if (minimum_multiple_births > 0) #adjust for the required minimum amount of mbd
    {
      Mk_N0.no_mbd <- create_A.no_mbd(lambda = lambda, mu = mu, nu = nu, q = q, k = soc, max_number_of_species = max_number_of_species, minimum_multiple_births = minimum_multiple_births)
      A2_v1.no_mbd <- A_operator(Q = Qi, transition_matrix = Mk_N0.no_mbd, time_interval = total_time,precision = 50L,methode=methode,A_abstol=abstol,A_reltol=reltol)
      A2_v1 <- A2_v1 - A2_v1.no_mbd
    }
    
    total_product <- A2_v1 * one_over_Cm * one_over_qm_binom
    Pc <- sum(total_product[tips_components])
  }
  return(list(Pc = Pc, A2_v1 = A2_v1))
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
alpha_analysis <- function(brts,
                           pars,
                           tips_interval,
                           cond,
                           soc,
                           alpha0,
                           max_k,
                           methode = 'expo',
                           abstol,
                           reltol,
                           minimum_multiple_births){
  deltaAlpha <- 1; count <- 0; same_result_count <- 0; Pc.notanumber <- 1;
  alpha <- alpha0
  while (Pc.notanumber)
  {
    Pc1 <- alpha_conditional_probability(brts = brts,
                                               pars = pars,
                                               tips_interval = tips_interval,
                                               cond = cond,
                                               soc = soc,
                                               alpha = alpha,
                                               methode = methode,
                                               abstol = abstol,
                                               reltol = reltol,
                                               minimum_multiple_births = minimum_multiple_births)$Pc
    # Pc1 <- calculate_conditional_probability(alpha = alpha, ...)$Pc
    Pc.notanumber <- is.nan(Pc1)
    alpha <- alpha - Pc.notanumber
  }
  while (deltaAlpha != 0 && count < 100 && same_result_count < 5)
  {
    Pc2 <- alpha_conditional_probability(brts = brts,
                                               pars = pars,
                                               tips_interval = tips_interval,
                                               cond = cond,
                                               soc = soc,
                                               alpha = alpha + deltaAlpha,
                                               methode = methode,
                                               abstol = abstol,
                                               reltol = reltol,
                                               minimum_multiple_births = minimum_multiple_births)$Pc
    # Pc2 <- calculate_conditional_probability(alpha = alpha + deltaAlpha, ...)$Pc
    if (is.nan(Pc2))
    {
      deltaAlpha <- deltaAlpha - 1
    }else if (Pc2 < Pc1)
    {
      deltaAlpha <- deltaAlpha + 1
      same_result_count <- same_result_count + 1
    }else
    {
      same_result_count <- 0
      deltaPc <- abs(Pc2 - Pc1)/Pc1;
      # deltaAlpha = floor(  10*(-1 + 2/( 1 + exp(-(1/2*deltaPc)) ))  )
      deltaAlpha <- floor(10 * deltaPc)
      alpha <- alpha + deltaAlpha
      Pc1 <- Pc2
    }
    
    count = count + 1
    # print(alpha)
  }
  if (max_k * alpha >= 2000)
  {#check to see whether alpha is too big to be handled without memory issues
    alpha <- floor(1500/max_k);
    Pc1 <- alpha_conditional_probability(brts = brts,
                                               pars = pars,
                                               tips_interval = tips_interval,
                                               cond = cond,
                                               soc = soc,
                                               alpha = alpha,
                                               methode = methode,
                                               abstol = abstol,
                                               reltol = reltol,
                                               minimum_multiple_births = minimum_multiple_births)$Pc
    # Pc1 <- calculate_conditional_probability(alpha = alpha, ...)$Pc
  }
  Pc <- Pc1
  if (count >= 100){alpha <- 10}
  if (Pc <= 0 | Pc == Inf | Pc == -Inf){Pc <- 1; print("there's a problem with Pc")}
  return(list(Pc = Pc, alpha = alpha))
}

#' ##' @title Internal mbd function
#' #' @description Internal mbd function.
#' #' @details This is not to be called by the user.
#' #' @export
#' calculate_conditional_probability <- function (brts, pars, tips_interval = c(0, Inf),
#'                                                cond = 1, soc = 2, alpha, methode = "expo",
#'                                                abstol = 1e-16, reltol = 1e-10,
#'                                                minimum_multiple_births = 0){
#'   
#'   lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4];
#'   min_tips <- tips_interval[1]; max_tips <- tips_interval[2];
#'   min_tips <- max(min_tips, soc * cond) #check this
#'   N0 <- soc
#'   total_time <- max(abs(brts));
#'   births <- c(0, brts2time_intervals_and_births(brts)$births)
#'   k_interval <- N0 + cumsum(births)
#'   max_k <- max(k_interval)
#'   max_number_of_species <- alpha * max_k; #alpha is the proportionality factor between max_k and the edge of the matrix
#'   
#'   if (!(cond == 1 | tips_interval[1] > 0 | tips_interval[2] < Inf))
#'   {
#'     Pc <- 1; A2_v1 <- c(1, rep(0, max_number_of_species))
#'   }else
#'   {
#'     m <- 0:max_number_of_species;
#'     one_over_Cm <- (3 * (m + 1))/(m + 3)
#'     one_over_qm_binom <- 1/choose((m + N0), N0)
#'     tips_components <- (1 + min_tips):(1 + min(max_tips, max_number_of_species)) #applying tips constrain
#'     if (cond == 1){tips_components <- tips_components - N0} #I am already considering the starting species to survive. I must not double count them!
#'     
#'     Qi <- c(1, rep(0, max_number_of_species))
#'     Mk_N0 <- create_A(lambda = lambda, mu = mu, nu = nu, q = q, k = soc,
#'                             max_number_of_species = max_number_of_species)
#'     A2_v1 <- A_operator(Q = Qi, transition_matrix = Mk_N0, time_interval = total_time,
#'                               precision = 50L, methode = methode, A_abstol = abstol, A_reltol = reltol)
#'     if (methode != "sexpm"){A2_v1 <- negatives_correction(A2_v1, pars)} #it removes some small negative values that can occurr as bugs from the integration process
#'     
#'     if (minimum_multiple_births > 0) #adjust for the required minimum amount of mbd
#'     {
#'       Mk_N0.no_mbd <- create_A.no_mbd(lambda = lambda, mu = mu, nu = nu, q = q, k = soc, max_number_of_species = max_number_of_species, minimum_multiple_births = minimum_multiple_births)
#'       A2_v1.no_mbd <- A_operator(Q = Qi, transition_matrix = Mk_N0.no_mbd, time_interval = total_time,precision = 50L,methode=methode,A_abstol=abstol,A_reltol=reltol)
#'       A2_v1 <- A2_v1 - A2_v1.no_mbd
#'     }
#'     
#'     total_product <- A2_v1 * one_over_Cm * one_over_qm_binom
#'     Pc <- sum(total_product[tips_components])
#'   }
#'   return(list(Pc = Pc, A2_v1 = A2_v1))
#' }
#' 
#' @title Internal mbd function
#' @description Internal mbd function.
#' @details This is not to be called by the user.
#' @export
#' calculate_conditional_probability2 <- function (brts, pars, missing_tips_interval = c(0, Inf),
#'                                                 soc = 2, alpha, methode = "expo",
#'                                                 abstol = 1e-16, reltol = 1e-10,
#'                                                 minimum_multiple_births = 0){
#'   
#'   lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4];
#'   min_tips <- missing_tips_interval[1]; max_tips <- missing_tips_interval[2];
#'   N0 <- soc
#'   total_time <- max(abs(brts))
#'   births <- c(0, brts2time_intervals_and_births(brts)$births)
#'   k_interval <- N0 + cumsum(births)
#'   max_k <- max(k_interval)
#'   max_number_of_species <- alpha * max_k; #alpha is the proportionality factor between max_k and the edge of the matrix
#'   
#'   m <- 0:max_number_of_species
#'   one_over_Cm <- (3 * (m + 1))/(m + 3)
#'   one_over_qm_binom <- 1/choose((m + N0), N0)
#'   tips_components <- (1 + min_tips):(1 + min(max_tips, max_number_of_species)) #applying tips constrain
#'   
#'   Qi    <- c(1, rep(0, max_number_of_species))
#'   Mk_N0 <- create_A(lambda = lambda, mu = mu, nu = nu, q = q, k = soc,
#'                           max_number_of_species = max_number_of_species)
#'   A2_v1 <- A_operator(Q = Qi, transition_matrix = Mk_N0, time_interval = total_time,
#'                             precision = 50L, methode = methode, A_abstol = abstol, A_reltol = reltol)
#'   
#'   # mbd::mbd_loglik(pars = pars, brts = c(total_time), soc = soc, cond = 0, missnumspec = 0)
#'   
#'   total_product <- A2_v1 * one_over_Cm * one_over_qm_binom
#'   Pc <- sum(total_product[tips_components]); Pc
#'   
#'   
#'   if (min_tips < soc)
#'     min_tips <- max(min_tips, soc * cond) #check this
#'   N0 <- soc
#'   total_time <- max(abs(brts));
#'   births <- c(0, brts2time_intervals_and_births(brts)$births)
#'   k_interval <- N0 + cumsum(births)
#'   max_k <- max(k_interval)
#'   max_number_of_species <- alpha * max_k; #alpha is the proportionality factor between max_k and the edge of the matrix
#'   
#'   if (!(cond == 1 | tips_interval[1] > 0 | tips_interval[2] < Inf))
#'   {
#'     Pc <- 1; A2_v1 <- c(1, rep(0, max_number_of_species))
#'   }else
#'   {
#'     m <- 0:max_number_of_species;
#'     one_over_Cm <- (3 * (m + 1))/(m + 3)
#'     one_over_qm_binom <- 1/choose((m + N0), N0)
#'     tips_components <- (1 + min_tips):(1 + min(max_tips, max_number_of_species)) #applying tips constrain
#'     if (cond == 1){tips_components <- tips_components - N0} #I am already considering the starting species to survive. I must not double count them!
#'     
#'     Qi <- c(1, rep(0, max_number_of_species))
#'     Mk_N0 <- create_A(lambda = lambda, mu = mu, nu = nu, q = q, k = soc,
#'                             max_number_of_species = max_number_of_species)
#'     A2_v1 <- A_operator(Q = Qi, transition_matrix = Mk_N0, time_interval = total_time,
#'                               precision = 50L, methode = methode, A_abstol = abstol, A_reltol = reltol)
#'     if (methode != "sexpm"){A2_v1 <- negatives_correction(A2_v1, pars)} #it removes some small negative values that can occurr as bugs from the integration process
#'     
#'     if (minimum_multiple_births > 0) #adjust for the required minimum amount of mbd
#'     {
#'       Mk_N0.no_mbd <- create_A.no_mbd(lambda = lambda, mu = mu, nu = nu, q = q, k = soc, max_number_of_species = max_number_of_species, minimum_multiple_births = minimum_multiple_births)
#'       A2_v1.no_mbd <- A_operator(Q = Qi, transition_matrix = Mk_N0.no_mbd, time_interval = total_time,precision = 50L,methode=methode,A_abstol=abstol,A_reltol=reltol)
#'       A2_v1 <- A2_v1 - A2_v1.no_mbd
#'     }
#'     
#'     total_product <- A2_v1 * one_over_Cm * one_over_qm_binom
#'     Pc <- sum(total_product[tips_components])
#'   }
#'   return(list(Pc = Pc, A2_v1 = A2_v1))
#' }
#' 
#' #' @title Internal mbd function
#' #' @description Internal mbd function.
#' #' @details This is not to be called by the user.
#' #' @export
#' calculate_conditional_probability3 <- function (brts, 
#'                                                 pars, 
#'                                                 lx = 200,
#'                                                 soc = 2){
#'   
#'   lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4];
#'   total_time <- max(abs(brts));
#'   
#'   m <- 0:lx; length(m)
#'   one_over_Cm <- (3 * (m + 1))/(m + 3); length(one_over_Cm)
#'   one_over_qm_binom <- 1/choose((m + soc), soc); length(one_over_qm_binom)
#'   Qi <- c(1, rep(0, lx)); length(Qi)
#'   
#'   TM <- create_A(lambda = lambda, mu = mu, nu = nu, q = q, k = soc,
#'                        max_number_of_species = lx); dim(TM); max(is.na(TM)); max(is.infinite(TM))
#'   
#'   # A2_v1 <- A_operator(Q = Qi, transition_matrix = TM, time_interval = total_time,
#'   # precision = 250L, methode = methode, A_abstol = abstol, A_reltol = reltol)
#'   
#'   A2_v1 <- try(expoRkit::expv(v = Qi, x = TM, t = total_time, m = 50L), silent = T)
#'   
#'   total_product <- A2_v1 * one_over_Cm * one_over_qm_binom
#'   Pc <- sum(total_product)
#'   
#'   return(list(Pc = Pc, A2_v1 = A2_v1))
#' }
#' 

