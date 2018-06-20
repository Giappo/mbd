#' @title Internal MBD function
#' @description Internal MBD function.
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
# if (.Platform$OS.type=="unix")   {matrix_builder = MBD:::hyperA_HannoX} #this is set to work with home directory of the cluster

#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
create_A0 <- function(max_number_of_species, lambda, mu, q, k, matrix_builder = hyperA_HannoX){
  nvec = 0:max_number_of_species
  M = lambda * matrix_builder(N = max_number_of_species,k = k,q = q)

  # diag(M) = (-lambda)*c( (1-(1-q)^(k+( nvec[1:max_number_of_species] ) )),0 ) - mu*(nvec+k)
  diag(M) = (-lambda) * ( 1-(1-q)^(k + nvec) ) - mu * (nvec + k) #new version to avoid the dumpster problem at the end of the matrix

  M[row(M) == col(M) - 1] = mu*nvec[2:(max_number_of_species+1)]
  return(M)
}

#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
create_B0 <- function(max_number_of_species, q, k, b, matrix_builder = hyperA_HannoX){#lambda * choose(k,b) * q^b  is going to be added in logB in the main script
  k2 <- k - b
  B <- matrix_builder(N = max_number_of_species, k = k2, q = q)
  return(B)
}

#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
create_A <- function(lambda, mu, nu, q, k, max_number_of_species){
  nvec <- 0:max_number_of_species
  M <- MBD:::create_A0(max_number_of_species = max_number_of_species,
                       lambda = nu, mu = mu, q = q, k = k)
  M[row(M) == col(M) + 1] <- M[row(M) == col(M) + 1] + lambda * (nvec[1:(max_number_of_species)]+2*k)

  # M[row(M) == col(M)] = M[row(M) == col(M)] - c(lambda*(nvec[-length(nvec)]+k),0)
  M[row(M) == col(M)] <- M[row(M) == col(M)] - lambda * (nvec + k) #new version to avoid the dumpster problem at the end of the matrix

  return(M)
}

#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
create_B <- function(lambda, nu, q, k, b, max_number_of_species){
  M <- MBD:::create_B0(max_number_of_species = max_number_of_species, q = q, k = k, b = b)
  B <- lambda * k * diag(max_number_of_species + 1) * (b == 1) + nu * choose(k, b) * (q^b) * M
  return(B)
}

#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
create_A.no_mbd = function(lambda,mu,nu,q,k,max_number_of_species,minimum_multiple_births){
  nvec <- 0:max_number_of_species
  M <- MBD:::create_A0(max_number_of_species = max_number_of_species,lambda = nu,mu = mu,q = q,k = k)
  M[row(M) == col(M) + 1] <- M[row(M) == col(M) + 1] + lambda * (nvec[1:(max_number_of_species)]+2*k)
  M[row(M) == col(M)]     <- M[row(M) == col(M)] - c(lambda * (nvec[-length(nvec)] + k), 0)
  M[row(M) > col(M) + minimum_multiple_births] <- 0
  return(M)
}

#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
create_B.no_mbd = function(lambda,nu,q,k,b,max_number_of_species,minimum_multiple_births){
  M <- MBD:::create_B0(max_number_of_species = max_number_of_species, q = q, k = k, b = b)
  B <- lambda * k * diag(max_number_of_species + 1) * (b == 1) + nu * choose(k,b) * (q^b) * M
  B[row(B) > col(B) + minimum_multiple_births] = 0
  return(B)
}

#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
calculate_conditional_probability <- function (brts, pars, tips_interval = c(0, Inf),
                                               cond = 1, soc = 2, alpha, methode = "expo",
                                               abstol = 1e-16, reltol = 1e-10,
                                               minimum_multiple_births = 0){

  lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4];
  min_tips <- tips_interval[1]; max_tips <- tips_interval[2];
  min_tips <- max(min_tips, soc * cond) #check this
  N0 <- soc
  total_time <- max(abs(brts));
  births <- c(0, MBD:::brts2time_intervals_and_births(brts)$births)
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
    Mk_N0 <- MBD:::create_A(lambda = lambda, mu = mu, nu = nu, q = q, k = soc,
                            max_number_of_species = max_number_of_species)
    A2_v1 <- MBD:::A_operator(Q = Qi, transition_matrix = Mk_N0, time_interval = total_time,
                              precision = 50L, methode = methode, A_abstol = abstol, A_reltol = reltol)
    if (methode != "sexpm"){A2_v1 <- MBD:::negatives_correction(A2_v1, pars)} #it removes some small negative values that can occurr as bugs from the integration process

    if (minimum_multiple_births > 0) #adjust for the required minimum amount of mbd
    {
      Mk_N0.no_mbd <- MBD:::create_A.no_mbd(lambda = lambda, mu = mu, nu = nu, q = q, k = soc, max_number_of_species = max_number_of_species, minimum_multiple_births = minimum_multiple_births)
      A2_v1.no_mbd <- MBD:::A_operator(Q = Qi, transition_matrix = Mk_N0.no_mbd, time_interval = total_time,precision = 50L,methode=methode,A_abstol=abstol,A_reltol=reltol)
      A2_v1 <- A2_v1 - A2_v1.no_mbd
    }

    total_product <- A2_v1 * one_over_Cm * one_over_qm_binom
    Pc <- sum(total_product[tips_components])
  }
  return(list(Pc = Pc, A2_v1 = A2_v1))
}

#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
calculate_conditional_probability2 <- function (brts, pars, missing_tips_interval = c(0, Inf),
                                                soc = 2, alpha, methode = "expo",
                                                abstol = 1e-16, reltol = 1e-10,
                                                minimum_multiple_births = 0){

  lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4];
  min_tips <- missing_tips_interval[1]; max_tips <- missing_tips_interval[2];
  N0 <- soc
  total_time <- max(abs(brts))
  births <- c(0, MBD:::brts2time_intervals_and_births(brts)$births)
  k_interval <- N0 + cumsum(births)
  max_k <- max(k_interval)
  max_number_of_species <- alpha * max_k; #alpha is the proportionality factor between max_k and the edge of the matrix

  m <- 0:max_number_of_species
  one_over_Cm <- (3 * (m + 1))/(m + 3)
  one_over_qm_binom <- 1/choose((m + N0), N0)
  tips_components <- (1 + min_tips):(1 + min(max_tips, max_number_of_species)) #applying tips constrain

  Qi    <- c(1, rep(0, max_number_of_species))
  Mk_N0 <- MBD:::create_A(lambda = lambda, mu = mu, nu = nu, q = q, k = soc,
                          max_number_of_species = max_number_of_species)
  A2_v1 <- MBD:::A_operator(Q = Qi, transition_matrix = Mk_N0, time_interval = total_time,
                            precision = 50L, methode = methode, A_abstol = abstol, A_reltol = reltol)

  # MBD::mbd_loglik(pars = pars, brts = c(total_time), soc = soc, cond = 0, missnumspec = 0)

  total_product <- A2_v1 * one_over_Cm * one_over_qm_binom
  Pc <- sum(total_product[tips_components]); Pc





  if (min_tips < soc)
  min_tips <- max(min_tips, soc * cond) #check this
  N0 <- soc
  total_time <- max(abs(brts));
  births <- c(0, MBD:::brts2time_intervals_and_births(brts)$births)
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
    Mk_N0 <- MBD:::create_A(lambda = lambda, mu = mu, nu = nu, q = q, k = soc,
                            max_number_of_species = max_number_of_species)
    A2_v1 <- MBD:::A_operator(Q = Qi, transition_matrix = Mk_N0, time_interval = total_time,
                              precision = 50L, methode = methode, A_abstol = abstol, A_reltol = reltol)
    if (methode != "sexpm"){A2_v1 <- MBD:::negatives_correction(A2_v1, pars)} #it removes some small negative values that can occurr as bugs from the integration process

    if (minimum_multiple_births > 0) #adjust for the required minimum amount of mbd
    {
      Mk_N0.no_mbd <- MBD:::create_A.no_mbd(lambda = lambda, mu = mu, nu = nu, q = q, k = soc, max_number_of_species = max_number_of_species, minimum_multiple_births = minimum_multiple_births)
      A2_v1.no_mbd <- MBD:::A_operator(Q = Qi, transition_matrix = Mk_N0.no_mbd, time_interval = total_time,precision = 50L,methode=methode,A_abstol=abstol,A_reltol=reltol)
      A2_v1 <- A2_v1 - A2_v1.no_mbd
    }

    total_product <- A2_v1 * one_over_Cm * one_over_qm_binom
    Pc <- sum(total_product[tips_components])
  }
  return(list(Pc = Pc, A2_v1 = A2_v1))
}

#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
A_operator <- function(Q, transition_matrix, time_interval, precision = 50L,
                       A_abstol = 1e-16, A_reltol = 1e-10, methode = "expo"){

  precision_limit <- 3000

  if (methode == "sexpm")
  {
    exp_matrix <- rsexpm:::sexpm(transition_matrix * time_interval)
    result     <- exp_matrix %*% Q
  }else if (methode == "expo")
  {
    result <- try(expoRkit:::expv(v = Q, x = transition_matrix, t = time_interval, m = precision), silent = T)
    while ( ( any(!is.numeric(result)) || any(is.nan(result)) ) && precision < precision_limit )
    {
      precision <- precision + 200
      result <- try(expoRkit:::expv(v = Q, x = transition_matrix, t = time_interval, m = precision), silent = T)}
    if (!any(is.nan(result)))
    {
      while (any(result < 0) && precision < precision_limit)
      {
        precision <- precision + 200
        result <- try(expoRkit:::expv(v = Q, x = transition_matrix, t = time_interval, m = precision), silent = T)
      }
    }
  }else if (methode=="lsoda")
  {
    times <- c(0, time_interval)
    ode_matrix <- transition_matrix
    # result<-deSolve::ode(y = Q, times = times, func = MBD:::mbd_loglik_rhs, parms = ode_matrix,atol=A_abstol,rtol=A_reltol)[2,-1]
    R.utils:::evalWithTimeout(result <- deSolve::ode(y = Q, times = times, func = MBD:::mbd_loglik_rhs, parms = ode_matrix,atol=A_abstol,rtol=A_reltol)[2,-1], timeout = 1000)
  }

  if (any(!is.numeric(result)) || any(is.nan(result))) #sometimes expoRkit gives weird negative values. In this case perform standard lsoda integration.
  {
    times <- c(0, time_interval)
    ode_matrix <- transition_matrix
    # result<-deSolve::ode(y = Q, times = times, func = MBD:::mbd_loglik_rhs, parms = ode_matrix,atol=A_abstol,rtol=A_reltol)[2,-1]
    R.utils:::evalWithTimeout(result <- deSolve::ode(y = Q, times = times, func = MBD:::mbd_loglik_rhs, parms = ode_matrix,atol=A_abstol,rtol=A_reltol)[2,-1], timeout = 1000)
  }
  else if (any(result < 0)) #sometimes expoRkit gives weird negative values. In this case perform standard lsoda integration.
  {
    times <- c(0, time_interval)
    ode_matrix <- transition_matrix
    # result=deSolve::ode(y = Q, times = times, func = MBD:::mbd_loglik_rhs, parms = ode_matrix,atol=A_abstol,rtol=A_reltol)[2,-1]
    R.utils:::evalWithTimeout(result <- deSolve::ode(y = Q, times = times, func = MBD:::mbd_loglik_rhs, parms = ode_matrix,atol=A_abstol,rtol=A_reltol)[2,-1], timeout = 1000)
  }

  # if ( ( any(!is.numeric(result)) || any(is.nan(result)) ) && methode!="sexpm"){
  #     # to run this you actually need sexpm to be installed which, on the cluster, might not be the case.
  #     exp_matrix = rsexpm:::sexpm(transition_matrix*time_interval)
  #     result = exp_matrix %*% Q
  # }

  return(result)
}

#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
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
