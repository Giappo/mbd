# MAIN LOGLIK FUNCTION---------------------------------
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
#' @param tips_interval It takes into account tips boundaries constrain on simulated dataset.
#' @param missnumspec The number of species that are in the clade but missing in the phylogeny.
#' @param methode Specifies how the integration must be performed: set "sexpm" if you want to use sexpm; set "expo" if you want to use expoRkit; set "lsoda" if you want to use the "lsoda" method with the "deSolve::ode" function.
#' @param safety_threshold It determines the precision on the parameters.
#' @return The function returns the natural logarithm of the likelihood for the process.
#'
#' @examples
#' set.seed(11)
#' simulated_data = MBD:::mbd_sim( pars=c(0.6,0.1,2.2,0.1),soc=2,age=10,cond=1 )
#' plot(simulated_data$tas)
#' mbd_loglik( pars=c(0.8,0.05,2.2,0.1),brts=simulated_data$brts,soc=2,cond=1,missnumspec=0 )
#'
#' @export

mbd_loglik <- function(pars, brts, soc = 2, cond = 1, tips_interval = c(0, Inf),
                       missnumspec = 0, safety_threshold = 1e-3,
                       methode = "expo", alpha = 10, minimum_multiple_births = 0, print_errors = 1){
  
  #Optional stuff that I might need to run the program one line at the time:
  #brts=sim_data[[1]];missnumspec=0;pars=sim_pars;missing_interval=c(1,Inf);methode="expo"
  
  #BASIC SETTINGS AND CHECKS
  lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4];
  min_tips <- tips_interval[1]; max_tips <- tips_interval[2]; abstol <- 1e-16; reltol <- 1e-10;
  
  condition1 <- (any(is.nan(pars)) != 0 | any(is.infinite(pars)) != 0)
  condition2 <- (lambda < 0 | mu < 0 | nu < 0 |
                   q <= 0 + safety_threshold | q >= 1 - safety_threshold |
                   minimum_multiple_births < 0)
  condition3 <- (length(pars) != 4)
  if       (condition1)
  {
    if (print_errors == 1){print("input parameters are either infinite or NaN")}
    loglik <- -Inf
  }else if (condition2)
  {
    if (print_errors == 1){print("input parameters have wrong values")}
    loglik <- -Inf
  }else if (condition3)
  {
    if (print_errors == 1){print("wrong number of input parameters")}
    loglik <- -Inf
  }else if (mu == 0 && cond == 0 && tips_interval == c(0, Inf) &&
            missnumspec == 0 && minimum_multiple_births == 0)
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
    lx <- max_number_of_species
    
    #ALPHA ANALYSIS
    deltaAlpha <- 1; count <- 0; same_result_count <- 0; Pc.notanumber <- 1;
    while (Pc.notanumber)
    {
      Pc1 <- MBD:::calculate_conditional_probability(brts = brts, pars = pars, tips_interval = tips_interval, cond = cond, soc = soc, alpha = alpha, methode = methode,
                                                     abstol = abstol, reltol = reltol, minimum_multiple_births = minimum_multiple_births)$Pc
      Pc.notanumber <- is.nan(Pc1)
      alpha <- alpha - Pc.notanumber
    }
    while (deltaAlpha != 0 && count < 100 && same_result_count < 5)
    {
      Pc2 <- MBD:::calculate_conditional_probability(brts = brts, pars = pars, cond = cond,
                                                     soc = soc, tips_interval = tips_interval,
                                                     alpha = alpha + deltaAlpha, methode = methode,
                                                     abstol = abstol, reltol = reltol,
                                                     minimum_multiple_births = minimum_multiple_births)$Pc
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
      Pc1 <- MBD:::calculate_conditional_probability(brts = brts, pars = pars, tips_interval = tips_interval, cond = cond, soc = soc, alpha = alpha, methode = methode,
                                                  abstol = abstol, reltol = reltol, minimum_multiple_births = minimum_multiple_births)$Pc
    }
    Pc <- Pc1
    if (count >= 100){alpha <- 10}
    if (Pc <= 0 | Pc == Inf | Pc == -Inf){Pc <- 1; print("there's a problem with Pc")}

    start_over_again <- 1; iterations <- 0;
    negative_values <- nan_values <- 0;
    while (start_over_again == 1 & iterations < 100)
    {
      #MATRIX DIMENSION SETUP
      max_number_of_species <- alpha * max_k; #alpha is the proportionality factor between max_k and the edge of the matrix
      nvec <- 0:max_number_of_species

      #SETTING INITIAL CONDITIONS (there's always a +1 because of Q0)
      Qi <- c(1,rep(0,max_number_of_species))
      Qt <- matrix(0,ncol = (max_number_of_species+1), nrow = length(time_intervals))
      Qt[1,] <- Qi
      dimnames(Qt)[[2]] <- paste0("Q", 0:max_number_of_species)
      k <- N0 #N0 is the number of species at t=1
      t <- 2  #t is starting from 2 so everything is ok with birth[t] and time_intervals[t] vectors
      C <- rep(1, (length(time_intervals))); D <- C
      logB <- 0;

      #EVOLVING THE INITIAL STATE TO THE LAST BRANCHING POINT
      while (t <= length(time_intervals))
      {
        #Applying A operator
        transition_matrix <- MBD:::create_A(lambda = lambda, mu = mu, nu = nu, q = q, k = k,max_number_of_species = max_number_of_species)
        Qt[t,] <- MBD:::A_operator(Q = Qt[(t-1),], transition_matrix = transition_matrix, time_interval = time_intervals[t], precision = 50L, methode = methode, A_abstol = abstol, A_reltol = reltol)
        if (methode != "sexpm"){Qt[t,] <- MBD:::negatives_correction(Qt[t,], pars)} #it removes some small negative values that can occurr as bugs from the integration process
        if ( any(is.nan(Qt[t,])) )
        {
          if (Sys.info()[['sysname']] == "Windows")
          {
            print(pars); print(Qt[t,])
          }
          nan_values <- 1; break
        }
        if (any(Qt[t,]<0)){negative_values <- 1; break}

        #Applying C operator (this is a trick to avoid precision issues)
        C[t] <- 1/(sum(Qt[t,])); Qt[t,] <- Qt[t,] * C[t]

        if (t < length(time_intervals))
        {
          #Applying B operator
          B <- MBD:::create_B(lambda = lambda, nu = nu, q = q, k = k, b = births[t],
                              max_number_of_species = max_number_of_species)
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
      if (nan_values == 1)     {alpha <- alpha - 5; if(alpha <= 0){alpha <- 5}}
      iterations <- iterations + 1
      start_over_again <- (nan_values) || (negative_values)
    }

    #Selecting the state I am interested in
    vm <- 1/choose((k+missnumspec), k)
    P  <- vm * Qt[t, (missnumspec + 1)] #I have to include +1 because of Q0

    #Removing C and D effects from the LL
    loglik <- log(P) + logB - sum(log(C)) - sum(log(D))

    #Various checks
    loglik <- as.numeric(loglik)
    if (is.nan(loglik) | is.na(loglik))
    {
      loglik <- -Inf
    }else
    {
      loglik <- loglik - log(Pc) #conditioned likelihood
    }

  }
  return(loglik)
}
# OLD CONDITIONING  ---------------------------------
# OLD WAY TO CONDITIONING THE LIKELIHOOD ON THE SURVIVAL OF CROWN SPECIES AND ADJUST ALPHA
# Pc=-1;deltaAlpha=1;A2_v1=1;iterations = 0
# while (( deltaAlpha!=0 | Pc<0 | any(A2_v1<0) ) & iterations < 100){
#   Pc1_list=MBD:::calculate_conditional_probability(brts=brts,pars=pars,tips_interval=tips_interval,cond=cond,soc=soc,alpha=alpha,methode=methode,
#                                                    abstol=abstol, reltol=reltol, minimum_multiple_births=minimum_multiple_births)
#   Pc1=Pc1_list$Pc
#   A2_v1=Pc1_list$A2_v1
#   if ( is.nan(Pc1) | any( is.nan(A2_v1) )){deltaAlpha = -1;Pc1=Pc2=-1}else
#   {
#     Pc2=0
#     Pc2=MBD:::calculate_conditional_probability(brts=brts,pars=pars,tips_interval=tips_interval,cond=cond,soc=soc,alpha=alpha+5,methode=methode,
#                                                 abstol=abstol, reltol=reltol, minimum_multiple_births=minimum_multiple_births)$Pc
#     if ( is.nan(Pc2) ){deltaAlpha = -21;Pc1=Pc2=-1}else{deltaPc=(Pc2-Pc1)/Pc2;deltaAlpha=floor( (deltaPc*10) )}
#   }
#
#   alpha=alpha+deltaAlpha; if (alpha<=0){alpha=5}
#   Pc=Pc1
#   iterations = iterations + 1
#   flush.console();print(alpha)
# }


# PROTOTYPE LOGLIK FUNCTION  ---------------------------------
#' @author Giovanni Laudanno
#' @title Calculates the likelihood for a multiple birth-death process
#' @description mbd_loglik0 provides the likelihood for a process in which multiple births (from different parents) at the same time are possible.
#' @param pars vector of parameters:
#' \itemize{
#'   \item pars[1] is the multiple speciation trigger rate;
#'   \item pars[2] is the extinction rate;
#'   \item pars[3] is the single-lineage speciation probability.
#' }
#' @param brts A set of branching times of a phylogeny.
#' @param soc Sets whether stem or crown age should be used (1 or 2)
#' @param cond Set 1 if you want to condition on stem or crown age and non-extinction of the phylogeny. Set 0 otherwise.
#' @param tips_interval It takes into account tips boundaries constrain on simulated dataset.
#' @param missnumspec The number of species that are in the clade but missing in the phylogeny.
#' @param methode Specifies how the integration must be performed: set "sexpm" if you want to use sexpm; set "expo" if you want to use expoRkit; set "lsoda" if you want to use the "lsoda" method with the "deSolve::ode" function.
#' @param safety_threshold It determines the precision on the parameters.
#' @return The function returns the natural logarithm of the likelihood for the process.
#'
#' @examples
#' set.seed(11)
#' simulated_data = MBD:::mbd_sim0( pars=c(2.5,0.1,0.1),soc=2,age=10,cond=1 )
#' plot(simulated_data$tas)
#' mbd_loglik0( pars=c(1.2,0.05,0.1),brts=simulated_data$brts,soc=2,cond=1,missnumspec=0 )
#'
#' @export

#MAIN
mbd_loglik0 <- function(pars, brts, soc = 2, cond = 0, tips_interval = c(0,Inf), missnumspec = 0,
                       safety_threshold = 1e-4,methode = "expo", debug_check = 0,alpha = 20){

  #Optional stuff that I might need to run the program one line at the time:
  #brts=sim_data[[1]];missnumspec=0;pars=sim_pars;missing_interval=c(1,Inf)

  #BASIC SETTINGS AND CHECKS
  lambda=pars[1]; mu=pars[2]; q=pars[3]; min_tips=tips_interval[1]; max_tips=tips_interval[2]; abstol=1e-16; reltol=1e-10;
  starting_alpha = alpha

  condition1 = ( any(is.nan(pars))!=0 | any(is.infinite(pars))!=0 )
  condition2 = ( mu<0 | lambda<=0+safety_threshold | q<=0+safety_threshold | q>=1-safety_threshold ) # mu>=lambda |
  if (condition1 | condition2){loglik = -Inf}
  else if(length(pars)!=3){print("input parameters are wrong");loglik = -Inf}
  else{
    Pc=-1
    while (Pc<0 && alpha <= max(starting_alpha,80)){

  #ADJUSTING DATA
  data=MBD:::brts2time_intervals_and_births(brts)
  time_intervals=c(0,data$time_intervals)
  births=c(0,data$births)

  #SET UP
  N0 = soc #number of starting species
  k_interval = N0 + cumsum(births)
  max_k = max(k_interval)
  max_number_of_species = alpha*max_k; #alpha is the proportionality factor between max_k and the edge of the matrix
  nvec=0:max_number_of_species

  #SETTING INITIAL CONDITIONS (there's always a +1 because of Q0)
  Qi=c(1,rep(0,max_number_of_species))
  Qt=matrix(0,ncol = (max_number_of_species+1), nrow = length(time_intervals))
  Qt[1,]=Qi
  dimnames(Qt)[[2]] = paste("Q",0:max_number_of_species,sep="")
  k=N0 #N0 is the number of species at t=1
  t=2  #t is starting from 2 so everything is ok with birth[t] and time_intervals[t] vectors
  C=rep(1,(length(time_intervals))); D=C
  logB=0;

  #EVOLVING THE INITIAL STATE TO THE LAST BRANCHING POINT
  while ( t<length(time_intervals) ){

    #Applying A operator
    transition_matrix=create_A0(max_number_of_species = max_number_of_species,lambda = lambda,mu = mu,q = q,k = k)
    Qt[t,]=MBD:::A_operator(Q = Qt[(t-1),],transition_matrix = transition_matrix,time_interval = time_intervals[t],precision = 50L,methode=methode,A_abstol=abstol,A_reltol=reltol)
    if (methode!="sexpm"){Qt[t,]=MBD:::negatives_correction(Qt[t,],pars)} #it removes some small negative values that can occurr as bugs from the integration process

    #Applying C operator (this is a trick to avoid precision issues)
    if (debug_check==1){print(head(Qt[t,]))}
    C[t]=1/(sum(Qt[t,]))
    Qt[t,]=Qt[t,]*C[t]

    #Applying B operator
    B <- create_B0(max_number_of_species = max_number_of_species, q = q, k = k, b = births[t])
    # B[(row(B)>(2*col(B)+k-births[t])) | col(B)>row(B) ]=0 #this is a constrain due to maximum number of speciations being (2*n+k-b); probably it is redundant
    # if (max(is.nan(B))>0){print(paste("NaN were produced in the B matrix at time=",t))}
    Qt[t,]=(B %*% Qt[t,])
    if (methode!="sexpm"){Qt[t,]=MBD:::negatives_correction(Qt[t,],pars)}
    logB = logB + log(lambda) + lchoose(k,births[t]) + births[t]*log(q)

    #Applying D operator (this works exactly like C)
    if (debug_check==1){print(head(Qt[t,]))}
    D[t]=1/(sum(Qt[t,]))
    Qt[t,]=Qt[t,]*D[t]

    #Updating running parameters
    k=k+births[t]
    t=t+1
  }

  #Applying A operator from the last branching time to the present
  transition_matrix=create_A0(max_number_of_species = max_number_of_species,lambda = lambda,mu = mu,q = q,k = k)
  Qt[t,]=MBD:::A_operator(Q = Qt[(t-1),],transition_matrix = transition_matrix,time_interval = time_intervals[t],precision = 50L,methode=methode,A_abstol=abstol,A_reltol=reltol)
  if (methode!="sexpm"){Qt[t,]=MBD:::negatives_correction(Qt[t,],pars)}
  if (debug_check==1){print(head(Qt[t,]))}

  #Selecting the state I am interested in
  vm = 1/choose((k+missnumspec),k)
  P  = vm * Qt[t,(missnumspec+1)] #I have to include +1 because of Q0

  #Removing C and D effects from the LL
  loglik = log(P) + logB - sum(log(C)) - sum(log(D))

  #Various checks
  loglik = as.numeric(loglik)
  if ( is.nan(loglik) | is.na(loglik) ){loglik=-Inf}

  #CONDITIONING THE LIKELIHOOD ON THE SURVIVAL OF CROWN SPECIES
  Pc=1
  if ( (cond==1 | tips_interval[1]>1 | tips_interval[2]<Inf ) & !is.infinite(loglik) ){ #difference between sexpm and expo are not here

    total_time=max(abs(brts));
    m=0:max_number_of_species
    one_over_Cm=(3*(m+1))/(m+3)
    one_over_qm_binom=1/choose((m+N0),N0)
    tips_components=(1+min_tips):(1+min(max_tips,max_number_of_species)) #applying tips constrain

    Mk_N0=MBD:::create_A0(max_number_of_species = max_number_of_species,lambda = lambda,mu = mu,q = q,k = N0)
    A2_v1=MBD:::A_operator(Q = Qt[1,],transition_matrix = Mk_N0,time_interval = total_time,precision = 50L,methode=methode,A_abstol=abstol,A_reltol=reltol)
    if (methode!="sexpm"){A2_v1=MBD:::negatives_correction(A2_v1,pars)} #it removes some small negative values that can occurr as bugs from the integration process
    if (debug_check==1){print(head(A2_v1,max_tips))}
    total_product=A2_v1*one_over_Cm*one_over_qm_binom
    Pc=sum(total_product[tips_components])

    if (Pc==0){#slowest and best accuracy
      # ode_matrix=as.matrix(Mk_N0) #use this only if you use sparsematrices
      ode_matrix=MK_N0
      times=c(0,total_time)
      A2_v1=deSolve::ode(y = Qt[1,], times = times, func = MBD:::mbd_loglik_rhs, parms = ode_matrix,atol=abstol,rtol=reltol)[2,-1] #evolving crown species to the present
      total_product=A2_v1*one_over_Cm*one_over_qm_binom
      Pc=sum(total_product[tips_components])
    }

  }
  alpha = alpha + 5
  }
  loglik=loglik-log(Pc) #conditioned likelihood

  }
  # loglik=-loglik #Rampal's optimizer uses loglik rather than -loglik
  return(loglik)
}
