# @Giappo: add doc
#' Does something K
#' @inheritParams default_params_doc
#' @export
pmb_loglik <- function(pars, brts, soc = 2){

#BOTH LAmbdA AND NU
#setup
test_pars=pars  # test_pars=c(0.4,0,0.1,0.2)
test_brts=brts  # test_brts=c(-10,-5,-5,-3,-3,-3,-2,-1,-1,-1)
#numerical loglik
# loglik=mbd_loglik(pars = test_pars, brts = test_brts, soc = soc, cond = 0, tips_interval = c(0,Inf))

#theoretical loglik
N0 <- soc
lambda <- test_pars[1]; mu <- test_pars[2]; nu <- test_pars[3]; q <- test_pars[4];
condition1 <- ( any(is.nan(test_pars)) != 0 | any(is.infinite(test_pars)) != 0 )
condition2 <- ( lambda < 0 | mu != 0 | nu < 0 | q <= 0 | q >= 1)
if (condition1 | condition2){th_loglik = -Inf}else
{
  data <- brts2time_intervals_and_births(test_brts)
  time_intervals <- data$time_intervals
  births <- data$births
  k <- N0 + cumsum(c(0, births))
  A_term <- rep(1, length(time_intervals))     #branches
  B_term <- rep(1, length(time_intervals) - 1) #nodes
  #calculating branches contribution
  i <- 0:1e6
  for (t in 1:length(time_intervals))
  {
    poisson_term <- stats::dpois(i, nu*time_intervals[t], log = FALSE)[
      stats::dpois(i, nu*time_intervals[t], log = FALSE) != 0] #(nu*(t_k-t_k-1))^i*exp(-nu*(t_k-t_k-1))/k!
    ii <- i[stats::dpois(i, nu*time_intervals[t], log = FALSE) !=0 ]
    A_term[t] <- sum( (1 - q)^(ii * k[t]) * poisson_term ) *    # nu contribution: (1-q)^(k*i) * (nu*(t_k-t_k-1))^i*exp(-nu*(t_k-t_k-1))/k!
                 exp(- k[t] * lambda * (time_intervals[t])) # lambda contribution: exp(-k*lambda*(t_k-t_k-1))
  }
  #calculating nodes contribution
  B_term <- (nu * choose(k[-length(k)], births) * q^births * (1 - q)^(k[-length(k)] - births)) + # nu contribution: nu*(k,b)*q^b*(1-q)^(k-b)
            lambda * k[-length(k)] * (births == 1)                                               # lambda contribution: lambda*k (only if b==1)

  th_loglik <- sum(log(A_term)) + sum(log(B_term))
}
return(th_loglik)
}

# @Giappo: add doc
#' Does something L
#' @inheritParams default_params_doc
#' @export
pmb_loglik_Qvector <- function(pars, brts, soc = 2){
  #I would like to get the entire Q vector out of this...

  #BOTH LAmbdA AND NU
  #setup
  test_pars=pars  # test_pars=c(0.4,0,0.1,0.2)
  test_brts=brts  # test_brts=c(-10,-5,-5,-3,-3,-3,-2,-1,-1,-1)
  #numerical loglik
  # loglik=mbd_loglik(pars = test_pars, brts = test_brts, soc = soc, cond = 0, tips_interval = c(0,Inf))

  #theoretical loglik
  N0 <- soc
  lambda <- test_pars[1]; mu <- test_pars[2]; nu <- test_pars[3]; q <- test_pars[4];
  condition1 <- ( any(is.nan(test_pars)) != 0 | any(is.infinite(test_pars)) != 0 )
  condition2 <- ( lambda < 0 | mu != 0 | nu < 0 | q <= 0 | q >= 1)
  if (condition1 | condition2){th_loglik = -Inf}else
  {
    data <- brts2time_intervals_and_births(test_brts)
    time_intervals <- data$time_intervals
    births <- data$births
    k <- N0 + cumsum(c(0, births))
    A_term <- rep(1, length(time_intervals)    ) #branches
    B_term <- rep(1, length(time_intervals) - 1) #nodes
    #calculating branches contribution
    i <- 0:1e6
    for (t in 1:length(time_intervals))
    {
      poisson_term <- stats::dpois(
          i, nu * time_intervals[t], log = FALSE
        )[stats::dpois(i, nu * time_intervals[t], log = FALSE) != 0] #(nu*(t_k-t_k-1))^i*exp(-nu*(t_k-t_k-1))/i!
      ii <- i[stats::dpois(i, nu * time_intervals[t], log = FALSE) != 0]
      A_term[t] <- sum( (1 - q)^(ii * k[t]) * poisson_term ) *  # nu contribution: (1-q)^(k*i) * (nu*(t_k-t_k-1))^i*exp(-nu*(t_k-t_k-1))/i!
                   exp(- k[t] * lambda * (time_intervals[t]))   # lambda contribution: exp(-k*lambda*(t_k-t_k-1))
    }
    #calculating nodes contribution
    B_term <- (nu * choose(k[-length(k)], births) * q^births * (1 - q)^(k[-length(k)] - births)) + # nu contribution: nu*(k,b)*q^b*(1-q)^(k-b)
              lambda * k[-length(k)] * (births == 1)                                               # lambda contribution: lambda*k (only if b==1)

    th_loglik <- sum(log(A_term)) + sum(log(B_term))
  }
  return(th_loglik)
}

# @Giappo: add doc
#' Does something M
#' @inheritParams default_params_doc
#' @export
pmb_loglik_choosepar <- function(trparsopt, trparsfix = 0, idparsopt = c(1,3,4),
                                 idparsfix = (1:4)[-idparsopt], brts, soc = 2,
                                 pars.transform = 0){

  #This function provides a likelihood for a subset of parameters. This is built to work inside mbd_minusLL_vs_single_parameter or any optimizer like optim or simplex
  #idparsopt are the ids of the parameters you want to analyze
  #trparsopt are the values for parameters you want to analyze
  #idparsfix are the ids of the parameters you want to fix
  #trparsfix are the values for parameters you want to fix

  cond <- 0
  namepars <- c("lambda","mu","nu","q"); Npars <- length(namepars);
  if (length(trparsopt) == 4 && missing(trparsfix)){trparsfix <- NULL}
  # idparsopt=(1:3)[-c(idparsfix)] #this argument is useless but I let the user specify it because Rampal also did it (for some reason)
  trpars1 <- rep(0, Npars)
  trpars1[idparsopt] <- trparsopt
  if(length(idparsfix) != 0)
  {
    trpars1[idparsfix] <- trparsfix
  }
  if( min(trpars1[1:Npars]) < 0 ){loglik <- -Inf}else
  {
    if (pars.transform == 1)
    {
      #Rampal's transformation
      pars1 = trpars1/(1 - trpars1)
    }else
    {
      pars1 <- trpars1
    }
    loglik <- pmb_loglik(pars = pars1, brts = brts)
  }
  if(is.nan(loglik) || is.na(loglik))
  {
    cat("There are parameter values used which cause numerical problems:",trpars1,"\n")
    loglik <- -Inf
  }
  return(loglik)
}

# @Giappo: add doc
#' Does something N
#' @inheritParams default_params_doc
#' @export
pmb_ML <- function(brts, initparsopt, soc = 2,
                 res = 10 * (1+length(brts)+missnumspec), tol = c(1E-3, 1E-4, 1E-6),
                 maxiter = 1000 * round((1.25)^length(idparsopt)),
                 changeloglikifnoconv = FALSE, optimmethod = 'simplex',
                 pars.transform = 1)
{# bracket#1
  # - tol = tolerance in optimization
  # - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
  # - maxiter = the maximum number of iterations in the optimization
  # - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
  # - optimmethod = 'simplex' (current default) or 'subplex' (default of previous versions)
  idparsfix <- 2
  parsfix   <- 0
  idparsopt <- c(1,3,4)
  # if (missing(parsfix)&&(length(idparsfix)==0)){parsfix=NULL}

  options(warn=-1)
  namepars <- c("lambda","mu","nu","q"); Npars <- length(namepars); #if you add more parameters to your model just change this
  failpars <- rep(-1, Npars); names(failpars) <- namepars; #those are the parameters that you get if something goes sideways
  if(is.numeric(brts) == FALSE)
  {# bracket#2
    cat("The branching times should be numeric.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
  } else {
    idpars <- sort(c(idparsopt, idparsfix))
    if ( (sum(idpars == (1:Npars)) != Npars) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)) )
    {# bracket#3
      cat("The parameters to be optimized and/or fixed are incoherent.\n")
      out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    } else {
      if (length(namepars[idparsopt]) == 0) { optstr <- "nothing" } else { optstr <- namepars[idparsopt] }
      cat("You are optimizing",optstr,"\n")
      if (length(namepars[idparsfix]) == 0) { fixstr <- "nothing" } else { fixstr <- namepars[idparsfix] }
      cat("You are fixing",fixstr,"\n")
      cat("Optimizing the likelihood - this may take a while.","\n")
      utils::flush.console()
      if (pars.transform == 1)
      {
        trparsopt = initparsopt/(1 + initparsopt)
        trparsopt[which(initparsopt == Inf)] = 1
        trparsfix = parsfix/(1 + parsfix)
        trparsfix[which(parsfix == Inf)] = 1
      }else
      {
        trparsopt <- initparsopt
        trparsfix <- parsfix
      }
      optimpars = c(tol, maxiter)
      initloglik <- pmb_loglik_choosepar(trparsopt = trparsopt, trparsfix = trparsfix, idparsopt = idparsopt,
                                               idparsfix = idparsfix, brts = brts, soc = soc,
                                               pars.transform = pars.transform) #there's no pars2 here and instead 3 more args at the end
      cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
      utils::flush.console()
      if(initloglik == -Inf)
      {# bracket#4
        cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
        out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
      } else {
        out <- DDD::optimizer(optimmethod = optimmethod, optimpars = optimpars,
                               fun = pmb_loglik_choosepar, trparsopt = trparsopt,
                               trparsfix = trparsfix, idparsopt = idparsopt,
                               idparsfix = idparsfix, brts = brts, soc = soc,
                               pars.transform = pars.transform)
        if(out$conv != 0)
        {# bracket#5
          cat("Optimization has not converged. Try again with different initial values.\n")
          out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
        } else {
          MLtrpars <- as.numeric(unlist(out$par))
          if (pars.transform == 1)
          {
            MLpars = MLtrpars/(1-MLtrpars)
          }else
          {
            MLpars <- MLtrpars
          }
          MLpars1 <- rep(0,Npars); names(MLpars1)=namepars
          MLpars1[idparsopt] <- MLpars
          if (length(idparsfix) != 0) {MLpars1[idparsfix] <- parsfix}
          ML <- as.numeric(unlist(out$fvalues))
          out2 <- data.frame(t(MLpars1), loglik = ML, df = length(initparsopt), conv = unlist(out$conv))

          tobeprint <- "Maximum likelihood parameter estimates:"
          for (ii in 1:Npars)
          {
            tobeprint <- paste(tobeprint,paste(names(MLpars1[ii]),":",sep = ""),MLpars1[ii])
          }
          s1 <- sprintf(tobeprint)

          if(out2$conv != 0 & changeloglikifnoconv == T) {out2$loglik <- -Inf}
          s2 <- sprintf('Maximum loglikelihood: %f',ML)
          cat("\n",s1,"\n",s2,"\n\n")
        }# bracket#5
      }# bracket#4
    }# bracket#3
  }# bracket#2
  invisible(out2)
}# bracket#1

