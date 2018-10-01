# @Giappo: add doc
#' Does something K
#' @inheritParams default_params_doc
#' @export
pmb_loglik <- function(
  pars,
  brts,
  soc = 2
) {

  #BOTH LAmbdA AND NU
  #setup
  test_pars <- pars
  test_brts <- brts
  #numerical loglik

  #theoretical loglik
  init_n_lineages <- soc
  lambda <- test_pars[1]
  mu <- test_pars[2]
  nu <- test_pars[3]
  q <- test_pars[4];
  condition1 <- (any(is.nan(test_pars)) != 0 |
    any(is.infinite(test_pars)) != 0
  )
  condition2 <- (lambda < 0 | mu != 0 | nu < 0 | q <= 0 | q >= 1)
  if (condition1 | condition2) {
    th_loglik <- -Inf
  } else {
    data <- brts2time_intervals_and_births(test_brts)  # nolint internal function
    time_intervals <- data$time_intervals
    births <- data$births
    k <- init_n_lineages + cumsum(c(0, births))
    a_term <- rep(1, length(time_intervals))     #branches
    b_term <- rep(1, length(time_intervals) - 1) #nodes
    #calculating branches contribution
    i <- 0:1e6
    for (t in 1:length(time_intervals)) {
      #(nu *(t_k-t_k-1))^i * exp(-nu * (t_k - t_k - 1)) / k!
      poisson_term <- stats::dpois(i, nu * time_intervals[t], log = FALSE)[
        stats::dpois(i, nu * time_intervals[t], log = FALSE) != 0]
      ii <- i[stats::dpois(i, nu * time_intervals[t], log = FALSE) != 0]
      # (1) nu contribution: (1-q)^(k * i) * (nu *(t_k-t_k-1))^i
      #                      * exp(-nu *(t_k-t_k-1)) / k!
      # (2) lambda contribution: exp(-k * lambda *(t_k-t_k-1))
      a_term[t] <- sum( (1 - q) ^ (ii * k[t]) * poisson_term ) * # (1)
                   exp(-k[t] * lambda * (time_intervals[t]))   # (2)
    }
    #calculating nodes contribution
    # (1) nu contribution: nu *(k, b)* q^b *(1-q)^(k-b)
    # (2) lambda contribution: lambda * k (only if b==1)
    b_term <- (
      nu * choose(k[-length(k)], births) * q^births *  # (1)
      (1 - q) ^ (k[-length(k)] - births)               # (1)
    ) + lambda * k[-length(k)] * (births == 1)         # (2)

    th_loglik <- sum(log(a_term)) + sum(log(b_term))
  }
  th_loglik
}

# @Giappo: add doc
#' Does something L
#' @inheritParams default_params_doc
#' @export
pmb_loglik_q_vector <- function(pars, brts, soc = 2){
  #I would like to get the entire Q vector out of this...

  #BOTH LAmbdA AND NU
  #setup
  test_pars <- pars
  test_brts <- brts
  #numerical loglik

  #theoretical loglik
  init_n_lineages <- soc
  lambda <- test_pars[1]
  mu <- test_pars[2]
  nu <- test_pars[3]
  q <- test_pars[4];
  condition1 <- (any(is.nan(test_pars)) != 0 | any(is.infinite(test_pars)) != 0)
  condition2 <- (lambda < 0 | mu != 0 | nu < 0 | q <= 0 | q >= 1)
  if (condition1 | condition2) {
    th_loglik <- -Inf
  } else {
    data <- brts2time_intervals_and_births(test_brts) # nolint internal function
    time_intervals <- data$time_intervals
    births <- data$births
    k <- init_n_lineages + cumsum(c(0, births))
    a_term <- rep(1, length(time_intervals)    ) #branches
    b_term <- rep(1, length(time_intervals) - 1) #nodes
    #calculating branches contribution
    i <- 0:1e6
    for (t in 1:length(time_intervals)) {
      #(nu *(t_k-t_k-1))^i * exp(-nu *(t_k-t_k-1)) / i!
      poisson_term <- stats::dpois(
          i, nu * time_intervals[t], log = FALSE
        )[stats::dpois(i, nu * time_intervals[t], log = FALSE) != 0]
      ii <- i[stats::dpois(i, nu * time_intervals[t], log = FALSE) != 0]
      # (1) nu contribution: (1-q)^(k * i) * (nu *(t_k-t_k-1))^i * exp(-nu *(t_k-t_k-1)) / i!
      # (2) lambda contribution: exp(-k * lambda *(t_k-t_k-1))
      a_term[t] <- sum( (1 - q) ^ (ii * k[t]) * poisson_term ) * # (1)
                   exp(-k[t] * lambda * (time_intervals[t]))     # (2)
    }
    #calculating nodes contribution
    # (1) nu contribution: nu *(k, b) * q ^ b * (1 - q) ^ (k - b)
    # (2) lambda contribution: lambda * k (only if b==1)
    b_term <- (
      nu * choose(k[-length(k)], births) * q ^ births *  # (1)
        (1 - q)^(k[-length(k)] - births)                 # (1)
    ) + lambda * k[-length(k)] * (births == 1)           # (2)

    th_loglik <- sum(log(a_term)) + sum(log(b_term))
  }
  th_loglik
}

# @Giappo: add doc
#' Does something M
#' @inheritParams default_params_doc
#' @export
pmb_loglik_choosepar <- function(
  trparsopt,
  trparsfix = 0,
  idparsopt = c(1, 3, 4),
  idparsfix = (1:4)[-idparsopt],
  brts,
  soc = 2,
  pars_transform = 0
) {
  #This function provides a likelihood for a subset of parameters.
  # This is built to work inside mbd_minusLL_vs_single_parameter
  # or any optimizer like optim or simplex
  #idparsopt are the ids of the parameters you want to analyze
  #trparsopt are the values for parameters you want to analyze
  #idparsfix are the ids of the parameters you want to fix
  #trparsfix are the values for parameters you want to fix

  namepars <- mbd::get_mbd_param_names()
  n_pars <- length(namepars)
  if (length(trparsopt) == 4 && missing(trparsfix)) {
    trparsfix <- NULL
  }
  trpars1 <- rep(0, n_pars)
  trpars1[idparsopt] <- trparsopt
  if (length(idparsfix) != 0) {
    trpars1[idparsfix] <- trparsfix
  }
  if (min(trpars1[1:n_pars]) < 0) {
    loglik <- -Inf
  } else {
    if (pars_transform == 1) {
      #Rampal's transformation
      pars1 <- trpars1 / (1 - trpars1)
    } else {
      pars1 <- trpars1
    }
    loglik <- mbd::pmb_loglik(pars = pars1, brts = brts)
  }
  if (is.nan(loglik) || is.na(loglik)) {
    cat("There are parameter values used which cause numerical problems:", trpars1, "\n")
    loglik <- -Inf
  }
  return(loglik)
}

# @Giappo: add doc
#' Does something N
#' @inheritParams default_params_doc
#' @export
pmb_ml <- function(
  brts,
  initparsopt,
  idparsopt = c(1, 3, 4), # RJCB: No idea why @Giappo uses that as a default
  soc = 2,
  res = 10 * (1 + length(brts) + missnumspec),
  tol = c(1E-3, 1E-4, 1E-6),
  maxiter = 1000 * round((1.25) ^ length(idparsopt)),
  changeloglikifnoconv = FALSE,
  optimmethod = "simplex",
  pars_transform = 1,
  missnumspec = 0
) {
  # - tol = tolerance in optimization
  # - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
  # - maxiter = the maximum number of iterations in the optimization
  # - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
  # - optimmethod = "simplex" (current default) or 'subplex' (default of previous versions)
  idparsfix <- 2
  parsfix   <- 0
  options(warn = -1)
  namepars <- mbd::get_mbd_param_names()
  n_pars <- length(namepars); #if you add more parameters to your model just change this
  failpars <- rep(-1, n_pars); names(failpars) <- namepars; #those are the parameters that you get if something goes sideways
  if (is.numeric(brts) == FALSE) {
    cat("The branching times should be numeric.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
  } else {
    idpars <- sort(c(idparsopt, idparsfix))
    if ((sum(idpars == (1:n_pars)) != n_pars) ||
        (length(initparsopt) != length(idparsopt)) ||
        (length(parsfix) != length(idparsfix))
    ) {
      cat("The parameters to be optimized and / or fixed are incoherent.\n")
      out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    } else {
      if (length(namepars[idparsopt]) == 0) {
        optstr <- "nothing"
      } else {
        optstr <- namepars[idparsopt]
      }
      cat("You are optimizing", optstr, "\n")
      if (length(namepars[idparsfix]) == 0) {
        fixstr <- "nothing"
      } else {
        fixstr <- namepars[idparsfix]
      }
      cat("You are fixing", fixstr, "\n")
      cat("Optimizing the likelihood - this may take a while.", "\n")
      utils::flush.console()
      if (pars_transform == 1) {
        trparsopt <- initparsopt / (1 + initparsopt)
        trparsopt[which(initparsopt == Inf)] <- 1
        trparsfix <- parsfix / (1 + parsfix)
        trparsfix[which(parsfix == Inf)] <- 1
      } else {
        trparsopt <- initparsopt
        trparsfix <- parsfix
      }
      optimpars <- c(tol, maxiter)
      # there's no pars2 here and instead 3 more args at the end
      initloglik <- pmb_loglik_choosepar(
        trparsopt = trparsopt, trparsfix = trparsfix, idparsopt = idparsopt,
        idparsfix = idparsfix, brts = brts, soc = soc,
        pars_transform = pars_transform
      )
      cat("The loglikelihood for the initial parameter values is", initloglik, "\n")
      utils::flush.console()
      if (initloglik == -Inf) {
        cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
        out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
      } else {
        out <- DDD::optimizer(optimmethod = optimmethod, optimpars = optimpars,
                               fun = pmb_loglik_choosepar, trparsopt = trparsopt,
                               trparsfix = trparsfix, idparsopt = idparsopt,
                               idparsfix = idparsfix, brts = brts, soc = soc,
                               pars_transform = pars_transform)
        if (out$conv != 0) {
          cat("Optimization has not converged. Try again with different initial values.\n")
          out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
        } else {
          mltrpars <- as.numeric(unlist(out$par))
          if (pars_transform == 1) {
            ml_pars <- mltrpars / (1 - mltrpars)
          } else {
            ml_pars <- mltrpars
          }
          ml_pars1 <- rep(0, n_pars)
          names(ml_pars1) <- namepars
          ml_pars1[idparsopt] <- ml_pars
          if (length(idparsfix) != 0) {
            ml_pars1[idparsfix] <- parsfix
          }
          max_lik <- as.numeric(unlist(out$fvalues))
          out2 <- data.frame(
            t(ml_pars1),
            loglik = max_lik,
            df = length(initparsopt),
            conv = unlist(out$conv)
          )

          tobeprint <- "Maximum likelihood parameter estimates:"
          for (ii in 1:n_pars) {
            tobeprint <- paste(
              tobeprint, paste(names(ml_pars1[ii]), ":", sep = ""), ml_pars1[ii]
            )
          }
          s1 <- sprintf(tobeprint)

          if (out2$conv != 0 & changeloglikifnoconv == T) {
            out2$loglik <- -Inf
          }
          s2 <- sprintf("Maximum loglikelihood: %f", max_lik)
          cat("\n", s1, "\n", s2, "\n\n")
        }# bracket#5
      }# bracket#4
    }# bracket#3
  }# bracket#2
  invisible(out2)
}
