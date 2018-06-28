# mbd_ML----------------
#' @author Giovanni Laudanno
#' @title Maximization of the loglikelihood under a multiple birth-death diversification model
#' @description mbd_ML computes the maximum likelihood estimates of the parameters of a multiple birth-death diversification model for a given set of phylogenetic branching times. It also outputs the corresponding loglikelihood that can be used in model comparisons. Differently from mbd_ML it can account for three kind of events: sympatric (single) speciation, multiple (allopatric) speciation and extinction.
#' @param brts A set of branching times of a phylogeny.
#' @param initparsopt The initial values of the parameters that must be optimized
#' @param idparsopt The ids of the parameters that must be optimized. The ids are defined as follows:
#' \itemize{
#'   \item pars[1] is lambda, the sympatric speciation rate;
#'   \item pars[2] is mu, the extinction rate;
#'   \item pars[3] is nu, the multiple allopatric speciation trigger rate;
#'   \item pars[4] is q, the single-lineage speciation probability.
#' }
#' @param idparsfix The ids of the parameters that should not be optimized. The default is to fix all parameters not specified in idparsopt.
#' @param parsfix The values of the parameters that should not be optimized.
#' @param missnumspec The number of species that are in the clade but missing in the phylogeny.
#' @param cond Set 1 if you want to condition on stem or crown age and non-extinction of the phylogeny. Set 0 otherwise.
#' @param soc Sets whether stem or crown age should be used (1 or 2).
#' @param tips_interval It takes into account tips boundaries constrain on simulated dataset.
#' @param res Sets the maximum number of species for which a probability must be computed, must be larger than 1 + length(brts).
#' @param tol Sets the tolerances in the optimization. Consists of:
#' \itemize{
#' \item reltolx = relative tolerance of parameter values in optimization
#' \item reltolf = relative tolerance of function value in optimization
#' \item abstolx = absolute tolerance of parameter values in optimization
#' }
#' @param max_iter Sets the maximum number of iterations in the optimization.
#' @param changeloglikifnoconv If TRUE the loglik will be set to -Inf if ML does not converge.
#' @param optimmethod Method used in optimization of the likelihood. Current default is 'subplex'. Alternative is 'simplex' (default of previous versions).
#' @param methode Set "sexpm" if you want to use sexpm. Set "expo" if you want to use expoRkit. Set "lsoda" if you want to use "lsoda".
#' @return The output is a dataframe containing estimated parameters and maximum
#' loglikelihood. The computed loglikelihood contains the factor q! m! / (q + m)!
#' where q is the number of species in the phylogeny and m is the number of
#' missing species, as explained in the supplementary material to Etienne et al. 2012.
#'
#' @examples
#' set.seed(10)
#' test_pars <- c(0.3, 0.1, 0.1, 0.15)
#' simulated_data = MBD:::mbd_sim(pars = test_pars, soc = 2, age = 10, cond = 1)
#' plot(simulated_data$tes)
#' MBD:::mbd_ML(brts = simulated_data$brts, initparsopt = 0.11 ,idparsopt = 4,
#' idparsfix = c(1,2,3), parsfix = test_pars[idparsfix], missnumspec = 0, cond = 1, soc = 2)
#' @export
mbd_ML <- function(brts, initparsopt, idparsopt, idparsfix = (1:4)[-idparsopt], parsfix,
                   missnumspec = 0, cond = 1, soc = 2, tips_interval = c(0, Inf),
                   res = 10 * (1 + length(brts) + missnumspec), tol = c(1E-3, 1E-4, 1E-6),
                   maxiter = 1000 * round((1.25)^length(idparsopt)),
                   changeloglikifnoconv = FALSE, optimmethod = 'subplex', methode = "expo",
                   alpha = 10, minimum_multiple_births = 0, pars.transform = 1, print_errors = 0)
  {# bracket#1
  # - tol = tolerance in optimization
  # - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
  # - maxiter = the maximum number of iterations in the optimization
  # - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
  # - optimmethod = 'subplex' (current default) or 'simplex' (default of previous versions)

  if (missing(parsfix) && (length(idparsfix) == 0)){parsfix <- NULL}

  options(warn=-1)
  namepars <- c("lambda", "mu", "nu", "q"); Npars <- length(namepars); #if you add more parameters to your model just change this
  failpars <- rep(-1, Npars); names(failpars) <- namepars; #those are the parameters that you get if something goes sideways
  if (is.numeric(brts) == FALSE)
  {# bracket#2
    cat("The branching times should be numeric.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
  }else
  {# bracket#2
    idpars <- sort(c(idparsopt,idparsfix))
    if ( (sum(idpars == (1:Npars)) != Npars) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)) )
    {# bracket#3
      cat("The parameters to be optimized and/or fixed are incoherent.\n")
      out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    }else
    {# bracket#3
      if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
      cat("You are optimizing",optstr,"\n")
      if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
      cat("You are fixing",fixstr,"\n")
      cat("Optimizing the likelihood - this may take a while.","\n")
      flush.console()
      if (pars.transform == 1)
      {
        #Rampal's transformation
        trparsopt = initparsopt/(1 + initparsopt)
        trparsopt[which(initparsopt == Inf)] = 1
        trparsfix = parsfix/(1 + parsfix)
        trparsfix[which(parsfix == Inf)] = 1
      }else
      {
        trparsopt  <- initparsopt
        trparsfix  <- parsfix
      }
      optimpars  <- c(tol, maxiter)
      initloglik <- MBD:::mbd_loglik_choosepar(trparsopt = trparsopt, trparsfix = trparsfix,
                                               idparsopt = idparsopt, idparsfix = idparsfix,
                                               brts = brts, missnumspec = missnumspec,
                                               cond = cond, soc = soc, tips_interval = tips_interval,
                                               methode = methode, alpha = alpha,
                                               minimum_multiple_births = minimum_multiple_births,
                                               pars.transform = pars.transform,
                                               print_errors = print_errors) #there's no pars2 here and instead 3 more args at the end
      cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
      flush.console()
      if (initloglik == -Inf)
      {# bracket#4
        cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
        out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
      }else
      {# bracket#4
        out <- DDD:::optimizer(optimmethod = optimmethod, optimpars = optimpars,
                               fun = MBD:::mbd_loglik_choosepar,
                               trparsopt = trparsopt, trparsfix = trparsfix,
                               idparsopt = idparsopt, idparsfix = idparsfix,
                               brts = brts, missnumspec = missnumspec, cond = cond,
                               soc = soc, tips_interval = tips_interval, methode = methode,
                               alpha = alpha, minimum_multiple_births = minimum_multiple_births,
                               pars.transform = pars.transform, print_errors = print_errors)
        if (out$conv != 0)
        {# bracket#5
          cat("Optimization has not converged. Try again with different initial values.\n")
          out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
        }else
        {# bracket#5
          MLtrpars <- as.numeric(unlist(out$par))
          if (pars.transform == 1)
          {
            #Rampal's transformation
            MLpars = MLtrpars/(1-MLtrpars)
          }else
          {
            MLpars <- MLtrpars
          }
          MLpars1 <- rep(0, Npars); names(MLpars1) <- namepars
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

# mbd_ML0----------------
#' @author Giovanni Laudanno
#' @title Maximization of the loglikelihood under a multiple birth-death diversification model
#' @description mbd_ML0 computes the maximum likelihood estimates of the parameters of a multiple birth-death diversification model for a given set of phylogenetic branching times. It also outputs the corresponding loglikelihood that can be used in model comparisons.
#' @param brts A set of branching times of a phylogeny.
#' @param initparsopt The initial values of the parameters that must be optimized
#' @param idparsopt The ids of the parameters that must be optimized. The ids are defined as follows:
#' \itemize{
#' \item id == 1 corresponds to lambda (multiple speciation trigger rate)
#' \item id == 2 corresponds to mu (extinction rate)
#' \item id == 3 corresponds to q (single-lineage speciation probability)
#' }
#' @param idparsfix The ids of the parameters that should not be optimized. The default is to fix all parameters not specified in idparsopt.
#' @param parsfix The values of the parameters that should not be optimized.
#' @param missnumspec The number of species that are in the clade but missing in the phylogeny.
#' @param cond Set 1 if you want to condition on stem or crown age and non-extinction of the phylogeny. Set 0 otherwise.
#' @param soc Sets whether stem or crown age should be used (1 or 2).
#' @param tips_interval It takes into account tips boundaries constrain on simulated dataset.
#' @param res Sets the maximum number of species for which a probability must be computed, must be larger than 1 + length(brts).
#' @param tol Sets the tolerances in the optimization. Consists of:
#' \itemize{
#' \item reltolx = relative tolerance of parameter values in optimization
#' \item reltolf = relative tolerance of function value in optimization
#' \item abstolx = absolute tolerance of parameter values in optimization
#' }
#' @param max_iter Sets the maximum number of iterations in the optimization.
#' @param changeloglikifnoconv If TRUE the loglik will be set to -Inf if ML does not converge.
#' @param optimmethod Method used in optimization of the likelihood. Current default is 'subplex'. Alternative is 'simplex' (default of previous versions).
#' @param methode Set "sexpm" if you want to use sexpm. Set "expo" if you want to use expoRkit. Set "lsoda" if you want to use "lsoda".
#' @return The output is a dataframe containing estimated parameters and maximum
#' loglikelihood. The computed loglikelihood contains the factor q! m! / (q + m)!
#' where q is the number of species in the phylogeny and m is the number of
#' missing species, as explained in the supplementary material to Etienne et al. 2012.
#'
#' @examples
#' set.seed(11)
#' test_pars = c(1.6,0.1,0.08)
#' simulated_data = MBD:::mbd_sim0( pars=test_pars,soc=2,age=10,cond=1 )
#' plot(simulated_data$tas)
#' MBD:::mbd_ML0(brts=simulated_data$brts, initparsopt = 0.11 ,idparsopt = 3,
#' idparsfix = 1:2 ,parsfix = test_pars[1:2],missnumspec=0,cond=1, soc = 2)
#' @export
mbd_ML0 <- function(brts, initparsopt, idparsopt, idparsfix = (1:3)[-idparsopt],
                    parsfix, missnumspec = 0, cond = 1, soc = 2, tips_interval=c(0,Inf),
                    res = 10 * (1 + length(brts) + missnumspec), tol = c(1E-3, 1E-4, 1E-6),
                    maxiter = 1000 * round((1.25)^length(idparsopt)), changeloglikifnoconv = FALSE,
                    optimmethod = 'subplex', methode = "expo", alpha = 20, pars.transform = 1)
{# bracket#1
  # - tol = tolerance in optimization
  # - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
  # - maxiter = the maximum number of iterations in the optimization
  # - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
  # - optimmethod = 'subplex' (current default) or 'simplex' (default of previous versions)

  if (missing(parsfix) && (length(idparsfix)==0)){parsfix <- NULL}

  options(warn=-1)
  namepars <- c("lambda","mu","q"); Npars <- length(namepars); #if you add more parameters to your model just change this
  failpars <- rep(-1,Npars); names(failpars) <- namepars; #those are the parameters that you get if something goes sideways
  if (is.numeric(brts) == FALSE)
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
      if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
      cat("You are optimizing",optstr,"\n")
      if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
      cat("You are fixing",fixstr,"\n")
      cat("Optimizing the likelihood - this may take a while.","\n")
      flush.console()
      if (pars.transform == 1)
      {
        #Rampal's transformation
        trparsopt = initparsopt/(1 + initparsopt)
        trparsopt[which(initparsopt == Inf)] = 1
        trparsfix = parsfix/(1 + parsfix)
        trparsfix[which(parsfix == Inf)] = 1
      }else
      {
        trparsopt  <- initparsopt
        trparsfix  <- parsfix
      }
      optimpars  <- c(tol, maxiter)
      initloglik <- MBD:::mbd_loglik_choosepar0(trparsopt = trparsopt, trparsfix = trparsfix,
                                                idparsopt = idparsopt, idparsfix = idparsfix,
                                                brts = brts, missnumspec = missnumspec,
                                                cond = cond, soc = soc,
                                                tips_interval = tips_interval, methode = methode,
                                                alpha = alpha, pars.transform = pars.transform) #there's no pars2 here and instead 3 more args at the end
      cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
      flush.console()
      if(initloglik == -Inf)
      {# bracket#4
        cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
        out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
      } else {
        out <- DDD:::optimizer(optimmethod = optimmethod, optimpars = optimpars,
                               fun = MBD:::mbd_loglik_choosepar0, trparsopt = trparsopt,
                               trparsfix = trparsfix, idparsopt = idparsopt,
                               idparsfix = idparsfix, brts = brts, missnumspec = missnumspec,
                               cond = cond, soc = soc, tips_interval = tips_interval,
                               methode = methode, alpha = alpha, pars.transform = pars.transform)
        if(out$conv != 0)
        {# bracket#5
          cat("Optimization has not converged. Try again with different initial values.\n")
          out2 = data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
        } else {
          MLtrpars = as.numeric(unlist(out$par))
          if (pars.transform == 1)
          {
            #Rampal's transformation
            MLpars = MLtrpars/(1-MLtrpars)
          }else
          {
            MLpars <- MLtrpars
          }
          MLpars1 <- rep(0, Npars); names(MLpars1) <- namepars
          MLpars1[idparsopt] <- MLpars
          if(length(idparsfix) != 0) {MLpars1[idparsfix] <- parsfix}
          ML <- as.numeric(unlist(out$fvalues))
          out2 <- data.frame(t(MLpars1), loglik = ML, df = length(initparsopt), conv = unlist(out$conv))

          tobeprint <- "Maximum likelihood parameter estimates:"
          for (ii in 1:Npars)
          {
            tobeprint <- paste(tobeprint,paste(names(MLpars1[ii]),":",sep = ""),MLpars1[ii])
          }
          s1 <- sprintf(tobeprint)

          if(out2$conv != 0 & changeloglikifnoconv == T) { out2$loglik = -Inf }
          s2 = sprintf('Maximum loglikelihood: %f',ML)
          cat("\n",s1,"\n",s2,"\n\n")
        }# bracket#5
      }# bracket#4
    }# bracket#3
  }# bracket#2
  invisible(out2)
}# bracket#1

# mbd_ML_cluster----------------
#' @author Giovanni Laudanno
#' @title Maximization of the loglikelihood under a multiple birth-death diversification model, wrapped for cluster usage.
#' @description mbd_ML_cluster computes the maximum likelihood estimates of the parameters of a multiple birth-death diversification model. Differently from mbd_ML it needs only the number "s" of the simulations, making it suitable to run on a cluster. You will need two files to make it work: "general_settings","sim_data"; they are both generated by mbd_sim_dataset. This second version supports both sympatric and allopatric speciations.
#' @param s The number of the simulation you want to evaluate.
#' @param initparsopt The initial values of the parameters that must be optimized
#' @param idparsopt The ids of the parameters that must be optimized. The ids are defined as follows:
#' \itemize{
#' \item id == 1 corresponds to lambda (speciation rate)
#' \item id == 2 corresponds to mu (extinction rate)
#' \item id == 3 corresponds to nu (multiple speciation trigger rate)
#' \item id == 4 corresponds to q (single-lineage speciation probability)
#' }
#' @param parsfix The values of the parameters that should not be optimized.
#' @return The output is saved on the document "mbd_MLE.txt".
#' \itemize{
#' \item First column contains ML estimates for lambda.
#' \item Second column contains ML estimates for mu.
#' \item Third column contains ML estimates for lambda.
#' \item Fourth column contains maximum likelihood.
#' \item Fifth column contains the number of additional species coming from multiple births in the evaluated tree.
#' }
#'
#' @examples
#' #You will need two files to make it work: "general_settings","sim_data".
#' MBD:::mbd_ML_cluster(1)
#'
#' @export
mbd_ML_cluster <- function(s, initparsopt = c(0.5, 0.1, 1.7, 0.15)){
  print(s)
  parnames <- c("lambda","mu","nu","q"); Npars <- length(parnames)
  idparsopt <- 1:Npars; parsfix <- NULL;

  #load general sim settings in order to make the estimation coherent with sim conditions
  # simpath=paste("sims/",sim_pars[1],"-",sim_pars[2],"-",sim_pars[3],"/",sep = '')
  simpath  <- getwd()
  datapath <- paste0(simpath,"/data")
  load(file = paste0(datapath,"/general_settings"))
  load(file = paste0(datapath,"/sim_data"))

  if (!file.exists(paste0(simpath,"/errors"))){dir.create(paste0(simpath,"/errors"))}
  sink(file = paste0(simpath,"/errors/mbd_MLE_errors",s,".txt"), append = T)

  res <- MBD:::mbd_ML(brts = sim_data[[s]],
                      initparsopt = initparsopt,
                      idparsopt = idparsopt,
                      idparsfix = (1:Npars)[-idparsopt],
                      parsfix = parsfix,
                      missnumspec = 0,
                      cond = cond,
                      soc = soc,
                      tips_interval = tips_interval,
                      res = 10 * (1 + length(brts) + missnumspec),
                      tol = c(1E-3, 1E-4, 1E-6),
                      maxiter = 1000 * round((1.25)^length(idparsopt)),
                      changeloglikifnoconv = FALSE,
                      optimmethod = 'subplex',
                      minimum_multiple_births = minimum_multiple_births)

  #additional tree info
  how_many_multiple <- percent_multiple <- -1;
  if (length(sim_data[[s]]) > 2)
  {
    test0 <- sim_data[[s]][-1]; #actual branching events
    test1 <- duplicated( test0 ); #additional species
    test2 <- test1; for (iii in 2:length(test1)){if(test1[iii] == T){test2[iii - 1] = T}} #considering also the first species at each multiple event
    how_many_multiple <- sum(test2);
    percent_multiple  <- how_many_multiple/length(test2);
  }
  # additional_species = sum( duplicated(sim_data[[s]]) );
  tips <- length(sim_data[[s]]) + 1;

  out  <- c(res[1:(Npars + 1)], how_many_multiple, tips, percent_multiple, s);
  out2 <- out; #names(out2)=c("lambda","mu","nu","q","LL","species born from multiple events","number of tips","percentage of species born from multiple events","tree id")
  names(out2) <- c(parnames,"LL","species born from multiple events","number of tips","percentage of species born from multiple events","tree id")
  print(out2)
  #out[1] = lambda
  #out[2] = mu
  #out[3] = nu
  #out[4] = q
  #out[5] = LL
  #out[6] = species born from multiple events
  #out[7] = number of tips
  #out[8] = percentage of species born from multiple events
  #out[9] = tree id
  sink()
  print(out2)

  write.table(matrix(out,ncol = length(out)),file = paste(simpath,"/mbd_MLE",s,".txt",sep = ''),append = T,row.names = F,col.names = F, sep = ",")
  if (res[1:4] != rep(-1, Npars)){suppressWarnings(  file.remove( paste(simpath,"/errors/mbd_MLE_errors",s,".txt",sep = '') )  )}
}

# pmb_ML_cluster----------------
#' @export
pmb_ML_cluster <- function(s, initparsopt = c(0.5, 0, 1.7, 0.15)){
  print(s)
  parnames <- c("lambda","mu","nu","q"); Npars <- length(parnames)
  idparsopt <- 1:Npars; parsfix <- NULL;
  
  #load general sim settings in order to make the estimation coherent with sim conditions
  # simpath=paste("sims/",sim_pars[1],"-",sim_pars[2],"-",sim_pars[3],"/",sep = '')
  simpath  <- getwd()
  datapath <- paste0(simpath,"/data")
  load(file = paste0(datapath,"/general_settings"))
  load(file = paste0(datapath,"/sim_data"))
  if (sim_pars[2] == 0 && cond == 0 && tips_interval == c(0, Inf))
  { # pure birth adjustments
    idparsopt2 <- idparsopt[idparsopt != 2]
    idparsfix2 <- (1:4)[-idparsopt2]
    parsfix2 <- initparsopt2 <- rep(NA,4);
    parsfix2[idparsfix] <- parsfix
    parsfix2[2] <- 0
    parsfix2 <- parsfix2[!is.na(parsfix2)]
    initparsopt2[idparsopt] <- initparsopt
    initparsopt2 <- initparsopt2[c(1,3,4)]
    initparsopt2 <- initparsopt2[!is.na(initparsopt2)]

    parsfix     <- parsfix2
    idparsfix   <- idparsfix2
    idparsopt   <- idparsopt2
    initparsopt <- initparsopt2
  }
  if (sim_pars[2] != 0) {stop("You should not use this function if mu > 0")}
  
  if (!file.exists(paste0(simpath,"/errors"))){dir.create(paste0(simpath,"/errors"))}
  sink(file = paste0(simpath,"/errors/mbd_MLE_errors",s,".txt"), append = T)
  
  res <- MBD:::mbd_ML(brts = sim_data[[s]],
                      initparsopt = initparsopt,
                      idparsopt = idparsopt,
                      idparsfix = (1:Npars)[-idparsopt],
                      parsfix = parsfix,
                      missnumspec = 0,
                      cond = cond,
                      soc = soc,
                      tips_interval = tips_interval,
                      res = 10 * (1 + length(brts) + missnumspec),
                      tol = c(1E-3, 1E-4, 1E-6),
                      maxiter = 1000 * round((1.25)^length(idparsopt)),
                      changeloglikifnoconv = FALSE,
                      optimmethod = 'subplex',
                      minimum_multiple_births = minimum_multiple_births)
  
  #additional tree info
  how_many_multiple <- percent_multiple <- -1;
  if (length(sim_data[[s]]) > 2)
  {
    test0 <- sim_data[[s]][-1]; #actual branching events
    test1 <- duplicated( test0 ); #additional species
    test2 <- test1; for (iii in 2:length(test1)){if(test1[iii] == T){test2[iii - 1] = T}} #considering also the first species at each multiple event
    how_many_multiple <- sum(test2);
    percent_multiple  <- how_many_multiple/length(test2);
  }
  # additional_species = sum( duplicated(sim_data[[s]]) );
  tips <- length(sim_data[[s]]) + 1;
  
  out  <- c(res[1:(Npars + 1)], how_many_multiple, tips, percent_multiple, s);
  out2 <- out; #names(out2)=c("lambda","mu","nu","q","LL","species born from multiple events","number of tips","percentage of species born from multiple events","tree id")
  names(out2) <- c(parnames,"LL","species born from multiple events","number of tips","percentage of species born from multiple events","tree id")
  print(out2)
  #out[1] = lambda
  #out[2] = mu
  #out[3] = nu
  #out[4] = q
  #out[5] = LL
  #out[6] = species born from multiple events
  #out[7] = number of tips
  #out[8] = percentage of species born from multiple events
  #out[9] = tree id
  sink()
  print(out2)
  
  write.table(matrix(out,ncol = length(out)),file = paste(simpath,"/mbd_MLE",s,".txt",sep = ''),append = T,row.names = F,col.names = F, sep = ",")
  if (res[1:4] != rep(-1, Npars)){suppressWarnings(  file.remove( paste(simpath,"/errors/mbd_MLE_errors",s,".txt",sep = '') )  )}
}

# mbd_ML_cluster0----------------
#' @author Giovanni Laudanno
#' @title Maximization of the loglikelihood under a multiple birth-death diversification model, wrapped for cluster usage.
#' @description mbd_ML_cluster0 computes the maximum likelihood estimates of the parameters of a multiple birth-death diversification model. Differently from mbd_ML it needs only the number "s" of the simulations, making it suitable to run on a cluster. You will need two files to make it work: "general_settings","sim_data"; they are both generated by mbd_sim_dataset0.
#' @param s The number of the simulation you want to evaluate.
#' @param initparsopt The initial values of the parameters that must be optimized
#' @param idparsopt The ids of the parameters that must be optimized. The ids are defined as follows:
#' \itemize{
#' \item id == 1 corresponds to lambda (multiple speciation trigger rate)
#' \item id == 2 corresponds to mu (extinction rate)
#' \item id == 3 corresponds to q (single-lineage speciation probability)
#' }
#' @param parsfix The values of the parameters that should not be optimized.
#' @return The output is saved on the document "mbd_MLE.txt".
#' \itemize{
#' \item First column contains ML estimates for lambda.
#' \item Second column contains ML estimates for mu.
#' \item Third column contains ML estimates for lambda.
#' \item Fourth column contains maximum likelihood.
#' \item Fifth column contains the number of additional species coming from multiple births in the evaluated tree.
#' }
#'
#' @examples
#' #You will need two files to make it work: "general_settings","sim_data".
#' MBD:::mbd_ML_cluster0(1)
#'
#' @export
mbd_ML_cluster0 <- function(s, initparsopt = c(1.8,0.3,0.15)){
  # initparsopt=c(1.8,0.3,0.15);
  parnames=c("lambda","mu","nu","q"); Npars = length(parnames)
  idparsopt=1:Npars;parsfix=NULL;

  # simpath=paste("sims/",sim_pars[1],"-",sim_pars[2],"-",sim_pars[3],"/",sep = '')
  simpath = getwd()

  datapath=paste(simpath,"/data",sep = '')
  load(file=paste(datapath,"/general_settings",sep = ''))
  load(file=paste(datapath,"/sim_data",sep = ''))
  print(s)

  if ( !file.exists(paste(simpath,"/errors",sep = '')) ){dir.create(paste(simpath,"/errors",sep = ''))}
  sink(file = paste(simpath,"/errors/mbd_MLE_errors",s,".txt",sep = ''), append = T)

  res <- MBD:::mbd_ML0(brts=sim_data[[s]],
         initparsopt=initparsopt,
         idparsopt=idparsopt,
         idparsfix = (1:Npars)[-idparsopt],
         parsfix = parsfix,
         missnumspec=0,
         cond= cond,
         soc = soc,
         tips_interval=tips_interval,
         res = 10*(1+length(brts)+missnumspec),
         tol = c(1E-3, 1E-4, 1E-6),
         maxiter = 1000 * round((1.25)^length(idparsopt)),
         changeloglikifnoconv = FALSE,
         optimmethod = 'subplex')

  #additional tree info
  how_many_multiple=percent_multiple=-1;
  if (length(sim_data[[s]])>2){
    test0=sim_data[[s]][-1]; #actual branching events
    test1=duplicated( test0 ); #additional species
    test2=test1;for (iii in 2:length(test1)){if(test1[iii]==T){test2[iii-1]=T}} #considering also the first species at each multiple event
    how_many_multiple=sum(test2);
    percent_multiple <- how_many_multiple/length(test2);
  }
  # additional_species = sum( duplicated(sim_data[[s]]) );
  tips=length(sim_data[[s]])+1;

  out=c(res[1:(Npars+1)],how_many_multiple,tips,percent_multiple,s);
  out2=out;#names(out2)=c("lambda","mu","q","LL","species born from multiple events","number of tips","percentage of species born from multiple events","tree id")
  names(out2)=c(parnames,"LL","species born from multiple events","number of tips","percentage of species born from multiple events","tree id")
  print(out2)
  #out[1] = lambda
  #out[2] = mu
  #out[3] = q
  #out[4] = LL
  #out[5] = species born from multiple events
  #out[6] = number of tips
  #out[7] = percentage of species born from multiple events
  #out[8] = tree id
  sink()
  print(out2)

  write.table(matrix(out,ncol = length(out)),file = paste(simpath,"/mbd_MLE",s,".txt",sep = ''),append = T,row.names = F,col.names = F, sep = ",")
  if (res[1:4]!=c(-1,-1,-1,-1)){suppressWarnings(  file.remove( paste(simpath,"/errors/mbd_MLE_errors",s,".txt",sep = '') )  )}
}

