#' @title Maximization of the loglikelihood
#'   under a multiple birth-death diversification model
#' @description mbd_ml computes the maximum likelihood estimates of
#'   the parameters of a multiple birth-death diversification model
#'   for a given set of phylogenetic branching times.
#'   It also outputs the corresponding loglikelihood
#'   that can be used in model comparisons.
#'   Differently from mbd_ml it can account for three kind of events:
#'   sympatric (single) speciation, multiple (allopatric) speciation
#'   and extinction.
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @param initparsopt The initial values of the parameters
#'   that must be optimized
#' @param idparsfix The ids of the parameters that should not be optimized.
#'   The default is to fix all parameters not specified in idparsopt.
#' @param parsfix The values of the parameters that should not be optimized.
#' @param res Sets the maximum number of species for which a probability
#'   must be computed, must be larger than 1 + length(brts).
#' @param tol Sets the tolerances in the optimization. Consists of:
#' \itemize{
#' \item reltolx = relative tolerance of parameter values in optimization
#' \item reltolf = relative tolerance of function value in optimization
#' \item abstolx = absolute tolerance of parameter values in optimization
#' }
#' @param changeloglikifnoconv If TRUE,
#'   the loglik will be set to -Inf if ML does not converge.
#' @param optimmethod Method used in optimization of the likelihood.
#'   Current default is 'simplex'.
#'   Alternative is 'subplex' (default of previous versions).
#' @param ... Something
#' @return The output is a dataframe containing
#'   estimated parameters and maximum loglikelihood.
#'   The computed loglikelihood contains the factor q! m! / (q + m)!
#'   where q is the number of species in the phylogeny and m is the number of
#'   missing species, as explained in the supplementary
#'   material to Etienne et al. 2012.
#'
#' @examples
#' set.seed(10)
#' test_pars <- c(0.3, 0.1, 0.1, 0.15)
#' simulated_data = mbd_sim(pars = test_pars, soc = 2, age = 10, cond = 1)
#' graphics::plot(simulated_data$tes)
#' # @Giappo: does not work
#' # mbd_ml(
#' #   brts = simulated_data$brts, initparsopt = 0.11 , idparsopt = 4,
#' #   idparsfix = c(1, 2, 3),
#' #   parsfix = test_pars[idparsfix],
#' #   missnumspec = 0, cond = 1, soc = 2
#' # )
#' @seealso use \link{pmb_ml} to perform a maximum likelihood estimation
#'   for a Pure Multiple Birth model.
#' @export
mbd_ml <- function(
  brts,
  initparsopt,
  idparsopt,
  idparsfix = (1:4)[-idparsopt],
  parsfix,
  missnumspec = 0,
  cond = 1,
  soc = 2,
  tips_interval = c(0, Inf),
  res = 10 * (1 + length(brts) + missnumspec),
  tol = c(1E-3, 1E-4, 1E-6),
  maxiter = 1000 * round((1.25)^length(idparsopt)),
  changeloglikifnoconv = FALSE,
  optimmethod = "simplex",
  methode = "expo",
  minimum_multiple_births = 0,
  pars_transform = 1,
  print_errors = 0,
  verbose = TRUE,
  ...
) 
{
  # - tol = tolerance in optimization
  # - changeloglikifnoconv = if TRUE the loglik
  #     will be set to -Inf if ML does not converge
  # - maxiter = the maximum number of iterations in the optimization
  # - changeloglikifnoconv = if TRUE
  #     the loglik will be set to -Inf if ML does not converge
  # - optimmethod = "simplex" (current default)
  #     or 'subplex' (default of previous versions)
  if (!is.numeric(brts)) {
    stop("'brts' must be numeric")
  }
  if (length(initparsopt) != length(idparsopt)) {
    stop("lengths of 'idparsopt' and'initparsopt' must match")
  }
  if (length(parsfix) != length(idparsfix)) {
    stop("lengths of 'idparsfix' and'parsfix' must match")
  }
  namepars <- mbd::get_mbd_param_names()
  n_pars <- length(namepars)
  if (!all(sort(c(idparsopt, idparsfix)) == 1:n_pars)) {
    stop(
      "IDs 1 to 4 must be present exactly once ",
      "in either 'idparsfix' or 'idparsopt'"
    )
  }

  if (length(idparsfix) == 0) {
    idparsfix <- NULL
  }
  if (missing(parsfix) && (length(idparsfix) == 0)) {
    parsfix <- idparsfix <- NULL
  }

  options(warn = -1)
  failpars <- rep(-1, n_pars)
  names(failpars) <- namepars
  if (length(namepars[idparsopt]) == 0) {
    optstr <- "nothing"
  } else {
    optstr <- namepars[idparsopt]
  }
  if (verbose == TRUE) {
    cat("You are optimizing", optstr, "\n")
  }
  if (length(namepars[idparsfix]) == 0) {
    fixstr <- "nothing"
  } else {
    fixstr <- namepars[idparsfix]
  }
  if (verbose == TRUE) {
    cat("You are fixing", fixstr, "\n")
    cat("Optimizing the likelihood - this may take a while.", "\n")
    utils::flush.console()
  }
  if (pars_transform == 1) {
    #Rampal's transformation
    trparsopt <- initparsopt / (1 + initparsopt)
    trparsopt[which(initparsopt == Inf)] <- 1
    trparsfix <- parsfix / (1 + parsfix)
    trparsfix[which(parsfix == Inf)] <- 1
  } else {
    trparsopt <- initparsopt
    trparsfix <- parsfix
  }
  optimpars  <- c(tol, maxiter)

  #there's no pars2 here and instead 3 more args at the end
  initloglik <- mbd_loglik_choosepar(
    trparsopt = trparsopt, trparsfix = trparsfix,
    idparsopt = idparsopt, idparsfix = idparsfix,
    brts = brts, missnumspec = missnumspec,
    cond = cond, soc = soc, tips_interval = tips_interval,
    methode = methode,
    minimum_multiple_births = minimum_multiple_births,
    pars_transform = pars_transform,
    print_errors = print_errors,
    ...
  )
  if (verbose == TRUE) {
    cat(
      "The loglikelihood for the initial parameter values is",
      initloglik, "\n"
    )
    utils::flush.console()
  }
  if (initloglik == -Inf) {
    warning(
      "The initial parameter values have a likelihood that is equal to 0",
      "or below machine precision. Try again with different initial values."
    )
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    return(invisible(out2))
  }
  if (verbose == TRUE) {
    sink(file = tempfile()) # Sink output here
  }
  out <- DDD::optimizer(
    optimmethod = optimmethod, optimpars = optimpars,
    fun = mbd_loglik_choosepar,
    trparsopt = trparsopt, trparsfix = trparsfix,
    idparsopt = idparsopt, idparsfix = idparsfix,
    brts = brts, missnumspec = missnumspec, cond = cond,
    soc = soc, tips_interval = tips_interval, methode = methode,
    minimum_multiple_births = minimum_multiple_births,
    pars_transform = pars_transform, print_errors = print_errors,
    ...
  )
  if (verbose == TRUE) {
    sink() # Give back the output
  }
  if (out$conv != 0) {
    warning("Optimization has not converged. ",
      "Try again with different initial values."
    )
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    return(invisible(out2))
  }
  mltrpars <- as.numeric(unlist(out$par))
  if (pars_transform == 1) {
    #Rampal's transformation
    ml_pars <- mltrpars / (1 - mltrpars)
  } else {
    ml_pars <- mltrpars
  }
  ml_pars1 <- rep(0, n_pars); names(ml_pars1) <- namepars
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
  if (verbose == TRUE) {
    s1 <- sprintf(tobeprint)
  }
  if (out2$conv != 0 & changeloglikifnoconv == TRUE) {
    out2$loglik <- -Inf
  }
  if (verbose == TRUE) {
    s2 <- sprintf("Maximum loglikelihood: %f", max_lik)
    cat("\n", s1, "\n", s2, "\n\n")
  }

  invisible(out2)
}

# mbd_ML_cluster----------------
#' @author Giovanni Laudanno
#' @title Maximization of the loglikelihood under a multiple birth-death diversification model, wrapped for cluster usage.
#' @description mbd_ML_cluster computes the maximum likelihood estimates of the parameters of a multiple birth-death diversification model. Differently from mbd_ML it needs only the number "s" of the simulations, making it suitable to run on a cluster. You will need two files to make it work: "general_settings","sim_data"; they are both generated by mbd_sim_dataset. This second version supports both sympatric and allopatric speciations.
#' This function is like Ulysses' bow: you cannot use it unless you are the owner of the package.
#' @inheritParams default_params_doc
#' @param s The number of the simulation you want to evaluate.
#' @param initparsopt The initial values of the parameters that must be optimized
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
#' # mbd_ML_cluster(1)
#'
#' @export
mbd_ML_cluster <- function(s, initparsopt = c(0.6, 0.1, 1.3, 0.16)){
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
  sink(file = paste0(simpath,"/errors/mbd_MLE_errors",s,".txt"), append = TRUE)
  
  res <- mbd::mbd_ML(brts = sim_data[[s]],
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
                     optimmethod = 'simplex',
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
  
  utils::write.table(matrix(out,ncol = length(out)),file = paste(simpath,"/mbd_MLE",s,".txt",sep = ''),append = T,row.names = F,col.names = F, sep = ",")
  if (res[1:4] != rep(-1, Npars)){suppressWarnings(  file.remove( paste0(simpath,"/errors/mbd_MLE_errors",s,".txt") )  )}
}

#' @title pmb_ML_cluster
#' @description Like mbd_ML_cluster but only for cases with mu=0.
#' This function is like Ulysses' bow: you cannot use it unless you are the owner of the package.
#' @inheritParams default_params_doc
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
  
  res <- mbd_ML(brts = sim_data[[s]],
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
                optimmethod = 'simplex',
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
  
  utils::write.table(matrix(out,ncol = length(out)),file = paste(simpath,"/mbd_MLE",s,".txt",sep = ''),append = T,row.names = F,col.names = F, sep = ",")
  if (res[1:4] != rep(-1, Npars)){suppressWarnings(  file.remove( paste(simpath,"/errors/mbd_MLE_errors",s,".txt",sep = '') )  )}
}