# mbd_ML----------------
#' @author Giovanni Laudanno
#' @title Maximization of the loglikelihood under a multiple birth-death diversification model
#' @description mbd_ML computes the maximum likelihood estimates of the parameters of a multiple birth-death diversification model for a given set of phylogenetic branching times. It also outputs the corresponding loglikelihood that can be used in model comparisons. Differently from mbd_ML it can account for three kind of events: sympatric (single) speciation, multiple (allopatric) speciation and extinction.
#' @inheritParams default_params_doc
#' @param initparsopt The initial values of the parameters that must be optimized
#' @param idparsfix The ids of the parameters that should not be optimized. The default is to fix all parameters not specified in idparsopt.
#' @param parsfix The values of the parameters that should not be optimized.
#' @param res Sets the maximum number of species for which a probability must be computed, must be larger than 1 + length(brts).
#' @param tol Sets the tolerances in the optimization. Consists of:
#' \itemize{
#' \item reltolx = relative tolerance of parameter values in optimization
#' \item reltolf = relative tolerance of function value in optimization
#' \item abstolx = absolute tolerance of parameter values in optimization
#' }
#' @param changeloglikifnoconv If TRUE the loglik will be set to -Inf if ML does not converge.
#' @param optimmethod Method used in optimization of the likelihood. Current default is 'simplex'. Alternative is 'subplex' (default of previous versions).
#' @param ... Something
#' @return The output is a dataframe containing estimated parameters and maximum
#' loglikelihood. The computed loglikelihood contains the factor q! m! / (q + m)!
#' where q is the number of species in the phylogeny and m is the number of
#' missing species, as explained in the supplementary material to Etienne et al. 2012.
#'
#' @examples
#' set.seed(10)
#' test_pars <- c(0.3, 0.1, 0.1, 0.15)
#' simulated_data = mbd_sim(pars = test_pars, soc = 2, age = 10, cond = 1)
#' graphics::plot(simulated_data$tes)
#' # @Giappo: does not work
#' # mbd_ML(
#' #   brts = simulated_data$brts, initparsopt = 0.11 ,idparsopt = 4,
#' #   idparsfix = c(1,2,3), parsfix = test_pars[idparsfix], 
#' #   missnumspec = 0, cond = 1, soc = 2
#' # )
#' @export
mbd_ML <- function(
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
  optimmethod = 'simplex', 
  methode = "expo",
  minimum_multiple_births = 0, 
  pars_transform = 1, 
  print_errors = 0, 
  verbose = TRUE, 
  ...
)
{# bracket#1
  # - tol = tolerance in optimization
  # - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
  # - maxiter = the maximum number of iterations in the optimization
  # - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
  # - optimmethod = 'simplex' (current default) or 'subplex' (default of previous versions)
  if (!is.numeric(brts))
  {
    stop("'brts' must be numeric")
  }
  if (length(idparsfix) == 0) {idparsfix <- NULL}
  if (missing(parsfix) && (length(idparsfix) == 0)){parsfix <- idparsfix <- NULL}

  options(warn=-1)
  namepars <- c("lambda", "mu", "nu", "q"); Npars <- length(namepars); #if you add more parameters to your model just change this
  failpars <- rep(-1, Npars); names(failpars) <- namepars; #those are the parameters that you get if something goes sideways
  if (is.numeric(brts) == FALSE)
  {
    cat("The branching times should be numeric.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    return(invisible(out2))
  }
  idpars <- sort(c(idparsopt, idparsfix))
  if ( (sum(idpars == (1:Npars)) != Npars) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)) )
  {
    cat("The parameters to be optimized and/or fixed are incoherent.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    return(out2)
  }
  
  if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
  if (verbose == TRUE) {
    cat("You are optimizing",optstr,"\n")
  }
  if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
  if (verbose == TRUE) {
    cat("You are fixing",fixstr,"\n")
    cat("Optimizing the likelihood - this may take a while.","\n")
    utils::flush.console()
  }
  if (pars_transform == 1)
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
  initloglik <- mbd_loglik_choosepar(trparsopt = trparsopt, trparsfix = trparsfix,
                                           idparsopt = idparsopt, idparsfix = idparsfix,
                                           brts = brts, missnumspec = missnumspec,
                                           cond = cond, soc = soc, tips_interval = tips_interval,
                                           methode = methode,
                                           minimum_multiple_births = minimum_multiple_births,
                                           pars_transform = pars_transform,
                                           print_errors = print_errors, ...) #there's no pars2 here and instead 3 more args at the end
  if (verbose == TRUE) {
    cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
    utils::flush.console()
  }
  if (initloglik == -Inf)
  {# bracket#4
    cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    return (invisible(out2))
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
  if (out$conv != 0)
  {# bracket#5
    cat("Optimization has not converged. Try again with different initial values.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    return(invisible(out2))
  }
  MLtrpars <- as.numeric(unlist(out$par))
  if (pars_transform == 1)
  {
    #Rampal's transformation
    MLpars = MLtrpars/(1 - MLtrpars)
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
  if (verbose == TRUE) {
    s1 <- sprintf(tobeprint)
  }
  if(out2$conv != 0 & changeloglikifnoconv == T) {out2$loglik <- -Inf}
  if (verbose == TRUE) {
    s2 <- sprintf('Maximum loglikelihood: %f',ML)
    cat("\n",s1,"\n",s2,"\n\n")
  }
  
  invisible(out2)
}

# mbd_ML_cluster----------------
# Moved to razzo


# pmb_ML_cluster----------------
# Moved to razzo
