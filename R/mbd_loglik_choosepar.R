#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
#' @export
mbd_loglik_choosepar  <- function(trparsopt, 
                                  trparsfix, 
                                  idparsopt = 1:4,
                                  idparsfix = (1:4)[-idparsopt], 
                                  brts, 
                                  cond = 1, 
                                  soc = 2,
                                  tips_interval = c(0, Inf), 
                                  missnumspec = 0, 
                                  methode = "expo",
                                  minimum_multiple_births = 0, 
                                  pars.transform = 0, 
                                  print_errors = 0){
  #This function provides a likelihood for a subset of parameters. This is built to work inside mbd_minusLL_vs_single_parameter or any optimizer like simplex, optim or subplex
  #idparsopt are the ids of the parameters you want to analyze
  #trparsopt are the values for parameters you want to analyze
  #idparsfix are the ids of the parameters you want to fix
  #trparsfix are the values for parameters you want to fix
  
  namepars <- c("lambda","mu","nu","q"); Npars <- length(namepars);
  if (length(trparsopt) == 4 && missing(trparsfix)){trparsfix <- NULL}
  trpars1 = rep(0, Npars)
  trpars1[idparsopt] <- trparsopt
  if (length(idparsfix) != 0)
  {
    trpars1[idparsfix] <- trparsfix
  }
  if ( min(trpars1[1:Npars]) < 0 ){loglik <- -Inf}else
  {
    if (pars.transform == 1)
    {
      #Rampal's transformation
      pars1 = trpars1/(1 - trpars1)
    }else
    {
      pars1 <- trpars1
    }
    loglik <- MBD:::mbd_loglik(pars = pars1, brts = brts, cond = cond, soc = soc,
                               tips_interval = tips_interval, methode = methode,
                               minimum_multiple_births = minimum_multiple_births, print_errors = print_errors)
  }
  if (is.nan(loglik) || is.na(loglik))
  {
    cat("There are parameter values used which cause numerical problems:",trpars1,"\n")
    loglik <- -Inf
  }
  return(loglik)
}