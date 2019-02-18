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
#' @return The output is a dataframe containing
#'   estimated parameters and maximum loglikelihood.
#'   The computed loglikelihood contains the factor q! m! / (q + m)!
#'   where q is the number of species in the phylogeny and m is the number of
#'   missing species, as explained in the supplementary
#'   material to Etienne et al. 2012.
#' @examples
#' set.seed(2)
#' lambda <- 0.2 # sympatric speciation rate
#' mu <- 0.15 # extinction rate;
#' nu <- 2.0 # multiple allopatric speciation trigger rate
#' q <- 0.1 # single-lineage speciation probability
#' sim_pars <- c(lambda, mu, nu, q)
#' crown_age <- 1
#' cond <- 1
#' n_0 <- 2
#' sim <- mbd_sim(
#'  pars = sim_pars,
#'  n_0 = n_0, # Use a crown age
#'  age = crown_age,
#'  cond = cond # Condition on non-extinction
#' )
#' start_pars <- c(0.2, 0.15, 2, 0.15)
#' optim_ids <- c(FALSE, FALSE, FALSE, TRUE)
#' graphics::plot(sim$reconstructed_tree)
#' out <- mbd::mbd_ml(
#'   start_pars = start_pars,
#'   optim_ids = optim_ids,
#'   brts = sim$brts,
#'   cond = cond,
#'   n_0 = n_0,
#'   verbose = FALSE
#' )
#' @export
mbd_ml <- function(
  loglik_function = mbd_loglik,
  brts,
  start_pars = c(0.5, 0.3, 0.5, 0.3),
  n_0 = 2,
  cond = 1,
  optim_ids = rep(TRUE, length(start_pars)),
  true_pars = start_pars,
  tips_interval = c(0, Inf),
  safety_threshold = 1e-3,
  verbose = TRUE,
  lx = 1 + 2 * (length(brts))
) {
  # setup and checks
  if (true_pars[3] == 0 | true_pars[4] == 0) {
    safety_threshold <- 0
  }
  par_names <- get_param_names() # nolint internal function
  testit::assert(length(optim_ids) == length(start_pars))
  testit::assert(length(true_pars) == length(start_pars))
  start_pars[!optim_ids] <- true_pars[!optim_ids]
  if (any(start_pars < 0)) {
    stop("You cannot start from negative parameters!")
  }
  out_names <- c(par_names, "loglik", "df", "conv")
  failpars <- rep(-1, length(start_pars))
  failout  <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
  colnames(failout) <- out_names

  # define function to optimize
  optim_fun <- function(tr_optim_pars) {
    pars2 <- rep(0, length(start_pars))
    optim_pars <- pars_transform_back(tr_optim_pars) # nolint internal function
    pars2[optim_ids] <- optim_pars
    pars2[!optim_ids] <- true_pars[!optim_ids]

    out <- -loglik_function(
      pars = pars2,
      brts = brts,
      cond = cond,
      n_0 = n_0,
      safety_threshold = safety_threshold,
      lx = lx
    )
    if (verbose == TRUE) {
      printed_values <- paste0(
        c(par_names, "loglik"),
        " = ",
        signif(c(pars2, -out), digits = 5)
      )
      print_this <- paste(printed_values, sep = ",")
      cat(print_this, "\n")
    }
    out
  }

  # initial likelihood
  tr_start_pars <- rep(0, length(start_pars))
  tr_start_pars <- pars_transform_forward(start_pars[optim_ids]) # nolint internal function
  initloglik <- -optim_fun(tr_start_pars)
  utils::flush.console()
  if (initloglik == -Inf) {
    cat(
      message = "The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n" # nolint
    )
    out2 <- failout
    return(invisible(out2))
  }

  # maximum likelihood
  out <- subplex::subplex(
    par = tr_start_pars,
    fn = function(x) optim_fun(x)
  )

  # report missed convergence
  if (out$conv > 0) {
    cat2(
      "Optimization has not converged. Try again with different initial values.\n", # nolint
      verbose = verbose
    )
    out2 <- data.frame(
      t(failpars),
      loglik = -1,
      df = -1,
      conv = unlist(out$conv)
    )
    names(out2) <- out_names
    return(invisible(out2))
  }

  # return mle results
  outpars <- rep(0, length(start_pars))
  outpars[optim_ids] <- pars_transform_back( # nolint internal function
    as.numeric(unlist(out$par))
  )
  outpars[!optim_ids] <- true_pars[!optim_ids]
  names(outpars) <- par_names

  out2 <- data.frame(
    row.names = NULL,
    outpars[1],
    outpars[2],
    outpars[3],
    outpars[4],
    -out$value,
    sum(optim_ids),
    unlist(out$conv)
  )
  names(out2) <- out_names
  return(out2)
}
