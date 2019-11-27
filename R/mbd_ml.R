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
#' sim <- mbd::mbd_sim(
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
#'   true_pars = sim_pars,
#'   cond = cond,
#'   n_0 = n_0,
#'   verbose = FALSE
#' )
#' @export
mbd_ml <- function(
  loglik_function = mbd::mbd_loglik,
  brts,
  start_pars = c(0.5, 0.3, 0.5, 0.3),
  n_0 = 2,
  cond = 1,
  missnumspec = 0,
  lx = min(1 + 3 * (length(brts) + max(missnumspec)), mbd::max_lx()),
  ml_steps = 2,
  optim_ids = rep(TRUE, length(start_pars)),
  true_pars = start_pars,
  tips_interval = c(n_0 * (cond > 0), Inf),
  q_threshold = 1e-3,
  verbose = TRUE,
  maxiter = 10000
) {

  # setup and checks
  if (true_pars[3] == 0 | true_pars[4] == 0) {
    q_threshold <- 0
  }
  par_names <- mbd::get_param_names()
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

  # set up each ml_step
  abstol_base <- 1e-6
  reltol_base <- 1e-5
  lx_steps <- sort(ceiling(lx * (3 / 2) ^ (1:ml_steps - 1)))
  abstol_steps <- abstol_base * 10 ^ - (1:ml_steps - 1)
  reltol_steps <- reltol_base * 10 ^ - (1:ml_steps - 1)
  maxiter_steps <- rep(maxiter, ml_steps)
  out_list <- vector("list", ml_steps)

  # run maximum likelihood for each ml_step
  for (nn in 1:ml_steps) {

    abstol_nn <- abstol_steps[nn]
    reltol_nn <- reltol_steps[nn]
    maxiter_nn <- maxiter_steps[nn]
    lx_nn <- lx_steps[nn]

    # define function to optimize
    optim_fun_nn <- function(tr_optim_pars) {
      pars2 <- rep(0, length(start_pars))
      optim_pars <- mbd::pars_transform_back(tr_optim_pars)
      pars2[optim_ids] <- optim_pars
      pars2[!optim_ids] <- true_pars[!optim_ids]

      out <- -loglik_function(
        pars = pars2,
        brts = brts,
        cond = cond,
        n_0 = n_0,
        q_threshold = q_threshold,
        lx = lx_nn,
        missnumspec = missnumspec
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
    if (nn == 1) {
      tr_start_pars <- rep(0, length(start_pars))
      tr_start_pars <- mbd::pars_transform_forward(start_pars[optim_ids]) # nolint internal function
      initloglik <- -optim_fun_nn(tr_start_pars)
      utils::flush.console()
      if (initloglik == -Inf) {
        cat(
          message = "The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n" # nolint
        )
        out_final <- failout
        return(invisible(out_final))
      }
      tr_start_pars_nn <- tr_start_pars
    }

    # maximum likelihood n-th
    out_nn <- subplex::subplex(
      par = tr_start_pars_nn,
      fn = optim_fun_nn,
      control = list(
        abstol = abstol_nn,
        reltol = reltol_nn,
        maxit = maxiter_nn
      )
    )
    out_list[[nn]] <- out_nn

    # report missed convergence
    if (out_nn$conv > 0) {
      cat2(
        "Optimization has not converged. Try again with different initial values.\n", # nolint
        verbose = verbose
      )
      out_final <- data.frame(
        t(failpars),
        loglik = -1,
        df = -1,
        conv = unlist(out_nn$conv)
      )
      names(out_final) <- out_names
      return(invisible(out_final))
    }

    # update starting parameters
    tr_start_pars_nn <- out_nn$par
  }

  # store output
  likelihoods <- rep(NA, ml_steps)
  for (nn in 1:ml_steps) {
    temp <- utils::capture.output(
      likelihoods[nn] <- -optim_fun_nn(out_list[[nn]]$par)
    )
  }
  rm(temp)
  the_max <- max(which(likelihoods == max(likelihoods)))
  out <- out_list[[the_max]]

  # return mle results
  outpars <- rep(0, length(start_pars))
  outpars[optim_ids] <- mbd::pars_transform_back( # nolint internal function
    as.numeric(unlist(out$par))
  )
  outpars[!optim_ids] <- true_pars[!optim_ids]
  names(outpars) <- par_names

  out_final <- data.frame(
    row.names = NULL,
    outpars[1],
    outpars[2],
    outpars[3],
    outpars[4],
    -out$value,
    sum(optim_ids),
    unlist(out$conv)
  )
  names(out_final) <- out_names
  out_final
}
