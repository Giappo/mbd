# mbd_ode_FORTRAN -----
#' Performs the integration of the ode using FORTRAN code
#' @inheritParams default_params_doc
#' @author Rampal S. Etienne
#' @export
#' @useDynLib mbd
mbd_ode_FORTRAN <- function(
  initprobs,
  tvec,
  parsvec,
  atol,
  rtol,
  methode,
  runmod = "mbd_runmod"
) {
  dimparsvec <- dim(parsvec)
  N <- length(initprobs)
  if (all(dimparsvec) > 1) {
    dim(parsvec) <- c(N ^ 2, 1)
  }
  if (runmod == "mbd_runmod") {
    initfunc <- "mbd_initmod"
  } else if (runmod == "mbd_runmodpcp" | runmod == "mbd_runmodpcq") {
    initfunc <- "mbd_initmodpc"
  }
  probs <- deSolve::ode(
    y = initprobs,
    parms = c(N + 0.),
    rpar = parsvec,
    times = tvec,
    func = runmod,
    initfunc = initfunc,
    ynames = c("SV"),
    dimens = N + 1,
    nout = 1,
    outnames = c("Sum"),
    dllname = "mbd",
    atol = atol,
    rtol = rtol,
    method = methode
  )[, 1:(N + 1)]

  return(probs)
}

#  mbd_integrate -----
#' performs integration of the ode
#' @inheritParams default_params_doc
#' @author Rampal S. Etienne
#' @export
mbd_integrate <- function(
  y,
  times,
  func,
  parms,
  atol,
  rtol,
  tcrit,
  methode = "lsoda"
) {
  func_name <- "no_name"
  if (is.character(func)) {
    func_name <- func
  }
  if (func_name == "no_name") {
    out <- deSolve::ode(
      y = y,
      times = times,
      func = func,
      parms = parms,
      atol = atol,
      rtol = rtol,
      tcrit = tcrit,
      method = methode
    )
  } else {
    out <- mbd_ode_FORTRAN(
      initprobs = y,
      tvec = times,
      parsvec = parms,
      atol = atol,
      rtol = rtol,
      methode = methode,
      runmod = func_name
    )
  }
  return(out)
}

#  mbd_solve -----
#' @title mbd ODE system integrator
#' @description Integrates "func" in the time interval
# *if* this function returns, the result doesn't contains
# any negative number
#' @inheritParams default_params_doc
#' @param func function for the right hand side of the ODE
#' @export
#' @author Hanno Hildenbrandt, adapted by Giovanni Laudanno
mbd_solve <- function(
  vector,
  time_interval,
  func = mbd_loglik_rhs,
  parms
) {

  y <- vector
  t1 <- time_interval

  g <- 10 # granularity
  t0 <- 0
  start_rtol <- 1e-8
  atol <- 1e-100 # realistically zero
  rtol <- start_rtol # something reasonable, hopefully
  while (TRUE) {
    tseq <- seq(t0, t1, length.out = g)
    out <- mbd_integrate(
      y = y,
      times = tseq,
      func = func,
      parms = parms,
      atol = atol,
      rtol = rtol,
      tcrit = t1
    )
    # it might be useful for debug istate = attributes(out)$istate
    # it might be useful for debug rstate = attributes(out)$rstate
    lkg <- 0 # last known good
    for (ff in 1:g) {
      a <- any(out[ff, -1] < 0)
      if (!is.na(a) && a) {
        break;
      }
      lkg <- lkg + 1
    }
    if (lkg == g) {
      break; # done and dusted
    }
    if (lkg > 1) {
      # trace back to last known good and try from there
      t0 <- as.numeric(out[lkg, 1])
      y <- as.numeric(out[lkg, -1])
      # relax tol to default
      rtol <- start_rtol
    } else {
      # no progress, make tol more strict
      rtol <- rtol / 100
    }
  }
  out[g, -1]
}
