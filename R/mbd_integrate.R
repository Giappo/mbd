# mbd_ode_fortran -----
#' Performs the integration of the ode using FORTRAN code
#' @inheritParams default_params_doc
#' @author Rampal S. Etienne
#' @export
#' @useDynLib mbd
mbd_ode_fortran <- function(
  initprobs,
  tvec,
  parsvec,
  atol,
  rtol,
  methode,
  runmod = "mbd_runmod"
) {
  dimparsvec <- dim(parsvec)
  nn <- length(initprobs)
  if (all(dimparsvec) > 1) {
    dim(parsvec) <- c(nn ^ 2, 1)
  }
  if (runmod == "mbd_runmod") {
    initfunc <- "mbd_initmod"
  } else if (runmod == "mbd_runmodpcp" | runmod == "mbd_runmodpcq") {
    initfunc <- "mbd_initmodpc"
  }
  probs <- deSolve::ode(
    y = initprobs,
    parms = c(nn + 0.),
    rpar = parsvec,
    times = tvec,
    func = runmod,
    initfunc = initfunc,
    ynames = c("SV"),
    dimens = nn + 1,
    nout = 1,
    outnames = c("Sum"),
    dllname = "mbd",
    atol = atol,
    rtol = rtol,
    method = methode
  )
  probs <- probs[, 1:(nn + 1)]
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
    out <- mbd::mbd_ode_fortran(
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
# *if* this function returns, the result doesn't contain
# any negative number
#' @inheritParams default_params_doc
#' @param func function for the right hand side of the ODE
#' @author Hanno Hildenbrandt, adapted by Giovanni Laudanno
#' @export
mbd_solve <- function(
  vector,
  time_interval,
  func = mbd::mbd_loglik_rhs,
  parms
) {
  methodes <- c("lsoda", "ode45", "lsodes")
  i <- 1
  out <- NULL
  while (is.null(out) && i <= length(methodes)) {
    methode <- methodes[i]
    temp <- my_try_catch(
      mbd::mbd_solve2(
        vector = vector,
        time_interval = time_interval,
        func = func,
        parms = parms,
        methode = methode
      )
    )
    out <- temp$value
    i <- i + 1
  }
  out
}

#' @title mbd ODE system integrator
#' @description Integrates "func" in the time interval
# *if* this function returns, the result doesn't contain
# any negative number
#' @inheritParams default_params_doc
#' @param func function for the right hand side of the ODE
#' @author Hanno Hildenbrandt, adapted by Giovanni Laudanno
#' @export
mbd_solve2 <- function(
  vector,
  time_interval,
  func = mbd::mbd_loglik_rhs,
  parms,
  methode = "lsoda"
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
    out <- mbd::mbd_integrate(
      y = y,
      times = tseq,
      func = func,
      parms = parms,
      atol = atol,
      rtol = rtol,
      tcrit = t1,
      methode = methode
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
