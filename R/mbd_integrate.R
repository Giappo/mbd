#' @useDynLib DDD
mbd_ode_FORTRAN <- function(
  initprobs,
  tvec,
  parsvec,
  atol,
  rtol,
  methode,
  runmod = "mbd_runmod"
)
{
  N <- length(initprobs)
  probs <- deSolve::ode(y = initprobs, parms = c(N + 0.), rpar = parsvec,
                        times = tvec, func = runmod, initfunc = "mbd_initmod",
                        ynames = c("SV"), dimens = N + 1, nout = 1, outnames = c("Sum"),
                        dllname = "mbd",atol = atol, rtol = rtol, method = methode)
  #[,1:(N + 1)]
  #if this fails then N + 2 instead of N + 1, or N?
  return(probs)
}

mbd_integrate <- function(
  y = y,
  times = tseq,
  func = func,
  parms = parms,
  atol = atol,
  rtol = rtol,
  tcrit = t1
)
{
  func_name <- 'no_name'
  if(is.character(func))
  {
    func_name <- func
  }
  if(func_name == 'no_name')
  {
    out <- deSolve::ode(
      y = y,
      times = times,
      func = func,
      parms = parms,
      atol = atol,
      rtol = rtol,
      tcrit = tcrit
    )
  } else
  {
    #func_name <- 'mbd_runmod'
    out <- mbd_ode_FORTRAN(
      initprobs = y,
      tvec = times,
      parsvec = parms,
      atol = atol,
      rtol = rtol,
      methode = 'ode45',
      runmod = func_name
    )
  }
  return(out)
}
