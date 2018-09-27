context("PureBirthSmallTrees")

test_that("PureBirth theoretical check", {

  skip("Fix @Giappo: cannot find function 'mbd_loglik0'")
  test_size <- 10
  soc <- 2
  cond <- 0

  mbd_theoretical_loglik <- function(pars, brts, soc = 2) {
    #this works only for mu=0
    lambda <- pars[1]
    mu <- pars[2]
    q <- pars[3]
    if (mu != 0 | lambda <= 0 | q <= 0 | q >= 1) {
      return(-Inf)
    }
    #SETTING TIME INTERVALS AND NUMBER OF BIRTHS
    #     data=brts2time_intervals_and_births(brts)
    #     time_intervals=data$time_intervals
    #     births=data$births
    data <- brts2time_intervals_and_births(brts)
    time_intervals <- data$time_intervals
    births <- data$births

    #LOGLIK COMPUTATION
    init_n_lineages <- soc
    k <- init_n_lineages + cumsum(c(0, births))
    a_term <- 1
    i <- 0:1e6
    for (t in 1:length(time_intervals)) {
      poisson_term <- stats::dpois(i, lambda*time_intervals[t], log = FALSE)
      ii <- i[poisson_term != 0]
      pois_not_zero <- which(poisson_term != 0)
      a_term <- a_term * sum(((1 - q) ^ (ii * k[t])) * poisson_term[pois_not_zero] )
    }

    b_term <- prod(lambda*choose(k[-length(k)], births)*q^births *
      (1 - q) ^ (k[-length(k)] - births))

    loglik <- log(a_term*b_term)
    # loglik=-loglik #Rampal's optimizer uses loglik rather than -loglik
    loglik
  }

  mbd_test_pure_birth_small_trees <- function(
    pars,
    brts,
    soc = 2,
    cond = 0
  ) {
    if (pars[2] != 0) {
      print("This is meant to work only for mu=0")
      break
    }

    theorethicalLL <- mbd_theoretical_loglik(pars = pars, brts = brts, soc = soc)
    lsodaLL <- mbd_loglik0(
      pars = pars, brts = brts, soc = soc, cond = cond,
      missnumspec = 0, methode = "lsoda"
    )
    expoLL <- mbd_loglik0(
      pars = pars, brts = brts, soc = soc, cond = cond,
      missnumspec = 0, methode = "expo"
    )
    #I need to divide to get rid of common factors
    ctrl_pars <- pars / 2
    ctrl_brts <- brts

    ctrl_theorethical=mbd_theoretical_loglik(pars=ctrl_pars, brts = ctrl_brts, soc = soc)
    ctrl_lsoda=mbd_loglik0(pars=ctrl_pars, brts = ctrl_brts, soc=soc, cond=cond, missnumspec=0, methode="lsoda")
    ctrl_expo=mbd_loglik0(pars=ctrl_pars, brts = ctrl_brts, soc=soc, cond=cond, missnumspec=0, methode="expo")

    theorethicalLL  =  theorethicalLL  - ctrl_theorethical
    lsodaLL =  lsodaLL - ctrl_lsoda
    expoLL = expoLL - ctrl_expo

    out=list(theorethicalLL=theorethicalLL, lsodaLL=lsodaLL, expoLL=expoLL)
    return(out)
  }

  lambda=sample(size = test_size, x = (1:18)*0.1, replace = T)
  q=sample(size = test_size, x = (1:18)*0.01, replace = T)
  res_theorethical=rep(0, test_size)
  res_lsoda=rep(0, test_size)
  res_expo=rep(0, test_size)

  brts = c(10)
  for (j in 1:test_size){
    testpars=c(lambda[j],0, q[j])
    test_PB = mbd_test_pure_birth_small_trees(testpars, brts, soc=soc, cond=cond)
    res_theorethical[j] = test_PB$theorethicalLL
    res_lsoda[j]        = test_PB$lsodaLL
    res_expo[j]         = test_PB$expoLL
  }
  # all.equal(res_theorethical, res_expo)
  # all.equal(res_theorethical, res_lsoda)



  for (i in 1:test_size){
  expect_equal(res_theorethical[i], res_expo[i])
  expect_equal(res_theorethical[i], res_lsoda[i])
  }
})

# warnings()


