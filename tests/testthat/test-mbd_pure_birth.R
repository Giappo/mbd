context("PureBirth")

test_that("PureBirth theoretical check", {
  skip("For @Giappo: cannot find function 'mbd_loglik0'")
  test_size=10; methode = "both" #"lsoda" "expo"
  #without the new Hanno's product I can check only up to some degree of accuracy
  #best possible accuracy is 0

  mbd_theoretical_loglik=function(pars,brts,soc=2){
    #this works only for mu=0
    lambda=pars[1]; mu=pars[2]; q=pars[3]
    if (mu!=0 | lambda<=0 | q<=0 | q>=1){
      return(-Inf)
    }
    #SETTING TIME INTERVALS AND NUMBER OF BIRTHS
    #     data=brts2time_intervals_and_births(brts)
    #     time_intervals=data$time_intervals
    #     births=data$births
    data=brts2time_intervals_and_births(brts)
    time_intervals=data$time_intervals
    births=data$births

    #LOGLIK COMPUTATION
    N0=soc
    k=N0+cumsum(c(0,births))
    A_term=1
    i=0:1e6
    for (t in 1:length(time_intervals)){
      poisson_term <- stats::dpois(i, lambda*time_intervals[t], log = FALSE)
      ii = i[poisson_term!=0]
      pois_not_zero = which(poisson_term!=0)
      A_term=A_term*sum( ((1-q)^(ii*k[t]))*poisson_term[pois_not_zero] )
      # poisson_term=stats::dpois(i, lambda*time_intervals[t], log = FALSE)[dpois(i, lambda*time_intervals[t], log = FALSE)!=0]
      # ii=i[stats::dpois(i, lambda*time_intervals[t], log = FALSE)!=0]
      # A_term=A_term*sum( ((1-q)^(ii*k[t]))*poisson_term )
    }

    B_term=prod( lambda*choose(k[-length(k)],births)*q^births*(1-q)^(k[-length(k)]-births) )

    loglik=log(A_term*B_term)
    # loglik=-loglik #Rampal's optimizer uses loglik rather than -loglik
    return(loglik)
  }

  mbd_test_pure_birth=function(pars,brts,soc=2,cond=0,methode="both"){

    if (pars[2]!=0){
      print("This is meant to work only for mu=0")
      break}

    theorethicalLL  = mbd_theoretical_loglik(pars=pars,brts = brts,soc = soc)
    computationalLL = mbd_loglik0(pars=pars,brts = brts,soc=soc, cond=cond,missnumspec=0,methode=methode)
    #I need to divide to get rid of common factors
    ctrl_pars=pars/2
    ctrl_brts=brts
    ctrl_theorethical=mbd_theoretical_loglik(pars=ctrl_pars,brts = ctrl_brts, soc = soc)
    ctrl_computational=mbd_loglik0(pars=ctrl_pars,brts = ctrl_brts,soc=soc, cond=cond,missnumspec=0,methode=methode)

    theorethicalLL  =  theorethicalLL  - ctrl_theorethical
    computationalLL =  computationalLL - ctrl_computational

    out=list(theorethicalLL=theorethicalLL,computationalLL=computationalLL)
    return(out)
  }

  lambda=sample(size = test_size*2,x = (1:18)*0.1,replace = T)
  q=sample(size = test_size*2,x = (1:18)*0.01,replace = T)
  theorethical_test=rep(0,test_size)
  computational_test=rep(0,test_size)

  for (j in 1:test_size){
    simpars = c(lambda[j],0,q[j]);
    brts = mbd:::mbd_sim0( pars=simpars,soc=2,age=10,cond=0 )$brts
    testpars=c(lambda[test_size+j],0,q[test_size+j])
    test_PB = mbd_test_pure_birth(testpars,brts,soc=2,cond=0,methode=methode)
    theorethical_test[j]=test_PB$theorethicalLL
    computational_test[j]=test_PB$computationalLL
  }
  # all.equal(theorethical_test, computational_test)
  # print(prod(theorethical_test))
  # print(prod(computational_test))

  for (j in 1:test_size) {
    expect_equal(theorethical_test[j], computational_test[j])
  }
})
