context("TestSims")

test_that("Sims are ok", {

  test_size=40;test_sim=vector("list")
  soc=2
  lambda=sample(size = test_size,x = (1:18)*0.1,replace = T)
  mu=sample(size = test_size,x = (1:18)*0.01,replace = T)
  q=sample(size = test_size,x = (1:16)*0.01,replace = T)
  
  for (j in 1:test_size){
    test_sim[[j]]=mbd::mbd_sim0(pars = c(lambda[j],mu[j],q[j]),soc = soc,age = 10,cond = 1)$L
  }
  
  #check if all the species are born before their parents die
  back_to_the_future_test = 1
  for (j in 1:test_size){
    L=test_sim[[j]]
    parents=abs(L[,2])
    back_to_the_future_test = back_to_the_future_test * prod( L[-c(1:soc),1]>=L[parents[-c(1:soc)],4] )
  }
  
  #check if the conditioning on survival works
  conditioning_test = 1
  for (j in 1:test_size){
    L=test_sim[[j]]
    check_signs = sort(unique(sign(L[,3])))
    conditioning_test = conditioning_test * prod ( check_signs == c(-1,1) )
  }

  expect_equal(back_to_the_future_test,1)
  expect_equal(conditioning_test,1)
  # warnings()
})

