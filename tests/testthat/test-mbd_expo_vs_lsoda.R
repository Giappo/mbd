context("ExpovsLsoda")
test_size=10;

soc=2;cond=1
lambda=sample(size = test_size*2,x = (1:10)*0.1,replace = T)
mu=sample(size = test_size*2,x = (1:10)*0.01,replace = T)
q=sample(size = test_size*2,x = (2:14)*0.01,replace = T)
test_expo=rep(0,test_size)
test_lsoda=rep(0,test_size)

for (j in 1:test_size){
simpars = c(lambda[j],mu[j],q[j]);
brts = MBD:::mbd_sim0( pars=simpars,soc=soc,age=10,cond=cond )$brts
testpars=c(lambda[test_size+j],mu[test_size+j],q[test_size+j])
test_expo[j]=MBD::mbd_loglik0(pars = testpars,brts = brts,soc = soc,cond = cond,methode="expo")
test_lsoda[j]=MBD::mbd_loglik0(pars = testpars,brts = brts,soc = soc,cond = cond,methode="lsoda")
}
# print(prod(test_lsoda))
# print(prod(test_expo) )

test_that("PureBirth theoretical check", {
  for (i in 1:test_size){
  expect_equal(test_lsoda[i],test_expo[i])
  }
})

# warnings()


