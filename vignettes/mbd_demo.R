## ------------------------------------------------------------------------
rm(list = ls()); options(warn = -1);
suppressWarnings( suppressPackageStartupMessages( library(MBD) ) )

## ------------------------------------------------------------------------
set.seed(14) #set.seed(42)
sim_pars <- c(0.2, 0.15, 2, 0.1); soc <- 2; cond <- 1; age <- 10;#sim_pars <- c(0.2, 0.1, 1.5, 0.1); soc = 2; cond = 1; age = 10;
sim <- MBD:::mbd_sim(pars = sim_pars, soc = soc, age = age, cond = cond, tips_interval = c(0, Inf))


## ---- fig.show='hold'----------------------------------------------------
plot(sim$tas)

## ---- fig.show='hold'----------------------------------------------------
plot(sim$tes)

## ------------------------------------------------------------------------
sim$L

## ------------------------------------------------------------------------
sim$brts

## ------------------------------------------------------------------------
test_pars <- c(0.3, 0.05, 1.0, 0.08)
mbd_loglik(pars = test_pars, brts = sim$brts, soc = soc, cond = cond, missnumspec = 0)

## ------------------------------------------------------------------------
# Uncomment if you are keen to wait.
# idparsopt <- 4
# ids <- 1:4; idparsfix <- ids[-idparsopt];
# parsfix <- sim_pars[idparsfix]; initparsopt <- 0.15;
# MBD:::mbd_ML(brts = sim$brts, initparsopt = initparsopt, idparsopt = idparsopt, parsfix = parsfix, idparsfix = idparsfix, soc = soc, cond = cond)


