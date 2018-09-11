# #Here you can find and test all the major functions in the package MBD.
# #Run it line by line.
# #Problem1: There's a problem with "mbd_theoretical_pure_birth". Results are very close but not equal. I guess this is due to the inaccuracy of the matrix builder. We have to try again with the new version that Hanno will provide us.
# #Problem2: I can't maximize the likelihood using either "simplex" or "subplex". With the new matrix builder that Hanno is going to provide us we should (hopefully) get rid of this problem.
#
# #install the package
# # install.packages(file.choose(), repos = NULL, type = "source")
#
# #load data and settings
# data_path <- system.file(
#   "data", "sim_data", package = "MBD"
# )
# settings_path <- system.file(
#   "data", "general_settings", package = "MBD"
# )
# load(data_path     , envir=globalenv())
# load(settings_path , envir=globalenv())
# rm(data_path,settings_path)
# # load(file.choose(), envir=globalenv()) #use it for "data"
# # load(file.choose(), envir=globalenv()) #use it for "general settings"
# #You can also generate your own dataset with:
# #mbd_sim_dataset0(sim_pars=c(2.5,0.1,0.10),soc=2,cond=1,age=10,max_sims=1000,edge=Inf)
# #NB: If you use it, it will delete the previous one.
#
# #select an method of optimization
# optimmethod = "simplex"
# optimmethod = "subplex"
#
# #setup1
# brts = sim_data[[1]];initparsopt = sim_pars; idparsopt = 1:3; parsfix = NULL;
# idparsfix = NULL; soc = 2; cond = 1;trparsfix=parsfix; trparsopt=sim_pars[idparsopt]
# check_pars=rep(0,length(idparsfix)+length(idparsopt));check_pars[idparsfix]=trparsfix; check_pars[idparsopt]=trparsopt;check_pars
#
# #setup2
# brts = sim_data[[1]];initparsopt = sim_pars[3]; idparsopt = 3; parsfix = sim_pars[1:2];
# idparsfix = 1:2; soc = 2; cond = 1;trparsfix=parsfix;trparsopt=sim_pars[idparsopt]
# check_pars=rep(0,length(idparsfix)+length(idparsopt));check_pars[idparsfix]=trparsfix; check_pars[idparsopt]=trparsopt;check_pars
#
# #test the functions
# mbd_sim0(pars=sim_pars,soc=2,age=10,cond=1,tips_interval=c(0,Inf))
# mbd_loglik0(pars=sim_pars,brts = brts,soc=2,cond=1,tips_interval=c(0,Inf) )
# mbd_loglik_choosepar0(trparsfix = trparsfix, trparsopt = trparsopt, idparsopt = idparsopt,idparsfix = idparsfix, brts = brts)
# mbd_ML0(brts = brts,initparsopt = initparsopt, idparsopt = idparsopt, parsfix = parsfix, idparsfix = idparsfix, soc = soc, cond = cond,optimmethod = optimmethod)
# mbd_theoretical_pure_birth(pars=c(sim_pars[1],0,sim_pars[3]),brts=brts,soc=2,cond=0)
#
#
#
#
