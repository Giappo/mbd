#' @export
mbd_cluster_main <- function(sim_pars = c(0.5,0.1,0.3,0.15), soc = 2, cond = 1, age = 10, max_sims = 1000,
                             tips_interval = c(0, 100), edge = Inf, minimum_multiple_births = 0){
# Set up the cores -----
list.of.packages <- c("parallel", "snow")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)){install.packages(new.packages)}
library("snow"); library("parallel")
# Set path to mpi include files and library files -----
if(.Platform$OS.type == "unix")
{
  mpipath = dirname(dirname(Sys.which('mpirun')))
  mpiinc = paste('--with-Rmpi-include=',mpipath,'/include',sep='')
  mpilib = paste('--with-Rmpi-libpath=',mpipath,'/lib',sep='')
  if(!require("Rmpi"))install.packages("Rmpi", repos="http://cran.rstudio.com/", configure.args=c(mpiinc, mpilib, "--with-Rmpi-type=OPENMPI"))
  library("Rmpi")
}

# Set up my package -----
if (.Platform$OS.type == "windows")
{
  libs_dir <- dirname(getwd())
  results_dir <- paste0(dirname(getwd()),"/results")
  dir.create(results_dir, showWarnings = FALSE)
  # devtools::install_github(repo = "Giappo/mbd")
}
if(.Platform$OS.type == "unix")
{
  home_dir <- substring(getwd(),1,13)
  libs_dir <- paste0(home_dir,'/mbd_like/libs')
  results_dir <- paste0(dirname(getwd()),"/sims")
  dir.create(results_dir, showWarnings = FALSE)
  # lib_files <- list.files(pattern=paste0('[.]tar'), path=libs_dir, full.names=TRUE)
  # install.packages(lib_files, repos = NULL, dependencies = T)
  # devtools::install_github(repo = "Giappo/mbd") #MAYBE THE BEST OPTION FOR CLUSTER TOO???
}

# Counts the amount of CPUs minus the master CPU -----
# if (.Platform$OS.type == "unix")
# {
#   cpu <- (mpi.universe.size() - 1)
#   cl  <- makeMPIcluster(cpu)
# }
# if (.Platform$OS.type == "windows")
# {
  cpu <- detectCores() # Number of cores requested.
  hosts <- rep("localhost", cpu)
  cl <- makeCluster(hosts, type = "SOCK")
# }

# Set up sims -----
if (sim_pars[2] == 0){cond <- 0; tips_interval <- c(0, Inf)} #this allows me to use the analytical formula
# if (cond == 1){tips_interval[1] <- max(soc, tips_interval[1])}#if the tree is conditioned on the survival of crown species the minimum amount of tips has to be raised!!!
sim_parameter.list <- vector("list", max_sims); for (i in 1:max_sims){sim_parameter.list[[i]] <- sim_pars}

# Simulate trees -----
sim_res <- parLapply(cl = cl, #cluster
                     X = sim_parameter.list, #dataset
                     fun = MBD:::mbd_sim, #function you want to apply
                     soc = soc,
                     cond = cond,
                     age = age,
                     tips_interval = tips_interval,
                     minimum_multiple_births = minimum_multiple_births)

sim_data    <- lapply(sim_res, FUN = function(x) x$brts)
ext_species <- lapply(sim_res, FUN = function(x) x$extinct_species)
sim_tes     <- lapply(sim_res, FUN = function(x) x$tes)
sim_tas     <- lapply(sim_res, FUN = function(x) x$tas)
max_k <- (soc - 1) + (is.list(sim_data)) * max(sapply(sim_data, length)) + (1 - is.list(sim_data)) * length(sim_data)
if (is.list(sim_data))
{
  all_the_births <- sapply(sim_data, FUN = function(brts){return(MBD:::brts2time_intervals_and_births(brts)$births)})
}else
{
  all_the_births <- MBD:::brts2time_intervals_and_births(sim_data)$births
}
max_b <- max(unlist(all_the_births))
additional_species <- tips <- rep(0,max_sims)
for (s in 1:max_sims)
{
  additional_species[s] <- sum( duplicated(sim_data[[s]]) )
  tips[s] <- length(sim_data[[s]]) + 1
}

# Save simulation results -----
simpath  <- getwd()
datapath <- paste0(simpath,"/data")
dir.create(datapath, showWarnings = FALSE)
sim_data_name         <- paste0(datapath,"/sim_data")
general_settings_name <- paste0(datapath,"/general_settings")
sim_trees_name        <- paste0(datapath,"/sim_trees")
if (file.exists(sim_data_name)){suppressWarnings( file.remove(sim_data_name) )}
save(sim_data, file = sim_data_name)
if (file.exists(general_settings_name)){suppressWarnings( file.remove(general_settings_name) )}
save(sim_pars, soc, age, cond, max_sims, tips_interval, max_k, max_b, ext_species,
     additional_species, tips, minimum_multiple_births, file = general_settings_name)
if (file.exists(sim_trees_name)){suppressWarnings( file.remove(sim_trees_name) )}
save(sim_tas, sim_tes, file = sim_trees_name)

# Set up MLE -----
initparsopt <- c(0.5, 0.1, 1.7, 0.15)
parnames <- c("lambda","mu","nu","q"); Npars <- length(parnames)
idparsopt <- 1:Npars; parsfix <- NULL;

# Perform MLE -----
MLE_res <- parLapply(cl = cl, #cluster
                     X = sim_data, #dataset
                     fun = MBD::mbd_ML, #function you want to apply
                     initparsopt = initparsopt,
                     idparsopt = idparsopt,
                     idparsfix = (1:Npars)[-idparsopt],
                     parsfix = parsfix,
                     missnumspec = 0,
                     cond = cond,
                     soc = soc,
                     tips_interval = tips_interval,
                     tol = c(1E-3, 1E-4, 1E-6),
                     maxiter = 1000 * round((1.25)^length(idparsopt)),
                     changeloglikifnoconv = FALSE,
                     optimmethod = 'subplex',
                     minimum_multiple_births = minimum_multiple_births)

# Save MLE results -----
test.string  <- toString(sim_pars)
test.string2 <- gsub(", ", "-", test.string)
for (i in 1:4) {test.string2 <- sub("[.]", "_", test.string2)}
results_file_name <- paste0(test.string2,".RData")
save(MLE_res, file = paste0(results_dir, results_file_name))

# Close Stuff -----
gc() #garbage collector
stopCluster(cl)
mpi.exit()

return(sim_data)
}
