# devtools::install_github("richelbilderbeek/tracerer", dependencies = TRUE)
# devtools::install_github("richelbilderbeek/beastier", dependencies = TRUE)
# devtools::install_github("richelbilderbeek/babette", dependencies = TRUE)
# devtools::install_github("richelbilderbeek/pirouette", dependencies = TRUE)
# devtools::install_github("richelbilderbeek/ribir", dependencies = TRUE)
# install.packages('xmlparsedata')
# install.packages('lintr')
# devtools::install_github("richelbilderbeek/raket", dependencies = TRUE)
install.packages("nLTT")
library("ape"); library("raket"); library("nLTT")

compare.MBD.with.alignments.BEAST.posterior <- function(sim_phylo,
                                                        chain_length = 1e+07,
                                                        sample_interval = 1e+03,
                                                        sequence_length = 1e+03) {
  
  #settings
  age <- max(abs(branching.times(sim_phylo)))
  esses <- 1
  quality_ratio <- esses/(chain_length/sample_interval)
  while (esses < 200 | quality_ratio < 0.2)
  {
    #run pirouette
    pirouette.out <- pirouette::pir_run( #test1 is a posterior of 10001 trees generated from sim_tes[[1]]
      phylogeny = sim_phylo,
      sequence_length = sequence_length,
      mutation_rate = 1,
      mcmc = create_mcmc(chain_length = chain_length, store_every = sample_interval), #store_every = -1
      site_models = beautier::create_jc69_site_model(),
      clock_models = beautier::create_strict_clock_model(),
      tree_priors = beautier::create_bd_tree_prior(),
      mrca_distr = create_normal_distr(mean = age, sigma = 0.001),
      alignment_rng_seed = 0,
      beast2_rng_seed = 1,
      verbose = FALSE,
      beast_jar_path = beastier::get_default_beast2_jar_path()
    )
    
    #check esses
    esses2 <- tracerer::calc_esses(pirouette.out$estimates, sample_interval = sample_interval) #effective sample size: has to be at least 200, try not to be over 1000 otherwise you're overkilling (try to stick to 1/5 rule of thumb)
    esses  <- min(esses2, na.rm = TRUE)
    quality_ratio <- esses/(chain_length/sample_interval)
    
    #this will apply at the next while run if you don't fulfill the quality condition
    chain_length    <- 2 * chain_length
    sample_interval <- 2 * sample_interval
  }
  
  #calculating nLTT statistics
  nLTT.diff <- rep(NA, length(pirouette.out$trees))
  for (i in 1:length(pirouette.out$trees))
  {
    nLTT.diff[i] <- nLTT::nLTTstat(sim_phylo, pirouette.out$trees[[i]])
  }
  mean.nLTT <- mean(nLTT.diff); mean.nLTT
  std.nLTT  <- sqrt(var(nLTT.diff)); std.nLTT
  df.nLTT <- data.frame(diff = nLTT.diff)
  
  return(list(posterior = pirouette.out, nLTT = df.nLTT))
}

comparison.MBD.dataset.with.BEAST.posterior <- function(sim_pars = c(0.2, 0.15, 2, 0.15),
                                                        cond = 1,
                                                        age = 10,
                                                        max_sims = 100,
                                                        tips_interval = c(0, Inf),) {
  
  #setting
  soc <- 2
  
  #sim dataset
  out <- sim_tes <- vector("list", max_sims)
  #sim_data <- sim_tas <- vector("list", max_sims); ext_species <- rep(NA, max_sims)
  for (s in 1:max_sims)
  {
    simulation <- MBD:::mbd_sim(pars = sim_pars,
                                soc = soc,
                                age = age,
                                cond = cond,
                                tips_interval = tips_interval)
    sim_tes[[s]]   <- simulation$tes
    # sim_tas[[s]]   <- simulation$tas
    # sim_data[[s]]  <- simulation$brts
    # ext_species[s] <- simulation$extinct_species
  } 
  
  #single tree analysis
  out[[s]] <- compare.MBD.with.alignments.BEAST.posterior(sim_phylo = sim_tes[[s]],
                                              chain_length = 1e+07,
                                              sample_interval = 1e+03,
                                              sequence_length = 1e+03)
  
  
  
  #show results
  nLTT.plot <- ggplot2::ggplot(df.nLTT, ggplot2::aes(x = diff)) + 
    ggplot2::geom_histogram(ggplot2::aes(y = ..count..), colour = "darkblue", fill = "lightblue", bins = 40) + 
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean(diff)), color = "red", linetype = "dashed", size = 1) +
    ggplot2::ggtitle(paste0("nLTT statistics for a MBD tree compared \nwith tree posterior generated from alignments and BD prior \nmean = ", signif(mean.nLTT, digits = 2), ", std = ", signif(std.nLTT, digits = 2))) +
    ggplot2::xlab("nLTT value")
  plot(nLTT.plot)
  
  #save results
  dir_name <- paste0(dirname(getwd()),"//results//nLTT//"); suppressWarnings(dir.create(dir_name))
  folder_name <- paste0(dir_name, toString(sim_pars))
  suppressWarnings(dir.create(folder_name))
  j <- 1
  while ((file_name <- paste0(folder_name, "//nLTT", j)) %in% list.files(folder_name)) {j <- j + 1}
  save(nLTT = df.nLTT, nLTT.plot = nLTT.plot, posterior = test1, file = file_name)
  
}

sim_pars = c(0.2, 0.15, 2, 0.15)

nLTT.MBD.test <- compare.MBD.with.alignments.BEAST.posterior(sim_pars = c(0.2, 0.15, 2, 0.15), sequence_length = 10^3)
nLTT.BD.test  <- compare.MBD.with.alignments.BEAST.posterior(sim_pars = c(0.2, 0.15, 0, 0.15), sequence_length = 10^3)