# @Giappo: add doc
#' Does something
#' @export
P_t1_t2 <- function(lambda, mu, t1, t2) {
  P_t1_t2 <- (1 - mu/lambda)/(1 - (mu/lambda) * exp(-(lambda - mu) *(t2 - t1)) )
  return(P_t1_t2)
}

# @Giappo: add doc
#' Does something
#' @export
BD.Nct <- function(lambda, mu, t, N0 = 2, age) {
  #average amount of species for a conditioned BD process
  #Nct <- N0 * Lt/(1 - ((mu * (Lt - 1))/(Lt * lambda - mu))^N0) #old version from etienne 2008
  tt  <- t
  TT  <- abs(age)
  Lt  <- exp((lambda - mu) * tt)
  Nct <- N0 * Lt * P_t1_t2(lambda = lambda, mu = mu, t1 = tt, t2 = TT) * 
    (P_t1_t2(lambda = lambda, mu = mu, t1 = 0, t2 = TT))^-1#new version from etienne & rosindell 2012, formula #15
  
  return(Nct)
}

# @Giappo: add doc
#' Does something
#' @export
BD.Nmutations <- function(lambda, mu, age, N0 = 2, sequence_length = 1000, mutation_rate = 1/age) {
  
  age <- abs(age)
  ft  <- function(t) {MBD::BD.Nct(t, lambda = lambda, mu = mu, N0 = N0, age = age)}
  Nmutations <- mutation_rate * sequence_length * integrate(f = ft, lower = 0, upper = age)[[1]]
  return(Nmutations)
}

# @Giappo: add doc
#' Does something
#' @export
BD.infer.lambda.from.mutations <- function(Nsubs, 
                                           MBD.lambda, 
                                           mu, 
                                           age, 
                                           N0 = 2, 
                                           sequence_length = 1000, 
                                           mutation_rate = 1/age, 
                                           Nsteps = 40) {
  lavec <- seq((min.lambda <- 0.5 * MBD.lambda), (max.lambda <- 6 * MBD.lambda), by = abs(max.lambda - min.lambda)/(Nsteps))
  NN <- rep(NA, length(lavec))
  for (i in 1:length(lavec))
  {
    NN[i] <- MBD::BD.Nmutations(lambda = lavec[i], mu = mu, age = age,
                           N0 = N0, sequence_length = sequence_length, mutation_rate = mutation_rate)
  }
  md <- lm(log(NN) ~ lavec)
  fit.test <- nls(NN ~ b*exp(m*lavec), start = list(m = coef(md)[2], b = coef(md)[1]))
  # plot(lavec, predict(fit.test), type = "l"); points(lavec, NN)
  # plot(coef(fit.test)[2] * exp(coef(fit.test)[1] * lavec) ~ lavec , type = "l"); points(lavec, NN)
  best.lambda <- (log(Nsubs) - log(coef(fit.test)[2]))/coef(fit.test)[1]; #best.lambda
  return(best.lambda)
}

# @Giappo: add doc
#' Does something
#' @export
alignments_comparison_single <- function(sim_phylo,
                                         chain_length = 1e+07,
                                         sample_interval = 1e+03,
                                         sequence_length = 1e+03,
                                         max_repetitions = 10,
                                         mutation_rate = 1/max(abs(ape::branching.times(sim_phylo)))) {
  
  #settings
  age <- max(abs(ape::branching.times(sim_phylo)))
  esses <- 1
  quality_ratio <- esses/(chain_length/sample_interval)
  repetitions <- 1
  while (esses < 200 | quality_ratio < 0.2)
  {
    #run pirouette
    pirouette.out <- pirouette::pir_run( #test1 is a posterior of 10001 trees generated from sim_tes[[1]]
      phylogeny = sim_phylo,
      sequence_length = sequence_length,
      mutation_rate = mutation_rate,
      mcmc = beautier::create_mcmc(chain_length = chain_length, store_every = sample_interval), #store_every = -1
      site_models = beautier::create_jc69_site_model(),
      clock_models = beautier::create_strict_clock_model(),
      tree_priors = beautier::create_bd_tree_prior(),
      mrca_distr = beautier::create_normal_distr(mean = beautier::create_mean_param(value = age), sigma = beautier::create_sigma_param(value = 0.001)),
      alignment_rng_seed = 0,
      beast2_rng_seed = 1,
      verbose = FALSE,
      beast_jar_path = beastier::get_default_beast2_jar_path()
    )
    
    #check esses
    esses2 <- tracerer::calc_esses(pirouette.out$estimates, sample_interval = sample_interval) #effective sample size: has to be at least 200, try not to be over 1000 otherwise you're overkilling (try to stick to 1/5 rule of thumb)
    esses  <- esses2$posterior
    quality_ratio <- esses/(nrow(pirouette.out$estimates))
    
    #this will apply at the next while run if you don't fulfill the quality condition
    chain_length    <- 2 * chain_length
    sample_interval <- 2 * sample_interval
    
    repetitions <- repetitions + 1
    if (repetitions >= max_repetitions) {break}
  }
  
  #calculating nLTT statistics
  nLTT.diff <- rep(NA, length(pirouette.out$trees))
  for (i in 1:length(pirouette.out$trees))
  {
    nLTT.diff[i] <- nLTT::nLTTstat(sim_phylo, pirouette.out$trees[[i]])
  }
  mean.nLTT <- mean(nLTT.diff); mean.nLTT
  std.nLTT  <- sqrt(var(nLTT.diff)); std.nLTT
  df.nLTT   <- data.frame(diff = nLTT.diff)
  
  return(list(alignment = pirouette.out$alignment, 
              trees = pirouette.out$trees,
              estimates = pirouette.out$estimates,
              nLTT = df.nLTT))
}

# @Giappo: add doc
#' Does something
#' @export
alignments_comparison_multiple <- function(sim_pars = c(0.2, 0.15, 2, 0.15),
                                           max_sims = 1e+01,
                                           chain_length = 1e+06,
                                           sample_interval = 1e+03,
                                           sequence_length = 1e+03,
                                           cond = 1,
                                           age = 10,
                                           tips_interval = c(0, Inf),
                                           mutation_rate = 1/age) {
  
  #setting
  soc <- 2 #it has to start with a crown to use pirouette
  
  #set filename to save results
  dir_name <- paste0(dirname(getwd()),"//results//nLTT//"); suppressWarnings(dir.create(dir_name))
  sim_pars_string <- gsub(", ", "-", toString(sim_pars))
  folder_name <- paste0(dir_name, sim_pars_string)
  suppressWarnings(dir.create(folder_name))
  j <- 1; while ((file_name <- paste0("nLTT", j)) %in% list.files(folder_name)) {j <- j + 1}
  
  #sim dataset
  MBD.estimates   <- MBD.trees <- MBD.alignment <- MBD.nLTT <- vector("list", max_sims)
  BD.simulations  <- BD.estimates  <- BD.trees  <- BD.alignment  <- BD.nLTT  <- vector("list", max_sims)
  Nsubstitutions  <- rep(NA, (max_sims2 <- 100 * max_sims))
  for (s in 1:max_sims)
  {
    MBD.simulation <- mbd_sim(pars = sim_pars,
                                    soc = soc,
                                    age = age,
                                    cond = cond,
                                    tips_interval = tips_interval)
    
    full_tree           <- MBD.simulation$tas#; plot(full_tree)
    total_branch_length <- sum(full_tree$edge.length) # total branch length
    Nsubstitutions[s]   <- sequence_length * mutation_rate * total_branch_length

    MBD.out <- alignments_comparison_single(sim_phylo = MBD.simulation$tes,
                                            chain_length = chain_length,
                                            sample_interval = sample_interval,
                                            sequence_length = sequence_length,
                                            mutation_rate = mutation_rate)
    
    MBD.nLTT[[s]]      <- MBD.out$nLTT
    MBD.alignment[[s]] <- MBD.out$alignment
    MBD.trees[[s]]     <- MBD.out$trees
    MBD.estimates[[s]] <- MBD.out$estimates
  # }
  
  # BD.lambda <- BD.infer.lambda.from.mutations(Nsubs = mean(Nsubstitutions),
    BD.lambda <- BD.infer.lambda.from.mutations(Nsubs = Nsubstitutions[s],
                                                MBD.lambda = sim_pars[1],
                                                mu = sim_pars[2],
                                                age = age,
                                                N0 = soc,
                                                sequence_length = sequence_length,
                                                mutation_rate = mutation_rate,
                                                Nsteps = 40)
  
  # for (s in 1:max_sims)
  # {
    BD.simulations[[s]] <- mbd_sim(pars = c(BD.lambda, sim_pars[2], 0, 0),
                                         soc = soc,
                                         age = age,
                                         cond = cond,
                                         tips_interval = tips_interval)
    
    BD.out <- alignments_comparison_single(sim_phylo = BD.simulations[[s]]$tes,
                                           chain_length = chain_length,
                                           sample_interval = sample_interval,
                                           sequence_length = sequence_length)
    
    BD.nLTT[[s]]      <- BD.out$nLTT
    BD.alignment[[s]] <- BD.out$alignment
    BD.trees[[s]]     <- BD.out$trees
    BD.estimates[[s]] <- BD.out$estimates 
  }
  
  #roba
  MBD.nLTT.total <- unname(unlist(MBD.nLTT))
  BD.nLTT.total  <- unname(unlist(BD.nLTT))
  df.MBD.nLTT <- data.frame(MBD.nLTT = MBD.nLTT.total)
  df.BD.nLTT  <- data.frame(BD.nLTT  = BD.nLTT.total)
  
  MBD.mean <- mean(df.MBD.nLTT$MBD.nLTT)
  MBD.std  <- sqrt(var(df.MBD.nLTT$MBD.nLTT))
  BD.mean  <- mean(df.BD.nLTT$BD.nLTT)
  BD.std   <- sqrt(var(df.BD.nLTT$BD.nLTT))
  
  #show results
  MBD.nLTT.plot <- ggplot2::ggplot(df.MBD.nLTT, ggplot2::aes(x = MBD.nLTT)) + 
    ggplot2::geom_histogram(ggplot2::aes(y = ..count..), colour = "darkblue", fill = "lightblue", bins = 40) + 
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean(MBD.nLTT)), color = "red", linetype = "dashed", size = 1) +
    ggplot2::ggtitle(paste0("nLTT statistics for a MBD tree compared \nwith tree posterior generated from alignments and BD prior \nmean = ", signif(MBD.mean, digits = 2), ", std = ", signif(MBD.std, digits = 2))) +
    ggplot2::xlab("nLTT value")
  plot(MBD.nLTT.plot)
  
  BD.nLTT.plot <- ggplot2::ggplot(df.BD.nLTT, ggplot2::aes(x = BD.nLTT)) + 
    ggplot2::geom_histogram(ggplot2::aes(y = ..count..), colour = "darkblue", fill = "lightblue", bins = 40) + 
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean(BD.nLTT)), color = "red", linetype = "dashed", size = 1) +
    ggplot2::ggtitle(paste0("nLTT statistics for a BD tree compared \nwith tree posterior generated from alignments and BD prior \nmean = ", signif(BD.mean, digits = 2), ", std = ", signif(BD.std, digits = 2))) +
    ggplot2::xlab("nLTT value")
  plot(BD.nLTT.plot)
  
  #save results
  save(MBD.nLTT = df.MBD.nLTT,
       MBD.alignment = MBD.alignment,
       MBD.trees = MBD.trees,
       MBD.estimates = MBD.estimates,
       BD.nLTT = df.BD.nLTT,
       BD.alignment = BD.alignment,
       BD.trees = BD.trees,
       BD.estimates = BD.estimates,
       file = paste0(folder_name,"//" ,file_name))
  
  png(filename = paste0(folder_name, "//" ,file_name, "_MBD_plot.png"))
  plot(MBD.nLTT.plot)
  dev.off()
        
  png(filename = paste0(folder_name, "//" ,file_name, "_BD_plot.png"))
  plot(BD.nLTT.plot)
  dev.off()
  
  return(list(MBD.nLTT = df.MBD.nLTT, BD.nLTT = df.BD.nLTT))
}


#test <- MBD::alignments_comparison_multiple(sim_pars = (sim_pars <- c(0.2, 0.15, 2, 0.15)), max_sims = 100)

##### DUMPSTER
# devtools::install_github("richelbilderbeek/tracerer", dependencies = TRUE)
# devtools::install_github("richelbilderbeek/beastier", dependencies = TRUE)
# devtools::install_github("richelbilderbeek/babette", dependencies = TRUE)
# devtools::install_github("richelbilderbeek/pirouette", dependencies = TRUE)
# devtools::install_github("richelbilderbeek/ribir", dependencies = TRUE)
# install.packages('xmlparsedata')
# install.packages('lintr')
# devtools::install_github("richelbilderbeek/raket", dependencies = TRUE)
# install.packages("nLTT")
# library("ape"); library("raket"); library("nLTT")