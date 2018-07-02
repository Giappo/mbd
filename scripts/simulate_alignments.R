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

#settings
sim_pars = c(0.2, 0.15, 2, 0.15)
soc = 2
cond = 1
age = 10
max_sims = 100
tips_interval <- c(0, Inf)
sequence_length <- 10^3

#sim trees
sim_data <- sim_tes <- sim_tas <- vector("list",max_sims)
ext_species <- rep(NA, max_sims)
for (s in 1:max_sims)
{
  simulation <- MBD:::mbd_sim(pars = sim_pars, soc = soc, age = age, cond = cond,
                              tips_interval = tips_interval)
  sim_data[[s]]  <- simulation$brts
  ext_species[s] <- simulation$extinct_species
  sim_tes[[s]]   <- simulation$tes
  sim_tas[[s]]   <- simulation$tas
}


# alignment <- pirouette::sim_alignment(
#   phylogeny = sim_tes[[1]],
#   sequence_length = 10^3,
#   mutation_rate = 1
# )

#run pirouette
test1 <- pirouette::pir_run( #test1 is a posterior of 10001 trees generated from sim_tes[[1]]
  phylogeny = sim_tes[[1]],
  sequence_length = sequence_length,
  mutation_rate = 1,
  mcmc = create_mcmc(chain_length = 1e+07, store_every = -1),
  site_models = beautier::create_jc69_site_model(),
  clock_models = beautier::create_strict_clock_model(),
  tree_priors = beautier::create_bd_tree_prior(),
  crown_age = age,
  mrca_distr = NA,
  alignment_rng_seed = 0,
  beast2_rng_seed = 1,
  verbose = FALSE,
  beast_jar_path = beastier::get_default_beast2_jar_path()
)

#analyze differences
plot(sim_tes[[1]])
test.brts <- sort(ape::branching.times(sim_tes[[1]]), decreasing = T); MBD:::brts2time_intervals_and_births(test.brts)$births
plot(test1$trees[[1]])

nLTT.diff <- rep(NA, length(test1$trees))
for (i in 1:length(test1$trees))
{
  nLTT.diff[i] <- nLTT::nLTTstat(sim_tes[[1]], test1$trees[[i]])
}
mean.nLTT <- mean(nLTT.diff); mean.nLTT
std.nLTT  <- sqrt(var(nLTT.diff)); std.nLTT

#show results
df.nLTT <- data.frame(diff = nLTT.diff)
ggplot2::ggplot(df.nLTT, ggplot2::aes(x = diff)) + 
  ggplot2::geom_histogram(ggplot2::aes(y = ..count..), colour = "darkblue", fill = "lightblue", bins = 40) + 
  ggplot2::geom_vline(ggplot2::aes(xintercept = mean(diff)), color = "red", linetype = "dashed", size = 1) +
  ggplot2::ggtitle(paste0("nLTT statistics for a MBD tree compared \n
                   with tree posterior generated from alignments and BD prior \n
                   mean = ",mean.nLTT, ", std = ", std.nLTT)) +
  ggplot2::xlab("nLTT value")
