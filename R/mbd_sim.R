# mbd_sim---------------------------------
#' @author Giovanni Laudanno
#' @title Creates simulated trees under the multiple birth death process, including both sympatric and allopatric speciation
#' @description mbd_sim produces simulated trees allowing for three kind of events: sympatric speciation, multiple allopatric speciations and extinction.
#' @inheritParams default_params_doc
#' @param pars vector of parameters:
#' \itemize{
#'   \item pars[1] is lambda, the sympatric speciation rate;
#'   \item pars[2] is mu, the extinction rate;
#'   \item pars[3] is nu, the multiple allopatric speciation trigger rate;
#'   \item pars[4] is q, the single-lineage speciation probability.
#' }
#' @param soc Sets whether stem or crown age should be used (1 or 2).
#' @param age The age of the tree.
#' @param cond Set 1 if you want to condition on stem or crown age and non-extinction of the phylogeny. Set 0 otherwise.
#' @param tips_interval Sets tips boundaries constrain on simulated dataset.
#' @return The function returns a list of L table, branching times and number of extinct species.
#' \itemize{
#'   \item brts are the branching times;
#'   \item tes is a tree-object corresponding to the reconstructed phylogenetic tree;
#'   \item tas is a tree-object corresponding to the full phylogenetic tree;
#'   \item extinct_species is the number of species gone extinct before the present time;
#'   \item L is a matrix of all species where the columns are:
#'   \itemize{
#'      \item first is the time at which a species is born;
#'      \item second is the label of the parent of the species; positive and negative values only indicate whether the species belongs to the left or right crown lineage;
#'      \item third is the label of the daughter species itself; positive and negative values only indicate whether the species belongs to the left or right crown lineage;
#'      \item fourth is the time of extinction of the species. If this is equal to -1, then the species is still extant.
#'   }
#' }
#'
#' @examples
#' out = mbd_sim( pars=c(0.6,0.1,0.4,0.1),soc=2,age=10,cond=1,tips_interval=c(0,Inf) )
#' plot(out$tas)
#' plot(out$tes)
#' out$L
#'
#' @export
mbd_sim <- function(pars, soc = 2, age = 10, cond = 1,
                    tips_interval = c(0, Inf), minimum_multiple_births = 0)
{
  if (tips_interval[2] < tips_interval[1] || any(pars < 0)){stop("ERROR! Check again your settings.")}
  lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4]; N0 <- soc
  tips <- -1; crown_species_dead <- cond; multiple_births_check <- 0;
  keep_the_sim <- 0
  while (keep_the_sim == 0 | multiple_births_check == 0)
  {
    total_count <- N0
    pool <- 1:N0
    while (total_count == N0 | length(pool) < N0)
    {
      total_count <- N0
      N <- N0
      pool <- c(-1,2)
      t <- age
      L <- matrix(0, nrow = 1e6, 4)
      L[,4] <- -1
      L[,3] <- 0
      L[1,1:4] <- c(t,0,-1,-1)
      L[2,1:4] <- c(t,-1,2,-1)
      while (t > 0)
      {
        N <- length(pool)
        total_rate <- N * (lambda + mu) + nu
        if (total_rate > 0)
        {
          deltaT <- rexp(1, rate = total_rate)
          outcome <- sample(c(-1,1,2), size = 1, prob = c(N * mu, N * lambda, nu))
          deltaN <- -1 * (outcome == -1) + 1 * (outcome == 1) + (outcome == 2) * rbinom(n = 1, size = N, prob = q)
          t <- t - deltaT
  
          if (deltaN > 0 & t > 0)
          {
            if (N > 1) {parents <- sample(pool, replace = FALSE, size = deltaN)}else{parents <- pool}
            new_interval <- (total_count + 1):(total_count + deltaN)
            L[new_interval,1] <- t#-(deltaN:1)*1e-5 add this if you need separate time points
            L[new_interval,2] <- parents
            L[new_interval,3] <- abs(new_interval) * sign(parents)
  
            pool <- c(pool, abs(new_interval) * sign(parents) )
            total_count <- total_count + deltaN
          }
          if (deltaN < 0 & t > 0)
          {
            if (N > 1) {dead <- sample(pool, replace = FALSE, size = 1)}else{dead <- pool}
            L[abs(dead),4] <- t
            pool <- pool[pool != dead]
          }
          # print(pool)
        }else
        {
          t <- 0 
        }
      }
    }
    L <- L[(1:total_count),]
    extinct_species <- sum(L[,4] != -1)
    #tips check
    tips <- length(L[,4][L[,4] == -1])
    #survival of crown check
    alive <- L[L[,4] == -1,]
    alive <- matrix(alive, ncol = 4)
    crown_species_dead <- ( length( unique(sign(alive[,3])) ) != 2 ) * cond #if cond == 0 they will always look like they're alive, because I don't care
    #multiple births check
    births.reconstructed_tree <- unlist(unname(sort(DDD::L2brts(L, dropextinct = TRUE), decreasing = TRUE)) )
    births.full_tree <- unlist(unname(sort(DDD::L2brts(L, dropextinct = FALSE), decreasing = TRUE)) )
    multiple_births.reconstructed_tree <- sum( duplicated(births.reconstructed_tree) )
    multiple_births.full_tree <- sum( duplicated(births.full_tree) )
    #should i consider the full tree or the reconstructed one???
    # multiple_births_check = (multiple_births.full_tree>=minimum_multiple_births)
    multiple_births_check <- (multiple_births.reconstructed_tree >= minimum_multiple_births)
    
    keep_the_sim <- (!crown_species_dead) & (tips >= tips_interval[1] & tips <= tips_interval[2]) #should i keep this simulation?
  }
  time_points <- unlist(unname(sort(DDD::L2brts(L, dropextinct = TRUE), decreasing = TRUE)) )
  brts <- -sort(abs(as.numeric(time_points)), decreasing = TRUE)
  tes <- DDD::L2phylo(L, dropextinct = TRUE)
  tas <- DDD::L2phylo(L, dropextinct = FALSE)
  #   plot(tas)
  #   plot(tes)
  out <- list(brts = brts, tes = tes, tas = tas, extinct_species = extinct_species, L = L, 
              minimum_multiple_births = multiple_births.full_tree)
  return(out)
}

# mbd_sim0---------------------------------
#' @author Giovanni Laudanno
#' @title Creates simulated trees under the multiple birth death process
#' @description mbd_sim0 produces simulated trees including multiple speciations at the same time.
#' @param pars vector of parameters:
#' \itemize{
#'   \item pars[1] is the multiple speciation trigger rate;
#'   \item pars[2] is the extinction rate;
#'   \item pars[3] is the single-lineage speciation probability.
#' }
#' @param soc Sets whether stem or crown age should be used (1 or 2).
#' @param age The age of the tree.
#' @param cond Set 1 if you want to condition on stem or crown age and non-extinction of the phylogeny. Set 0 otherwise.
#' @param tips_interval Sets tips boundaries constrain on simulated dataset.
#' @return The function returns a list of L table, branching times and number of extinct species.
#' \itemize{
#'   \item brts are the branching times;
#'   \item tes is a tree-object corresponding to the reconstructed phylogenetic tree;
#'   \item tas is a tree-object corresponding to the full phylogenetic tree;
#'   \item extinct_species is the number of species gone extinct before the present time;
#'   \item L is a matrix of all species where the columns are:
#'   \itemize{
#'      \item first is the time at which a species is born;
#'      \item second is the label of the parent of the species; positive and negative values only indicate whether the species belongs to the left or right crown lineage;
#'      \item third is the label of the daughter species itself; positive and negative values only indicate whether the species belongs to the left or right crown lineage;
#'      \item fourth is the time of extinction of the species. If this is equal to -1, then the species is still extant.
#'   }
#' }
#'
#' @examples
#' out = mbd_sim0( pars=c(2.5,0.1,0.1),soc=2,age=10,cond=1,tips_interval=c(0,Inf) )
#' plot(out$tas)
#' plot(out$tes)
#' out$L
#'
#' @export
mbd_sim0 <- function(pars, soc = 2, age = 10, cond = 1,
                     tips_interval = c(soc * (cond == 1), Inf), minimum_multiple_births = 0)
{
  if (tips_interval[2]<tips_interval[1]){stop("ERROR! Check again your settings.")}
  lambda=pars[1]; mu=pars[2]; q=pars[3]; N0=soc
  tips=-1; conditioning_on_survival=cond; multiple_births_check = 0;
  while ( tips<tips_interval[1] | tips>tips_interval[2] | conditioning_on_survival | multiple_births_check == 0)
  {
    total_count=N0
    pool=1:N0
    while (total_count==N0 | length(pool)<N0)
    {
      total_count=N0
      N=N0
      pool=c(-1,2)
      t=age
      L = matrix(0,nrow=1e6,4)
      L[,4]=-1
      L[,3]=0
      L[1,1:4] = c(t,0,-1,-1)
      L[2,1:4] = c(t,-1,2,-1)
      while (t>0)
      {
        N=length(pool)
        deltaT=rexp(1,rate=(lambda+N*mu))
        outcome=sample(c(-1,1),size=1,prob = c(N*mu,lambda))
        deltaN=rbinom(n=1,size=N,prob=q)*(outcome==1)-1*(outcome==-1)
        t=t-deltaT
        if (deltaN>0 & t>0)
        {
          if (N>1) {parents=sample(pool,replace = F,size=deltaN)}else{parents=pool}
          new_interval=(total_count+1):(total_count+deltaN)
          L[new_interval,1]=t#-(deltaN:1)*1e-5 add this if you need separate time points
          L[new_interval,2]=parents
          L[new_interval,3]=abs(new_interval)*sign(parents)

          pool=c(pool, abs(new_interval)*sign(parents) )
          total_count=total_count+deltaN
        }
        if (deltaN<0 & t>0)
        {
          if (N>1) {dead=sample(pool,replace = F,size=1)}else{dead=pool}
          # dead=abs(dead)
          L[abs(dead),4]=t
          pool=pool[pool!=dead]
        }
        # print(pool)
      }
    }
    L=L[(1:total_count),]
    extinct_species = sum(L[,4]!=-1)
    #tips check
    tips=length(L[,4][L[,4]==-1])
    #survival of crown check
    alive=L[L[,4]==-1,]
    alive=matrix(alive,ncol=4)
    conditioning_on_survival = ( length( unique(sign(alive[,3])) )!=2 )*cond
    #multiple births check
    #multiple births check
    births.reconstructed_tree=unlist(unname(sort(DDD::L2brts(L,dropextinct = T),decreasing = T)) )
    births.full_tree=unlist(unname(sort(DDD::L2brts(L,dropextinct = F),decreasing = T)) )
    multiple_births.reconstructed_tree = sum( duplicated(births.reconstructed_tree) )
    multiple_births.full_tree = sum( duplicated(births.full_tree) )
    #should i consider the full tree or the reconstructed one???
    # multiple_births_check = (multiple_births.full_tree>=minimum_multiple_births)
    multiple_births_check = (multiple_births.reconstructed_tree>=minimum_multiple_births)
  }
  time_points=unlist(unname(sort(DDD::L2brts(L,dropextinct = T),decreasing = T)) )
  brts = -sort(abs(as.numeric(time_points)),decreasing = TRUE)
  tes = DDD::L2phylo(L,dropextinct = T)
  tas = DDD::L2phylo(L,dropextinct = F)
  #   plot(tas)
  #   plot(tes)
  out = list(brts=brts,tes = tes, tas = tas, extinct_species=extinct_species,L = L, minimum_multiple_births = minimum_multiple_births)
  return(out)
}

# mbd_sim_dataset---------------------------------
#' @author Giovanni Laudanno
#' @title Creates a full simulated dataset of trees under the multiple birth death process
#' @description mbd_sim_dataset produces a full dataset of max_sims simulated trees including multiple speciations at the same time. This second version takes into account the possibility both of allopatric and sympatric speciation.
#' @param pars vector of parameters:
#' \itemize{
#' \item id == 1 corresponds to lambda (speciation rate)
#' \item id == 2 corresponds to mu (extinction rate)
#' \item id == 3 corresponds to nu (multiple speciation trigger rate)
#' \item id == 4 corresponds to q (single-lineage speciation probability)
#' }
#' @param soc stands for stem or crown. Set 1 for stem or 2 for crown.
#' @param age is the age of the tree.
#' @param cond conditions on the survival of the crown (or the stem). Set to 1 if you want to condition on crown survival. Else set 0.
#' @param edge the program automatically detects the (estimated) average number of tips. "edge" defines the width of the spread around the mean.
#' @param tips_interval You can also define the tips_interval as you can usually do with a standard usage of mbd_sim.
#' @return The function returns a list of brts vectors, one for each of the max_sims simulations.
#' It also saves those in a file called "sim_data" and saves all the settings in a "general_settings" file.
#' N.B.: At each call of the function you overwrite the previous files.
#'
#' @examples
#' # @Giappo: does not work
#' # out <- mbd_sim_dataset( pars=c(0.4,0.1,0.2,0.15),soc=2,age=10,cond=1,edge=Inf )
#'
#' @export
mbd_sim_dataset <- function(sim_pars = c(0.5,0.1,0.3,0.15), soc = 2, cond = 1, age = 10, max_sims = 1000,
                           tips_interval = c(0, 100), edge = Inf, minimum_multiple_births = 0){
  # mbd_sim_dataset creates a full simulated dataset of "max_sims" trees
  #"edge" gives the extent of fluctuations around the mean that i want to consider

  if (sim_pars[2] == 0){cond = 0; tips_interval = c(0, Inf)} #this allows me to use the analytical formula CHECK THIS!

  N0 = soc
  # if (cond == 1){tips_interval[1] = max(N0, tips_interval[1])}#if the tree is conditioned on the survival of crown species the minimum amount of tips has to be raised!!!
  if (edge != Inf && tips_interval == c(0, Inf))
  {
    estimation <- mbd_P_eq(test_parameters = sim_pars, age = age, max_number_of_species = 3000, precision = 50L,soc=soc,output=0);
    max_tips <- round(estimation$avg_n*(1 + edge));
    min_tips <- max(0,round(estimation$avg_n * (1 - edge)));
    tips_interval <- c(min_tips, max_tips) #c(10,45) #setting the upper limit over 50 may be a problem #min and max number of tips for simulated trees
  }

  #simulate trees
  sim_data <- sim_tes <- sim_tas <- vector("list", max_sims)
  ext_species <- rep(NA, max_sims)
  for (s in 1:max_sims)
  {
    set.seed(s)
    simulation <- mbd_sim(pars = sim_pars, soc = soc, age = age, cond = cond,
                                tips_interval = tips_interval, minimum_multiple_births = minimum_multiple_births)
    sim_data[[s]]  <- simulation$brts
    ext_species[s] <- simulation$extinct_species
    sim_tes[[s]]   <- simulation$tes
    sim_tas[[s]]   <- simulation$tas
  }

  max_k <- (N0 - 1) + (is.list(sim_data)) * max(sapply(sim_data, length)) + (1 - is.list(sim_data)) * length(sim_data)
  if (is.list(sim_data))
  {
    all_the_births=sapply(sim_data, FUN = function(brts){return(brts2time_intervals_and_births(brts)$births)})
  }else{
    all_the_births=brts2time_intervals_and_births(sim_data)$births
  }
  max_b <- max(unlist(all_the_births))

  additional_species <- tips <- rep(0, max_sims)
  for (s in 1:max_sims)
  {
    additional_species[s] <- sum(duplicated(sim_data[[s]]))
    tips[s] <- length(sim_data[[s]]) + 1
  }
  #saving sims and settings
  # simpath=paste("sims/",sim_pars[1],"-",sim_pars[2],"-",sim_pars[3],"-",sim_pars[4],"/",sep = '')
  simpath <- getwd()
  datapath <- paste0(simpath,"/data")
  sim_data_name <- paste0(datapath,"/sim_data")
  general_settings_name <- paste0(datapath,"/general_settings")
  sim_trees_name <- paste0(datapath,"/sim_trees")
  if (file.exists(sim_data_name)){suppressWarnings( file.remove(sim_data_name) )}
  save(sim_data, file = sim_data_name)
  if (file.exists(general_settings_name)){suppressWarnings( file.remove(general_settings_name) )}
  save(sim_pars,soc,age,cond,max_sims,tips_interval,max_k,max_b,ext_species,additional_species,tips,minimum_multiple_births,file=general_settings_name)
  if (file.exists(sim_trees_name)){suppressWarnings( file.remove(sim_trees_name) )}
  save(sim_tas,sim_tes, file = sim_trees_name)
  return(sim_data)
}

# mbd_sim_dataset0---------------------------------
#' @author Giovanni Laudanno
#' @title Creates a full simulated dataset of trees under the multiple birth death process
#' @description mbd_sim_dataset0 produces a full dataset of max_sims simulated trees including multiple speciations at the same time.
#' @param pars vector of parameters:
#' \itemize{
#'   \item pars[1] is the multiple speciation trigger rate;
#'   \item pars[2] is the extinction rate;
#'   \item pars[3] is the single-lineage speciation probability.
#' }
#' @param soc stands for stem or crown. Set 1 for stem or 2 for crown.
#' @param age is the age of the tree.
#' @param cond conditions on the survival of the crown (or the stem). Set to 1 if you want to condition on crown survival. Else set 0.
#' @param edge the program automatically detects the (estimated) average number of tips. "edge" defines the width of the spread around the mean.
#' @param tips_interval You can also define the tips_interval as you can usually do with a standard usage of mbd_sim0.
#' @return The function returns a list of brts vectors, one for each of the max_sims simulations.
#' It also saves those in a file called "sim_data" and saves all the settings in a "general_settings" file.
#' N.B.: At each call of the function you overwrite the previous files.
#'
#' @examples
#' # @Giappo: does not work
#' # out <- mbd_sim_dataset0( pars=c(2.5,0.1,0.1),soc=2,age=10,cond=1,edge=Inf )
#'
#' @export
mbd_sim_dataset0 <- function(sim_pars = c(2.5,0.1,0.10), soc = 2, cond = 1, age = 10,
                             max_sims = 1000, tips_interval = c(0,Inf),edge = Inf,
                             minimum_multiple_births = 0){
  # mbd_sim_dataset0 creates a full simulated dataset of "max_sims" trees, using only three parameters.
  #"edge" gives the extent of fluctuations around the mean that i want to consider

  N0 <- soc
  if (cond == 1){tips_interval[1]=max(N0,tips_interval[1])}#if the tree is conditioned on the survival of crown species the minimum amount of tips has to be raised!!!
  if (edge!=Inf && tips_interval == c(0,Inf) ){
    estimation=mbd_P_eq(test_parameters=sim_pars,age=age,max_number_of_species = 3000, precision = 50L,soc=soc,output=0);
    max_tips=round(estimation$avg_n*(1+edge));
    min_tips=max(0,round(estimation$avg_n*(1-edge)));
    tips_interval=c(min_tips,max_tips) #c(10,45) #setting the upper limit over 50 may be a problem #min and max number of tips for simulated trees
  }

  #simulate trees
  sim_data=sim_tes=sim_tas=vector("list",max_sims)
  ext_species=rep(NA,max_sims)
  for (s in 1:max_sims){
    simulation <- mbd_sim0(pars = sim_pars, soc = soc, age = age, cond = cond,
                              tips_interval = tips_interval, minimum_multiple_births = minimum_multiple_births)
    sim_data[[s]]=simulation$brts
    ext_species[s]=simulation$extinct_species
    sim_tes[[s]] = simulation$tes
    sim_tas[[s]] = simulation$tas
  }

  max_k=(N0-1)+(is.list(sim_data))*max(sapply(sim_data, length))+(1-is.list(sim_data))*length(sim_data)
  if (is.list(sim_data))
  {
    all_the_births=sapply(sim_data, FUN = function(brts){return(brts2time_intervals_and_births(brts)$births)})
  }else{
    all_the_births=brts2time_intervals_and_births(sim_data)$births
  }
  max_b=max(unlist(all_the_births))

  additional_species = tips = rep(0,max_sims)
  for (s in 1:max_sims){
  additional_species[s] = sum( duplicated(sim_data[[s]]) )
  tips[s] = length(sim_data[[s]])+1
  }
  #saving sims and settings
  # simpath=paste("sims/",sim_pars[1],"-",sim_pars[2],"-",sim_pars[3],"/",sep = '')
  simpath = getwd()
  datapath=paste(simpath,"/data",sep = '')
  sim_data_name = paste(datapath,"/sim_data",sep = '')
  general_settings_name = paste(datapath,"/general_settings",sep = '')
  sim_trees_name = paste(datapath,"/sim_trees",sep = '')
  if (file.exists(sim_data_name)){suppressWarnings( file.remove(sim_data_name) )}
  save(sim_data,file=sim_data_name)
  if (file.exists(general_settings_name)){suppressWarnings( file.remove(general_settings_name) )}
  save(sim_pars,soc,age,cond,max_sims,tips_interval,max_k,max_b,ext_species,additional_species,tips,file=general_settings_name)
  if (file.exists(sim_trees_name)){suppressWarnings( file.remove(sim_trees_name) )}
  save(sim_tas,sim_tes,file=sim_trees_name)
  return(sim_data)
}

# mbd_sim3 - experimental: DON'T USE------------------------------
#' #' @author Giovanni Laudanno
#' #' @title Creates simulated trees under the multiple birth death process, including both sympatric and allopatric speciation
#' #' @description mbd_sim3 produces simulated trees in the same way as mbd_sim. The only difference is that you can only allow multiple births event to occur.
#' #' @export
#' mbd_sim3=function(pars,soc=2,age=10,cond=1,tips_interval = c(soc*(cond==1),Inf))
#' {
#'   if (tips_interval[2]<tips_interval[1]){print("ERROR! Check again your settings.");break}
#'   lambda = pars[1]; mu = pars[2]; nu = pars[3]; q=pars[4]; N0=soc
#'   if (lambda!=0){print("This is supposed to work only with multiple events. Lambda is not allowed.");break}
#'   tips=-1; conditioning_on_survival=cond;
#'   while ( tips<tips_interval[1] | tips>tips_interval[2] | conditioning_on_survival )
#'   {
#'     total_count=N0
#'     pool=1:N0
#'     while (total_count==N0 | length(pool)<N0 )
#'     {
#'       total_count=N0
#'       N=N0
#'       pool=c(-1,2)
#'       t=age
#'       L = matrix(0,nrow=1e6,4)
#'       L[,4]=-1
#'       L[,3]=0
#'       L[1,1:4] = c(t,0,-1,-1)
#'       L[2,1:4] = c(t,-1,2,-1)
#'       while (t>0)
#'       {
#'         N=length(pool)
#'         total_rate = N*(lambda + mu) + nu
#'         deltaT=rexp(1,rate = total_rate)
#'         outcome=sample(c(-1,1,2),size=1,prob = c(N*mu,N*lambda,nu))
#'         deltaN=-1*(outcome==-1)+1*(outcome==1)+rbinom(n=1,size=N,prob=q)*(outcome==2)
#'         t=t-deltaT
#'
#'         if (deltaN>1 & t>0)
#'         {
#'           if (N>1) {parents=sample(pool,replace = F,size=deltaN)}else{parents=pool}
#'           new_interval=(total_count+1):(total_count+deltaN)
#'           L[new_interval,1]=t#-(deltaN:1)*1e-5 add this if you need separate time points
#'           L[new_interval,2]=parents
#'           L[new_interval,3]=abs(new_interval)*sign(parents)
#'
#'           pool=c(pool, abs(new_interval)*sign(parents) )
#'           total_count=total_count+deltaN
#'         }
#'         if (deltaN<0 & t>0)
#'         {
#'           if (N>1) {dead=sample(pool,replace = F,size=1)}else{dead=pool}
#'           # dead=abs(dead)
#'           L[abs(dead),4]=t
#'           pool=pool[pool!=dead]
#'         }
#'         # print(pool)
#'       }
#'     }
#'     L=L[(1:total_count),]
#'     extinct_species = sum(L[,4]!=-1)
#'     #tips check
#'     tips=length(L[,4][L[,4]==-1])
#'     #survival of crown check
#'     alive=L[L[,4]==-1,]
#'     alive=matrix(alive,ncol=4)
#'     conditioning_on_survival = ( length( unique(sign(alive[,3])) )!=2 )*cond
#'   }
#'   time_points=unlist(unname(sort(DDD::L2brts(L,dropextinct = T),decreasing = T)) )
#'   brts = -sort(abs(as.numeric(time_points)),decreasing = TRUE)
#'   tes = DDD::L2phylo(L,dropextinct = T)
#'   tas = DDD::L2phylo(L,dropextinct = F)
#'   #   plot(tas)
#'   #   plot(tes)
#'   out = list(brts=brts,tes = tes, tas = tas, extinct_species=extinct_species,L = L)
#'   return(out)
#' }
