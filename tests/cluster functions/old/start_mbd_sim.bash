#!/bin/bash

echo "sbatch install_packages.bash"

echo "library(MBD)" > zzz_sim.R
echo "c_sim_pars=c(2.5,0.1,0.1)" >> zzz_sim.R
echo "c_soc=2" >> zzz_sim.R
echo "c_cond=1" >> zzz_sim.R
echo "c_age=10" >> zzz_sim.R
echo "c_max_sims=1000" >> zzz_sim.R
echo "c_multiple_births_interval=c(1,Inf)" >> zzz_sim.R
echo "c_tips_interval=c(0,50)" >> zzz_sim.R
echo "MBD:::mbd_sim_dataset(sim_pars=c_sim_pars,soc=c_soc,cond=c_cond,age=c_age,max_sims=c_max_sims,tips_interval = c_tips_interval,multiple_births_interval=c_multiple_births_interval)" >> zzz_sim.R

echo "module load R/3.3.1-foss-2016a"
echo "Rscript zzz_sim.R" 

