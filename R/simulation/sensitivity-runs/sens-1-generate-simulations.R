#### SETUP ---------------------------------------------------------------------
library(finalsize)
library(dplyr)
library(reshape2)
library(ggplot2)
library(deSolve)

source("./R/simulation/0-final-size-functions.R")
source("./R/simulation/sensitivity-runs/sens-0-parameters.R")


#### GENERATE SIMULATIONS ------------------------------------------------------
# run simulation
sim_out <- vector("list", n_reps)
for(i in 1:n_reps){
  print(i)
  sim_out[[i]] <- full_sim(n_locations = n_loc, n_models = n_models,
                      vax_cov_S1 = vax_scenarios[1], vax_cov_S2 = vax_scenarios[2],
                      lengthout_scenarios = 2, # only run for S1, S2, and T (don't need intermediate values here)
                      R0_lwr = 2, R0_upr = 3, model_bias_ind_sd = 0.05,
                      seed = seed_id[i], final_size_method = "simulation", 
                      alpha_lwr = 0.95, alpha_upr = 1, alpha_sd = 0.01)
}

saveRDS(sim_out, "output/simulation/sim_out_sens.rds")

