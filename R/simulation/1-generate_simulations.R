#### SETUP ---------------------------------------------------------------------
library(finalsize)
library(dplyr)
library(reshape2)
library(ggplot2)
library(deSolve)

source("./R/simulation/0-final-size-functions.R")
source("./R/simulation/0-parameters.R")


#### GENERATE SIMULATIONS ------------------------------------------------------
reps = 1
set.seed(10)
seed_id = sample(1:1E6, reps)

# run simulation
sim_out <- full_sim(n_locations = n_loc, n_models = n_models,
                    vax_cov_S1 = vax_scenarios[1], vax_cov_S2 = vax_scenarios[2],
                    R0_lwr = 2, R0_upr = 3, model_bias_ind_sd = 0.05,
                    seed = seed_id, final_size_method = "simulation", 
                    alpha_lwr = 0.95, alpha_upr = 1, alpha_sd = 0.01)

saveRDS(sim_out, "output/simulation/sim_out.rds")

#### A FEW PLOTS ---------------------------------------------------------------
left_join(sim_out$model_sims, 
          sim_out$true_sims) %>%
  mutate(final_size_plot = ifelse(model_id == "T", true_final_size, final_size)) %>%
  ggplot() + 
  geom_line(aes(x = vax_cov, y = final_size_plot, color = model_id)) +
  geom_line(data = sim_out$true_sims, aes(x = vax_cov, y = true_final_size), color = "black", size = 1) + 
  facet_wrap(vars(paste0(round(location_R0,2), " (", location_id, ")")), scales = "free") + 
  theme_bw() + 
  theme(legend.position = "none")

ggplot(data = sim_out$errors %>% filter(scenario_id == "T") %>% filter(scenario_id == "T"), aes(x = vax_cov, y = error)) +
  geom_hline(yintercept = 0) +
  geom_line(data = sim_out$errors, aes(group = location_id), alpha = 0.1) +
  geom_point() +
  facet_wrap(vars(model_id)) +
  theme_bw()

