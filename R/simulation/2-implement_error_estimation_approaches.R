library(dplyr)
library(reshape2)
library(mgcv)
library(ggplot2)

#### SETUP ---------------------------------------------------------------------
source("R/simulation/0-parameters.R")
source("R/simulation/0-helper-functions.R")

sim_out = readRDS("output/simulation/sim_out.rds")

obs_df = sim_out$true_sims %>%
  filter(scenario_id == "T")

pred_mat_covar = expand.grid(vax_cov = new_vax_cov, 
                             location_id = unique(obs_df$location_id)) %>%
  left_join(unique(obs_df %>% dplyr::select(location_id, location_R0)))

#### IMPLEMENT APPROACH 1 ------------------------------------------------------
# step 1: find most plausible scenario for each location 
approach1_errors = sim_out$model_sims %>% 
  filter(scenario_id %in% c("S1", "S2")) %>%
  left_join(sim_out$true_sims %>% filter(scenario_id == "T") %>% rename(true_vax_cov = vax_cov) %>% dplyr::select(-scenario_id)) %>%
  mutate(vax_cov_diff = abs(vax_cov - true_vax_cov)) %>%
  mutate(min_vax_cov_diff = min(vax_cov_diff), .by = c("model_id", "location_id")) %>%
  mutate(most_plausible_scenario = ifelse(vax_cov_diff == min_vax_cov_diff, TRUE, FALSE)) %>%
  filter(most_plausible_scenario)
  
# step 2: calculate error by subtracting observation from projection in th is scenario
approach1_errors = approach1_errors %>%
  mutate(est_error = final_size - true_final_size)

#### IMPLEMENT APPROACH 2 (NO COVARIATES) --------------------------------------
# step 1: get errors in realized scenarios
error_df = sim_out$errors %>% filter(scenario_id == "T")

# step 2: fit errors in realized scenarios
gam_errors_nocov <- vector("list", n_models)
for(i in 1: n_models){
  gam_errors_nocov[[i]] <- gam(error ~ s(vax_cov), data = error_df %>% filter(model_id == paste0("M", i)))
}

gam_errors_nocov = lapply(gam_errors_nocov, get_gam_PIs, xvals = data.frame(vax_cov = new_vax_cov), invfn2 = function(x){return(x)})

approach2_nocov_errors_across_locs = bind_rows(gam_errors_nocov, .id = "model_id") %>%
  mutate(model_id = paste0("M", model_id)) %>%
  mutate(scenario_id = ifelse(vax_cov == vax_scenarios[1], "S1", ifelse(vax_cov == vax_scenarios[2], "S2", "E")))
saveRDS(approach2_nocov_errors_across_locs, "output/simulation/estimated_errors_approach2_full_across_locs.rds")

#### IMPLEMENT APPROACH 2 (WITH COVARIATES) ------------------------------------
# step 1: get errors in realized scenarios
# (above)

# step 2: fit errors in realized scenarios (with location cov)
gam_errors_cov <- vector("list", n_models)
approach2_obs_covariate <- vector("list", length(n_models))
for(i in 1: n_models){
  gam_errors_cov[[i]] <- gam(error ~ s(vax_cov) + location_R0, data = error_df %>% filter(model_id == paste0("M", i)))
  tmp <- vector("list", length(alphas))
  for(j in 1:nrow(pred_mat_covar)){
    tmp[[j]] =  get_gam_PIs(gam_errors_cov[[i]], xvals = pred_mat_covar[j, -2]) %>%
      mutate(location_id = pred_mat_covar[j, "location_id"])
  }
  approach2_obs_covariate[[i]] = bind_rows(tmp) %>%
    mutate(model_id = paste0("M", i))
}

# estimated errors for each location
approach2_cov_errors_all_locs = bind_rows(approach2_obs_covariate) %>%
  rename(vax_cov = xval.vax_cov, location_R0 = xval.location_R0) %>% 
  mutate(scenario_id = ifelse(vax_cov == vax_scenarios[1], "S1", ifelse(vax_cov == vax_scenarios[2], "S2", "E")))
saveRDS(approach2_cov_errors_all_locs, "output/simulation/estimated_errors_approach2_full_location_specific.rds")

# distribution of errors across locations
approach2_cov_errors_across_locs = approach2_cov_errors_all_locs %>%
  reframe(est_error_samp = get_samps(quantile, value, 1e4),
          draw_id = 1:1e4, .by = c("vax_cov", "location_id", "model_id")
  ) %>% 
  reframe(quantile = quantiles, 
            value = quantile(est_error_samp, quantiles), .by = c("vax_cov", "model_id")) %>%
  mutate(scenario_id = ifelse(vax_cov == vax_scenarios[1], "S1", ifelse(vax_cov == vax_scenarios[2], "S2", "E")))

#### IMPLEMENT APPROACH 3 (NO COVARIATES) --------------------------------------
# step 1: infer observations across locations
gam_obs = gam(true_final_size ~ s(vax_cov), data = obs_df)
approach3_obs_nocovariate = get_gam_PIs(gam_obs, xvals = data.frame(vax_cov = new_vax_cov)) %>%
  mutate(scenario_id = ifelse(vax_cov == new_vax_cov[1], "S1", "S2"))
# fit_obs_noR0 = fit_obs_noR0_long %>%
#   mutate(quantile = paste0("Q", quantile*100)) %>%
#   dcast(scenario_id + vax_cov ~ quantile)

# step 2: calculate error
n_samp = 1e4
# first get samples from the distribution of observations
approach3_samp_nocovariate = approach3_obs_nocovariate %>%
  reframe(est_obs_samp = get_samps(quantile, value, n_samp),
          draw_id = 1:n_samp, .by = c("vax_cov", "scenario_id")
  )

# join samples with model projections for each scenario and calculate error
approach3_nocov_errors = approach3_samp_nocovariate %>% 
  filter(vax_cov %in% c(0.3, 0.5)) %>%
  left_join(sim_out$model_sims %>% filter(scenario_id %in% c("S1", "S2")), 
            relationship = "many-to-many") %>%
  mutate(est_error = final_size - est_obs_samp)

approach3_nocov_errors_all_locs = approach3_nocov_errors %>% 
  reframe(quantile = quantiles, 
          value = quantile(est_error, quantiles), .by = c("vax_cov", "scenario_id", "model_id", "location_id"))

approach3_nocov_errors_across_locs = approach3_nocov_errors %>% 
  reframe(quantile = quantiles, 
          value = quantile(est_error, quantiles), .by = c("vax_cov", "scenario_id", "model_id"))

#### IMPLEMENT APPROACH 3 (WITH COVARIATES) ----------------------------------
# step 1: infer observation for each location
gam_obs_covar = gam(true_final_size ~ s(vax_cov) + location_R0, 
                    data = obs_df) 
approach3_obs_covariate <- vector("list", length(alphas))
for(i in 1:nrow(pred_mat_covar)){
  approach3_obs_covariate[[i]] =  get_gam_PIs(gam_obs_covar, xvals = pred_mat_covar[i, -2]) %>%
    mutate(location_id = pred_mat_covar[i, "location_id"])
}

approach3_obs_covariate = bind_rows(approach3_obs_covariate) %>%
  rename(vax_cov = xval.vax_cov, location_R0 = xval.location_R0) %>% 
  mutate(scenario_id = ifelse(vax_cov == new_vax_cov[1], "S1", "S2"))
saveRDS(approach3_obs_covariate, "output/simulation/estimated_obs_approach3.rds")

# step 2: calculate error
approach3_samp_covariate = approach3_obs_covariate %>%
  reframe(est_obs_samp = get_samps(quantile, value, 1e4),
          draw_id = 1:1e4, .by = c("vax_cov", "location_id")
  )

approach3_cov_errors = approach3_samp_covariate %>% 
  filter(vax_cov %in% c(0.3, 0.5)) %>%
  mutate(scenario_id = ifelse(vax_cov == vax_scenarios[1], "S1", ifelse(vax_cov == vax_scenarios[2], "S2", NA))) %>%
  left_join(sim_out$model_sims %>% filter(scenario_id %in% c("S1", "S2")), 
            relationship = "many-to-many") %>%
  mutate(est_error = final_size - est_obs_samp)

approach3_cov_errors_all_locs = approach3_cov_errors %>% 
  reframe(quantile = quantiles, 
          value = quantile(est_error, quantiles), .by = c("vax_cov", "scenario_id", "model_id", "location_id"))

approach3_cov_errors_across_locs = approach3_cov_errors %>% 
  reframe(quantile = quantiles, 
          value = quantile(est_error, quantiles), .by = c("vax_cov", "scenario_id", "model_id"))

#### COMBINE OUTCOMES ----------------------------------------------------------
# first estimates of error distribution across all locations
all_ests_across_locs = approach1_errors %>%
  reframe(quantile = quantiles, 
          value = quantile(est_error, quantiles), .by = c("vax_cov", "scenario_id", "model_id")) %>%
  mutate(approach = "1") %>%
  # add approach 2
  bind_rows(
    approach2_cov_errors_across_locs %>% mutate(approach = "2-covariates") %>% filter(scenario_id %in% c("S1", "S2"))
  ) %>%
  bind_rows(
    approach2_nocov_errors_across_locs %>% mutate(approach = "2-nocovariates") %>% filter(scenario_id %in% c("S1", "S2"))
  ) %>% 
  # add approach 3
  bind_rows(
      approach3_cov_errors_across_locs %>% mutate(approach = "3-covariates")
  ) %>%
  bind_rows(
    approach3_nocov_errors_across_locs %>% mutate(approach = "3-nocovariates")
  ) %>%
  # add truth
  bind_rows(
    sim_out$errors %>% 
      reframe(quantile = quantiles, 
              value = quantile(error, quantiles), .by = c("vax_cov", "model_id", "scenario_id")) %>%
      mutate(approach = "truth") %>% filter(scenario_id %in% c("S1", "S2"))
  )
saveRDS(all_ests_across_locs, "output/simulation/estimated_errors_across_locations.rds")
  
# second, estimates of error distribution in each location specifically
all_ests_all_locs = approach1_errors %>%
  rename(value = est_error) %>%
  select(vax_cov, scenario_id, model_id, location_id, value) %>%
  mutate(quantile = 0.5, 
         approach = "1") %>%
  # add approach 2
  bind_rows(
    approach2_cov_errors_all_locs %>% select(-location_R0) %>% mutate(approach = "2-covariates")
  ) %>%
  # add approach 3
  bind_rows(
    approach3_cov_errors_all_locs %>% mutate(approach = "3-covariates")
  ) %>%
  # add truth
  bind_rows(
    sim_out$errors %>% select(model_id, location_id, scenario_id, vax_cov, error) %>%
      mutate(quantile = 0.5, approach = "truth") %>%
      rename(value = error)  %>% filter(scenario_id %in% c("S1", "S2"))
  )
saveRDS(all_ests_all_locs, "output/simulation/estimated_errors_location_specific.rds")

