library(dplyr)
library(reshape2)

#### SETUP ---------------------------------------------------------------------
source("R/simulation/sensitivity-runs/sens-0-parameters.R")
source("R/simulation/0-helper-functions.R")

# load estimates
all_ests_across_locs = readRDS("output/simulation/estimated_errors_across_locations_sens.rds") 
all_ests_all_locs = readRDS("output/simulation/estimated_errors_location_specific_sens.rds")

sim_out = readRDS("output/simulation/sim_out_sens.rds")
true_errors = bind_rows(lapply(sim_out, function(i){i$errors}), .id = "rep_id") %>%
  mutate(rep_id = as.integer(rep_id))

# in each case, we will use three metrics to evaluate performance:
# (1) mean absolute error, (2) KS statistic

#### EVALUATE PERFORMANCE OF DISTRIBUTION ACROSS LOCATIONS ---------------------
set.seed(0)
method_comparsion_results = expand.grid(model_id = paste0("M", 1:n_models), 
                                        scenario_id = c("S1", "S2"), 
                                        rep_id = 1:n_reps,
                                        approach = names(approach_labs[which(approach_labs != "truth")]), 
                                        n_df = NA, n_truth = NA,
                                        ks_test_stat = NA, ks_test_p = NA, mae = NA)
for(i in 1:nrow(method_comparsion_results)){
  # print(i)
  if(i %%100 == 0){print(paste0(i, "/", nrow(method_comparsion_results)))}
  tmp_scenario_id = method_comparsion_results[i, "scenario_id"]
  tmp_model_id = method_comparsion_results[i, "model_id"]
  tmp_method = method_comparsion_results[i, "approach"]
  tmp_rep = method_comparsion_results[i, "rep_id"]
  # pull error distribution of interest
  tmp_df = all_ests_across_locs %>% 
    filter(scenario_id == tmp_scenario_id, 
           model_id == tmp_model_id, 
           approach == tmp_method,
           rep_id == tmp_rep)
  # get samples
  tmp_df = tmp_df %>%
    reframe(est_error = get_samps(quantile, value, 1e4),
            draw_id = 1:1e4, .by = c("vax_cov", "model_id")
    )
  # get true error
  tmp_truth = true_errors %>% filter(scenario_id == as.character(tmp_scenario_id), model_id == tmp_model_id, rep_id == tmp_rep)
  # perform KS test
  tmp_ks = ks.test(tmp_df$est_error, tmp_truth$error)
  # save output
  method_comparsion_results[i, "n_df"] = length(tmp_df$est_error)
  method_comparsion_results[i, "n_truth"] = length(tmp_truth$error)
  method_comparsion_results[i, "ks_test_stat"] = tmp_ks$statistic
  method_comparsion_results[i, "ks_test_p"] = tmp_ks$p.value
  method_comparsion_results[i, "mae"] = abs(mean(tmp_df$est_error) - mean(tmp_truth$error))
  # method_comparsion_results[i, "ks_p5_lvl"] = 1.358*sqrt((n_df + n_truth)/(n_df*n_truth)) # level of ks statistic s.t., p-value = 0.05
}

saveRDS(method_comparsion_results, "output/simulation/performance_evaluation_results_across_locations_sens.rds")


#### EVALUATE PERFORMANCE OF DISTRIBUTION FOR EACH LOCATION --------------------
set.seed(0)
method_comparsion_results = expand.grid(model_id = paste0("M", 1:n_models),
                                        scenario_id = c("S1", "S2"),
                                        location_id = 1:n_loc,
                                        rep_id = 1:n_reps,
                                        approach = names(approach_labs[c(2, 4, 6)]),
                                        n_df = NA, n_truth = NA,
                                        ks_test_stat = NA, ks_test_p = NA, mae = NA)
for(i in 1:nrow(method_comparsion_results)){
  if(i %%100 == 0){print(paste0(i, "/", nrow(method_comparsion_results)))}
  browser()
  tmp_scenario_id = method_comparsion_results[i, "scenario_id"]
  tmp_model_id = method_comparsion_results[i, "model_id"]
  tmp_method = method_comparsion_results[i, "approach"]
  tmp_location_id = method_comparsion_results[i, "location_id"]
  tmp_rep = method_comparsion_results[i, "rep_id"]
  # get true error
  tmp_truth = true_errors %>%
    filter(scenario_id == as.character(tmp_scenario_id), model_id == tmp_model_id,
           location_id == tmp_location_id, rep_id == tmp_rep)
  # pull error distribution of interest
  tmp_df = all_ests_all_locs %>%
    filter(scenario_id == tmp_scenario_id, model_id == tmp_model_id,
           approach == tmp_method, location_id == tmp_location_id, rep_id == tmp_rep)
  if(nrow(tmp_df) == 0){next}
  if(nrow(tmp_df) == 1){method_comparsion_results[i, "mae"] = abs(tmp_df$value - tmp_truth$error);next}
  # get samples
  tmp_df = tmp_df %>%
    reframe(est_error = get_samps(quantile, value, 1e4),
            draw_id = 1:1e4, .by = c("vax_cov", "model_id", "location_id", "rep_id")
    )
  # perform KS test
  tmp_ks = ks.test(tmp_df$est_error, tmp_truth$error)
  # save output
  method_comparsion_results[i, "n_df"] = length(tmp_df$est_error)
  method_comparsion_results[i, "n_truth"] = length(tmp_truth$error)
  method_comparsion_results[i, "ks_test_stat"] = tmp_ks$statistic
  method_comparsion_results[i, "ks_test_p"] = tmp_ks$p.value
  method_comparsion_results[i, "mae"] = abs(mean(tmp_df$est_error) - mean(tmp_truth$error))
  # method_comparsion_results[i, "ks_p5_lvl"] = 1.358*sqrt((n_df + n_truth)/(n_df*n_truth)) # level of ks statistic s.t., p-value = 0.05
}

saveRDS(method_comparsion_results, "output/simulation/performance_evaluation_results_location_specific_sens.rds")

