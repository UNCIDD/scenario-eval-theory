# vaccination scenarios
vax_scenarios = c(0.3, 0.5)
scenario_labs = c("low vax scenario", "high vax scenario")
names(scenario_labs) = c("S1", "S2")
# vaccination coverage values to estimate 
new_vax_cov = seq(0.3, 0.5, 0.05)

# number of locations
n_loc = 50
# number of models
n_models = 10

# quantiles of interest for summary error distribution
quantiles = c(0.001, 0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99, 0.999)
# and related alpha values 
alphas = sort(unique(sapply(quantiles, function(i){round(ifelse(i < 0.5, 1-2*i, 1-2*(1-i)),3)})))
# alphas = alphas[-which(alphas == 0)]

# plotting utilities
approach_labs = c("truth", "approach 1: most plausible scenario", 
                  "approach 2: estimate observations (no covariates)", "approach 2: estimate observations (covariates)", 
                  "approach 3: estimate error (no covariates)", "approach 3: estimate error (covariates)")
names(approach_labs) = c("truth", "1", "2-nocovariates", "2-covariates", "3-nocovariates", "3-covariates")

model_labs = paste0("model ", 1:n_models)
names(model_labs) = paste0("M", 1:n_models)
