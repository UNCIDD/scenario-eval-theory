#### FUNCTIONS -----------------------------------------------------------------
filter_posterior <- function(posterior, filter_criteria){
  if(length(dim(posterior)) == 2) {posterior_filtered <- posterior[,grepl(filter_criteria, dimnames(posterior)[[2]])]}
  else if(length(dim(posterior)) == 3) {posterior_filtered <- posterior[,,grepl(filter_criteria, dimnames(posterior)$parameters)]}
  posterior_long <- as.data.frame(posterior_filtered) %>%
    mutate(draw_id = 1:nrow(as.data.frame(posterior))) %>%
    melt(c("draw_id")) %>%
    mutate(index_id = as.integer(substr(as.character(variable), nchar(filter_criteria), nchar(as.character(variable))-1)))
  return(posterior_long)
}

boxcox_transform <- function(y, lambda) {
  if (lambda == 0) log(y) else (y^lambda - 1) / lambda
} 

inverse_boxcox <- function(y_trans, lambda) {
  if (lambda == 0) exp(y_trans) else (lambda * y_trans + 1)^(1/lambda)
}

get_preds <- function(fit, new_vax_cov){
  y_exclude = filter_posterior(as.data.frame(fit), "y_exclude\\[") %>%
    rename(exclude = value) %>%
    dplyr::select(-variable)
  y_nominal = filter_posterior(as.data.frame(fit), "y_new_nominal\\[") %>%
    left_join(y_exclude) %>%
    left_join(data.frame(index_id = 1:length(new_vax_cov), 
                         vax_cov = new_vax_cov)) %>%
    filter(exclude == 0)
  return(y_nominal)
}

summarize_predints <- function(fit, new_vax_cov){
  y_nominal = get_preds(fit, new_vax_cov) %>%
    summarize(Q5 = quantile(value, 0.05), 
              Q10 = quantile(value, 0.10),
              Q15 = quantile(value, 0.15), 
              Q20 = quantile(value, 0.20), 
              Q25 = quantile(value, 0.25),
              Q30 = quantile(value, 0.30), 
              Q35 = quantile(value, 0.35), 
              Q40 = quantile(value, 0.40), 
              Q45 = quantile(value, 0.45), 
              Q50 = quantile(value, 0.5), 
              Q55 = quantile(value, 0.55), 
              Q60 = quantile(value, 0.60), 
              Q65 = quantile(value, 0.65), 
              Q70 = quantile(value, 0.70), 
              Q75 = quantile(value, 0.75),
              Q80 = quantile(value, 0.80), 
              Q85 = quantile(value, 0.85), 
              Q90 = quantile(value, 0.90), 
              Q95 = quantile(value, 0.95), .by = c("vax_cov"))
}

get_samps <- function(quantile, value, n_samps = 1e4, seed = 1002){
  set.seed(seed)
  approx(quantile, value, runif(n_samps), yleft = min(value), yright = max(value))$y
}


#### SETUP ---------------------------------------------------------------------
library(finalsize)
library(dplyr)
library(reshape2)
library(ggplot2)
library(quantreg)
library(gamlss)
library(bestNormalize)
library(MASS)
library(rstan)
library(rstanarm)
library(bayesplot)
library(readr)

source("./R/final-size-functions.R")

scenario_labs = c("low vax scenario", "high vax scenario")
names(scenario_labs) = c("S1", "S2")

#### GENERATE SIMULATIONS TO TEST ----------------------------------------------
reps = 1
set.seed(10)
seed_id = sample(1:1E6, reps)
n_models = 10

# quantiles of interest for summary error distribution
quantiles = c(0.001, 0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99, 0.999)
# and related alpha values 
alphas = sort(unique(sapply(quantiles, function(i){round(ifelse(i < 0.5, 1-2*i, 1-2*(1-i)),3)})))
# alphas = alphas[-which(alphas == 0)]

# new vaccination coverage
new_vax_cov = seq(0.3, 0.5, 0.05)

# run simulation - using 500 locations for now
t_large <- full_sim(n_locations = 500, n_models = n_models,
                    vax_cov_S1 = 0.3, vax_cov_S2 = 0.5,
                    R0_lwr = 2, R0_upr = 3,
                    model_bias_ind_sd = 0.05,
                    seed = seed_id, fit_outcomes = FALSE)

# note: scenario_id = T returns only observed errors (not true errors to test against)
error_df = t_large$errors %>% filter(scenario_id == "T") 

# do box-cox ahead of time (for now)
# bx = boxcox(error ~ vax_cov, data = error_df_sub, plotit = FALSE)
# lambda = with(bx, x[which.max(y)])
# 
# error_df_sub$error_boxcox = ifelse(error_df_sub < 0, boxcox_transform(-error_df_sub$error, lambda), boxcox_transform(error_df_sub$error, lambda))
# 
# ggplot(data = error_df_sub, aes(x = vax_cov, y = error_boxcox)) + 
#   geom_point()
# 
# # try with lm
# fit_lm <- lm(error_boxcox ~ vax_cov, error_df_sub)

#### FIT WITH STAN -------------------------------------------------------------
mixture_model <- stan_model("R/sim-experiment-final_size/mixture-distribution/mixture_model.stan")

fit_errors <- vector("list", n_models)

for(i in 1:n_models){
  warning(paste0("fitting M", i))
  print(paste0("fitting M", i))
  error_df_sub = error_df %>% 
    filter(model_id == paste0("M", i)) #%>%
    #mutate(error_boxcox = ifelse(error < 0, boxcox_transform(-error, lambda), boxcox_transform(error, lambda)))
  # fit with STAN
  fit_errors[[i]] <- sampling(
    mixture_model,
    data = list(
      N = nrow(error_df_sub),
      x = error_df_sub$vax_cov,
      y = error_df_sub$error, 
      N_new = length(new_vax_cov), 
      x_new = new_vax_cov
    ),
    seed = 7, 
    iter = 10000,
    chain = 4, 
    cores = 4, 
    control=list(max_treedepth = 12)
  )
}
error <- names(warnings())
out <- file("R/sim-experiment-final_size/mixture-distribution/warnings_noshrinkage.txt")
writeLines(error, out)
close(out)

write_rds(fit_errors, "R/sim-experiment-final_size/mixture-distribution/fit_noshrinkage.rda")


beepr::beep()

#### FIT WITH SHRINKAGE PARAMETERS ---------------------------------------------
mixture_model_shrinkage <- stan_model("R/sim-experiment-final_size/mixture-distribution/mixture_model_shrinkage.stan")

fit_errors_shrinkage <- vector("list", n_models)

for(i in 1:n_models){
  warning(paste0("fitting M", i))
  print(paste0("fitting M", i))
  error_df_sub = error_df %>% 
    filter(model_id == paste0("M", i)) #%>%
  #mutate(error_boxcox = ifelse(error < 0, boxcox_transform(-error, lambda), boxcox_transform(error, lambda)))
  # fit with STAN
  fit_errors_shrinkage[[i]] <- sampling(
    mixture_model_shrinkage,
    data = list(
      N = nrow(error_df_sub),
      x = error_df_sub$vax_cov,
      y = error_df_sub$error, 
      N_new = length(new_vax_cov), 
      x_new = new_vax_cov
    ),
    seed = 7, 
    iter = 10000,
    chain = 4, 
    cores = 4, 
    control=list(max_treedepth = 12)
  )
}

error <- names(warnings())
out <- file("R/sim-experiment-final_size/mixture-distribution/warnings_shrinkage.txt")
writeLines(error, out)
close(out)

write_rds(fit_errors_shrinkage, "R/sim-experiment-final_size/mixture-distribution/fit_shrinkage.rda")

# some summary/diagnostics
posterior <- as.array(fit_errors_hierparams[[4]])
posterior_df <- as.data.frame(fit_errors_hierparams[[4]])
np <- nuts_params(fit_errors_hierparams[[4]])

mcmc_pairs(fit_errors_hierparams[[4]], pars = c("alpha_neg", "beta_neg","lambda_neg", "sigma_neg","lambda_pos", "alpha_pos", "beta_pos", "sigma_pos", "p"), np = np)

#### GET PARAMETER ESTIMATES ---------------------------------------------------
pars_to_extract = c("alpha_pos", "alpha_neg", "beta_pos", "beta_neg", "lambda_pos", "lambda_neg", "p")

pars <- lapply(fit_errors, 
                         function(i){bind_cols(rstan::extract(i, pars_to_extract)) %>% mutate(draw_id = seq_len(n()))}) %>%
  bind_rows(.id = "model_id") %>%
  mutate(model_id = paste0("M", model_id))

pars_shrinkage <- lapply(fit_errors_shrinkage, 
                         function(i){bind_cols(rstan::extract(i, pars_to_extract)) %>% mutate(draw_id = seq_len(n()))}) %>%
  bind_rows(.id = "model_id") %>%
  mutate(model_id = paste0("M", model_id))

# plot distributions with and without shrinkage
bind_rows(pars %>% mutate(fit = "no shrinkage"),
          pars_shrinkage %>% mutate(fit = "shrinkage")) %>%
  pivot_longer(!c("fit", "model_id", "draw_id"), names_to = "variable") %>%
  separate(variable, into = c("variable", "sign")) %>%
  mutate(sign = ifelse(is.na(sign), "pos", sign)) %>%
  ggplot(aes(x = value, color = sign, linetype = fit, group = interaction(fit, sign))) +
  geom_density(linewidth = 0.8) + 
  # facet_grid(cols = vars(model_id), rows = vars(variable), scales = "free") +
  facet_wrap(vars(variable,model_id), ncol = 10, scales = "free") +
  scale_color_brewer(palette = "Set1") + 
  scale_linetype_manual(values = c("dotted", "solid")) + 
  theme_bw() +
  theme(legend.position = "bottom", 
        panel.grid.minor = element_blank())

# get mean estimates
mean_ests = lapply(fit_errors, function(i){unlist(lapply(rstan::extract(i, pars_to_extract), mean))}) %>%
  bind_rows(.id = "model_id") %>%
  mutate(model_id = paste0("M", model_id))

mean_ests_shrink = lapply(fit_errors_shrinkage, function(i){unlist(lapply(rstan::extract(i, pars_to_extract), mean))}) %>%
  bind_rows(.id = "model_id") %>%
  mutate(model_id = paste0("M", model_id))

#### GET PREDICTION INTERVALS --------------------------------------------------
pred_intervals <- lapply(fit_errors, summarize_predints, new_vax_cov = new_vax_cov) %>%
  bind_rows(.id = "model_id") %>%
  mutate(model_id = paste0("M", model_id))

pred_intervals_shrink <- lapply(fit_errors_shrinkage, summarize_predints, new_vax_cov = new_vax_cov) %>%
  bind_rows(.id = "model_id") %>%
  mutate(model_id = paste0("M", model_id))

# nominal scale
# ggplot(data = pred_intervals_shrink, aes(x = vax_cov)) + 
#   geom_point(data = error_df, aes(x = vax_cov, y = error), color = "black", shape = 21) + 
#   geom_ribbon(aes(ymin = Q5, ymax = Q95, fill = model_id), alpha = 0.4) + 
#   geom_ribbon(aes(ymin = Q25, ymax = Q75, fill = model_id), alpha = 0.6) + 
#   geom_line(aes(y = Q50, color = model_id), size = 1) + 
#   facet_wrap(vars(model_id), scales = "free") + 
#   theme_bw() + 
#   theme(legend.position = "none")

# predictions with and without shrinkage
bind_rows(pred_intervals %>% mutate(fit = "no shrinkage"),
          pred_intervals_shrink %>% mutate(fit = "shrinkage")) %>%
  ggplot(aes(x = vax_cov, color = model_id)) + 
  geom_point(data = error_df, aes(x = vax_cov, y = error), color = "darkgray", shape = 21) + 
  geom_line(aes(y = Q5, linetype = fit), linewidth = 1, alpha = 0.7) +
  geom_line(aes(y = Q95, linetype = fit), linewidth = 1, alpha = 0.7) +
  geom_line(aes(y = Q25, linetype = fit), linewidth = 1, alpha = 0.7) +
  geom_line(aes(y = Q75, linetype = fit), linewidth = 1, alpha = 0.7) +
  geom_line(aes(y = Q50, linetype = fit), linewidth = 2, alpha = 0.7) +
  facet_wrap(vars(model_id), scales = "free") + 
  guides(color = FALSE) + 
  theme_bw() + 
  theme(legend.position = "bottom")

# check coverage
cov = calculate_coverage(
  # do some reshaping to match expected format for quant_fits
  quant_fits = pred_intervals %>%
    melt(c("model_id", "vax_cov")) %>% 
    filter(variable != "mean") %>%
    mutate(quantile = as.integer(gsub("Q", "", variable))/100) %>%
    dplyr::select(-variable), 
  error_df = t_large$errors,  # %>% mutate(error = abs(error))
  vax_cov_S1 = 0.3, vax_cov_S2 = 0.5, 
  summarize_by = "model"
)

cov_shrink = calculate_coverage(
  # do some reshaping to match expected format for quant_fits
  quant_fits = pred_intervals_shrink %>%
    melt(c("model_id", "vax_cov")) %>% 
    filter(variable != "mean") %>%
    mutate(quantile = as.integer(gsub("Q", "", variable))/100) %>%
    dplyr::select(-variable), 
  error_df = t_large$errors,  # %>% mutate(error = abs(error))
  vax_cov_S1 = 0.3, vax_cov_S2 = 0.5, 
  summarize_by = "model"
)

# ggplot(data = cov, aes(x = alpha, y = cov, group = model_id)) + 
#   geom_line(aes(color = model_id)) + 
#   geom_abline(size = 1) + 
#   facet_wrap(vars(scenario_id), labeller = labeller(scenario_id = scenario_labs), ncol = 1) + 
#   labs(x = "expected coverage", y = "actual coverage") + 
#   theme_bw() + 
#   theme(legend.position = "none")

bind_rows(cov %>% mutate(fit = "no shrinkage"),
          cov_shrink %>% mutate(fit = "shrinkage")) %>%
  ggplot(aes(x = alpha, y = cov, group = interaction(model_id, fit))) + 
  geom_line(aes(color = model_id, linetype = fit), alpha = 0.7) + 
  geom_abline(size = 1) + 
  facet_wrap(vars(scenario_id), labeller = labeller(scenario_id = scenario_labs), ncol = 1) + 
  labs(x = "expected coverage", y = "actual coverage") + 
  theme_bw() + 
  theme(legend.position = "none")


#### REPEAT WITH FEWER LOCATIONS -----------------------------------------------
# run simulation - using 500 locations for now
t_small <- full_sim(n_locations = 50, n_models = n_models,
                    vax_cov_S1 = 0.3, vax_cov_S2 = 0.5,
                    R0_lwr = 2, R0_upr = 3,
                    model_bias_ind_sd = 0.05,
                    seed = seed_id, fit_outcomes = FALSE)

# note: scenario_id = T returns only observed errors (not true errors to test against)
error_df_small = t_small$errors %>% filter(scenario_id == "T")

#### FIT WITH SHRINKAGE PARAMETERS ---------------------------------------------
fit_errors_shrinkage_small <- vector("list", n_models)

for(i in 1:n_models){
  warning(paste0("fitting M", i))
  print(paste0("fitting M", i))
  error_df_sub = error_df_small %>% 
    filter(model_id == paste0("M", i)) #%>%
  #mutate(error_boxcox = ifelse(error < 0, boxcox_transform(-error, lambda), boxcox_transform(error, lambda)))
  # fit with STAN
  fit_errors_shrinkage_small[[i]] <- sampling(
    mixture_model_shrinkage,
    data = list(
      N = nrow(error_df_sub),
      x = error_df_sub$vax_cov,
      y = error_df_sub$error, 
      N_new = length(new_vax_cov), 
      x_new = new_vax_cov
    ),
    seed = 7, 
    iter = 10000,
    chain = 4, 
    cores = 4, 
    control=list(max_treedepth = 12)
  )
}

error <- names(warnings())
out <- file("R/sim-experiment-final_size/mixture-distribution/small_warnings_shrinkage.txt")
writeLines(error, out)
close(out)

write_rds(fit_errors_shrinkage_small, "R/sim-experiment-final_size/mixture-distribution/small_fit_shrinkage.rda")

fit_errors_shrinkage_small <- read_rds("R/sim-experiment-final_size/mixture-distribution/small_fit_shrinkage.rda")

#### GET PREDICTION INTERVALS --------------------------------------------------
pred_intervals_shrink_small <- lapply(fit_errors_shrinkage_small, summarize_predints, new_vax_cov = new_vax_cov) %>%
  bind_rows(.id = "model_id") %>%
  mutate(model_id = paste0("M", model_id))

# nominal scale
ggplot(data = pred_intervals_shrink_small, aes(x = vax_cov)) +
  geom_point(data = error_df_small, aes(x = vax_cov, y = error), color = "black", shape = 21) +
  geom_ribbon(aes(ymin = Q5, ymax = Q95, fill = model_id), alpha = 0.4) +
  geom_ribbon(aes(ymin = Q25, ymax = Q75, fill = model_id), alpha = 0.6) +
  geom_line(aes(y = Q50, color = model_id), size = 1) +
  facet_wrap(vars(model_id), scales = "free") +
  theme_bw() +
  theme(legend.position = "none")

cov_shrink_small = calculate_coverage(
  # do some reshaping to match expected format for quant_fits
  quant_fits = pred_intervals_shrink_small %>%
    melt(c("model_id", "vax_cov")) %>% 
    filter(variable != "mean") %>%
    mutate(quantile = as.integer(gsub("Q", "", variable))/100) %>%
    dplyr::select(-variable), 
  error_df = t_small$errors,  # %>% mutate(error = abs(error))
  vax_cov_S1 = 0.3, vax_cov_S2 = 0.5, 
  summarize_by = "model"
)

ggplot(data = cov_shrink_small, aes(x = alpha, y = cov, group = model_id)) +
  geom_line(aes(color = model_id)) +
  geom_abline(size = 1) +
  facet_wrap(vars(scenario_id), labeller = labeller(scenario_id = scenario_labs), ncol = 2) +
  labs(x = "expected coverage", y = "actual coverage") +
  theme_bw() +
  theme(legend.position = "none", 
        legend.title = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank())
ggsave("R/sim-experiment-final_size/mixture-distribution/coverage_nlocationssmall.pdf", width = 5, height = 3)

#### MAKE A FIGURE TO ILLUSTRATE THE APPROACH ----------------------------------
locs_to_plot = c(29, 5, 45)
loc_labs = c(15, 41, 32)
loc_cols = RColorBrewer::brewer.pal(4, "Set1")[2:4]
model_to_plot = "M4"

# overall results
mod_results = ggplot(data = pred_intervals_shrink_small %>% filter(model_id == model_to_plot), 
                    aes(x = vax_cov)) +
  geom_point(data = error_df_small %>% filter(model_id == model_to_plot), 
             aes(x = vax_cov, y = error), color = "black", shape = 21, size = 2.5) +
  geom_ribbon(aes(ymin = Q5, ymax = Q95, fill = model_id), alpha = 0.15, fill = "darkgray") +
  geom_ribbon(aes(ymin = Q25, ymax = Q75, fill = model_id), alpha = 0.3, fill = "darkgray") +
  geom_line(aes(y = Q50, color = model_id), alpha = 0.8, size = 1.5, color = "darkgray") +
  geom_point(data = error_df_small %>% filter(model_id == model_to_plot, location_id %in% locs_to_plot), 
            aes(y = error, color = as.factor(location_id)), size = 2.5) + 
  geom_text(data = error_df_small%>% 
              filter(model_id == model_to_plot, location_id %in% locs_to_plot) %>% 
              mutate(loc_lab = loc_labs[which(locs_to_plot == location_id)], .by = "location_id"),
            aes(y = error, label = loc_lab), size = 1.8, color = "white") +
  facet_wrap(vars(model_id), scales = "free") +
  labs(x = "realized vaccine uptake\n(scenario axis)") +
  scale_color_manual(values = loc_cols) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.minor =  element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_blank())
mod_results

# a few specific locations
loc_results = ggplot(data = t_small$model_sims %>% 
                       filter(model_id %in% c(model_to_plot, "T"), location_id %in% locs_to_plot, 
                              scenario_id != "E") %>%
                       mutate(location_id = factor(location_id, levels = locs_to_plot))) + 
  geom_line(data = t_small$model_sims %>% 
              filter(model_id == model_to_plot, location_id %in% locs_to_plot) %>%
              mutate(location_id = factor(location_id, levels = locs_to_plot)), 
             aes(x = vax_cov, y = final_size), alpha = 0.2) + 
  geom_point(aes(x = vax_cov, y = final_size, shape = scenario_id, fill = scenario_id), color = 'black', size = 2) +
  geom_text(data = data.frame(location_id = locs_to_plot, 
                              location_lab = loc_labs) %>%
              mutate(location_id = factor(location_id, levels = locs_to_plot)),
            aes(x = Inf, y = Inf, label = paste0("location ", location_lab), color = as.factor(location_id)), 
            hjust = 1, vjust = 1, size = 3.5) +
  geom_point(data = t_small$true_sims %>% filter(location_id %in% locs_to_plot, scenario_id == "T"), 
             aes(x = vax_cov, y = true_final_size), color = "red", size = 2) + 
  geom_segment(data = left_join(
    t_small$model_sims %>% filter(model_id == model_to_plot, location_id %in% locs_to_plot, scenario_id == "T"),
    t_small$true_sims %>% filter(location_id %in% locs_to_plot, scenario_id == "T")
  ), aes(x = vax_cov, xend = vax_cov, y = true_final_size, yend = final_size), arrow = arrow(length = unit(0.05, "npc")), size = 0.4) +
  facet_wrap(vars(location_id), ncol = 1) +
  labs(x = "realized vaccine uptake\n(scenario axis)", 
       y = "cumulative hospitalizations\n(projection axis)") +
  scale_color_manual(values = loc_cols[sapply(sort(locs_to_plot), function(i){which(locs_to_plot == i)})]) +
  scale_fill_manual(values = c("black", "black", "white")) +
  scale_shape_manual(values = c(16, 16, 21)) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank())
loc_results

cowplot::plot_grid(loc_results, mod_results, rel_widths = c(0.32, 0.68), 
                   labels = c("A", "B"))

ggsave("R/sim-experiment-final_size/mixture-distribution/approach_illustration.pdf", width = 7, height = 3.75)

#### EXAMPLE SCENARIO ERROR CALCULATION ----------------------------------------
# total error = calibration error + scenario error
# calibration error: estimated by the STAN fits
# total error (can be calculated directly)
# thus, scenario error = total error - calibration error
calibration_error <- lapply(fit_errors_shrinkage_small, get_preds, new_vax_cov = new_vax_cov)
total_error = t_small$model_sims %>% 
  filter(scenario_id %in% c("S1", "S2")) %>%
  dplyr::select(model_id, location_id, scenario_id, final_size) %>%
  left_join(t_small$true_sims %>% filter(scenario_id == c("T")) %>%
              dplyr::select(-scenario_id)) %>%
  mutate(total_error = final_size - true_final_size) %>% 
  dplyr::select(model_id, location_id, scenario_id, total_error)

scenario_error = calibration_error[as.integer(substr(model_to_plot,2,2))][[1]] %>%
  dplyr::select(-draw_id, -variable, -exclude) %>%
  rename(calibration_error = value) %>%
  mutate(scenario_id = ifelse(vax_cov == new_vax_cov[1], "S1", "S2"), 
         model_id = model_to_plot) %>%
  left_join(total_error %>% filter(model_id == model_to_plot), 
            relationship = "many-to-many", by = join_by(model_id, scenario_id)) %>%
  mutate(scenario_error = total_error - calibration_error)

scenario_error %>% 
  filter(!is.na(scenario_error)) %>%
  summarize(median = median(scenario_error), 
            lwr = quantile(scenario_error, 0.05), 
            upr = quantile(scenario_error, 0.95), 
            .by = c("scenario_id", "location_id", "model_id")) %>%
  mutate(true_slope = true_slope[location_id]) %>%
  ggplot(aes(x = true_slope, color = as.factor(location_id))) + 
  geom_hline(aes(yintercept = 0)) + 
  geom_point(aes(y = median), size = 3) + 
  geom_segment(aes(xend = true_slope, y = lwr, yend = upr), size = 1) + 
  geom_text(aes(y = median, label = location_id), size = 2.5, color = "white") + 
  facet_grid(cols = vars(model_id), rows = vars(scenario_id)) + 
  labs(x = "true slope of relationship in location", y = "scenario error") +
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid.minor = element_blank())


#### ESTIMATED ERROR DISTRIBUTION FOR EACH MODEL -------------------------------
# QUESTION: are we implicitly using different "observation distributions" for
# each model? 
calibration_error %>%
  bind_rows(.id = "model_id") %>%
  mutate(model_id = paste0("M", model_id)) %>%
  left_join(t_small$model_sims, 
            relationship = "many-to-many") %>%
  mutate(obs_est = final_size - value) %>%
  filter(scenario_id %in% c("S1", "S2"), obs_est < 1, obs_est > 0) %>%
  ggplot() + 
  geom_histogram(aes(x = obs_est, y = after_stat(density), fill = model_id)) + 
  geom_density(data = t_small$true_sims %>% filter(scenario_id %in% c("S1", "S2")), 
               aes(x = true_final_size), color = "black", size = 1) +
  facet_grid(cols = vars(scenario_id), rows = vars(model_id), scales = "free") +
   theme_bw() + 
  theme(legend.position = 'none')

# use KS test (vs. 50 random draws from the distribution)
true_obs_S1 = t_small$true_sims %>% 
  filter(scenario_id  == "S1") %>%
  pull(true_final_size)
true_obs_S2 = t_small$true_sims %>% 
  filter(scenario_id  == "S2") %>%
  pull(true_final_size)

est_obs_S1 = calibration_error %>%
  bind_rows(.id = "model_id") %>%
  mutate(model_id = paste0("M", model_id)) %>%
  left_join(t_small$model_sims, 
            relationship = "many-to-many") %>%
  mutate(obs_est = final_size - value) %>%
  filter(scenario_id == "S1",  obs_est < 1, obs_est > 0) # think about implications of this exclusion
est_obs_S2 = calibration_error %>%
  bind_rows(.id = "model_id") %>%
  mutate(model_id = paste0("M", model_id)) %>%
  left_join(t_small$model_sims, 
            relationship = "many-to-many") %>%
  mutate(obs_est = final_size - value) %>%
  filter(scenario_id == "S2", obs_est < 1, obs_est > 0)


mod_ks_rslt = vector("list", n_models)
for(i in 1:n_models){
  print(i)
  mod_ks_rslt[[i]] = data.frame(scenario_id = c("S1", "S2"), 
                         ks_statistic = NA, 
                         ks_pvalue = NA)
  ks_S1 = ks.test(est_obs_S1 %>% filter(model_id == paste0("M", i)) %>% pull(obs_est), 
                  true_obs_S1)
  ks_S2 = ks.test(est_obs_S2 %>% filter(model_id == paste0("M", i)) %>% pull(obs_est), 
                  true_obs_S2)
  mod_ks_rslt[[i]][, "ks_statistic"] = c(ks_S1$statistic, ks_S2$statistic)
  mod_ks_rslt[[i]][, "ks_pvalue"] = c(ks_S1$p.value, ks_S2$p.value)
 }
mod_ks_rslt = bind_rows(mod_ks_rslt, .id = "model_id") %>%
  mutate(model_id = paste0("M", model_id))

set.seed(222)
n_rand_samp = 100
rand_ks_rslt = vector("list", n_rand_samp)
for(j in 1:n_rand_samp){
  rand_ks_rslt[[j]] = data.frame(scenario_id = c("S1", "S2"), 
                                ks_statistic = NA, 
                                ks_pvalue = NA)
  rand_samp_S1 = sample(true_obs_S1, length(true_obs_S1), replace = TRUE)
  rand_samp_S2 = sample(true_obs_S2, length(true_obs_S2), replace = TRUE)
  ks_S1 = ks.test(rand_samp_S1, true_obs_S1)
  ks_S2 = ks.test(rand_samp_S2, true_obs_S2)
  rand_ks_rslt[[j]][, "ks_statistic"] = c(ks_S1$statistic, ks_S2$statistic)
  rand_ks_rslt[[j]][, "ks_pvalue"] = c(ks_S1$p.value, ks_S2$p.value)
}
rand_ks_rslt = bind_rows(rand_ks_rslt, .id = "rand_id")

# plot random results and individual model results
ggplot(data = rand_ks_rslt, aes(x = scenario_id, y = ks_statistic)) +
  geom_violin(alpha = 0.25, fill = "gray", color = NA) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), shape = 21) + 
  geom_point(data = mod_ks_rslt, aes(color = model_id), 
             position = position_jitter(seed = 1, width = 0.2), size = 2) + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid = element_blank())


### FIT OBSERVATIONS INSTEAD OF ERRORS (NO COVARIATES) -------------------------
obs_df_small = t_small$true_sims %>%
  filter(scenario_id == "T")

ggplot(data = obs_df_small, aes(x = vax_cov,  y = true_final_size)) + 
  geom_point()

hist(obs_df_small$true_final_size)

qqnorm(obs_df_small$true_final_size)
qqline(obs_df_small$true_final_size)

# log transform doesn't work
qqnorm(log(obs_df_small$true_final_size))
qqline(log(obs_df_small$true_final_size))

# sqrt transform doesn't work
qqnorm(sqrt(obs_df_small$true_final_size))
qqline(sqrt(obs_df_small$true_final_size))

# normalization doesn't work
qqnorm(scale(obs_df_small$true_final_size))
qqline(scale(obs_df_small$true_final_size))

# try box-cox, maybe slightly better
qqnorm(bestNormalize::boxcox(obs_df_small$true_final_size)$x.t)
qqline(bestNormalize::boxcox(obs_df_small$true_final_size)$x.t)

bx = boxcox(true_final_size ~ vax_cov, data = obs_df_small, plotit = TRUE)
lambda = with(bx, x[which.max(y)])
lm_obs = lm(boxcox_transform(true_final_size, lambda) ~ vax_cov, data = obs_df_small)
fit_obs_noR0 <- vector("list", length(alphas))
for(i in 1:length(alphas)){
  fit_obs_noR0[[i]] <- predict(lm_obs, newdata = data.frame(vax_cov = new_vax_cov), 
                  level = alphas[i], interval = "prediction") %>%
    as.data.frame() %>% 
    mutate(vax_cov = new_vax_cov, 
           alpha = alphas[i],
           fit = ifelse(is.na(inverse_boxcox(fit, lambda)), 0, inverse_boxcox(fit, lambda)),
           lwr = ifelse(is.na(inverse_boxcox(lwr, lambda)), 0, inverse_boxcox(lwr, lambda)), 
           upr = ifelse(is.na(inverse_boxcox(upr, lambda)), 0, inverse_boxcox(upr, lambda)))
}
fit_obs_noR0_long = bind_rows(fit_obs_noR0) %>% 
  reshape2::melt(c("vax_cov", "fit", "alpha")) %>%
  mutate(quantile = ifelse(variable == "lwr", (1-alpha)/2, 1-(1-alpha)/2), 
  )

# plot(lm_obs)

# let's see the fit
ggplot(data = fit_obs_noR0[[which(alphas == 0.95)]], aes(x = vax_cov)) +
  geom_line(aes(y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) + 
  geom_point(data = obs_df_small, aes(y = true_final_size)) + 
  ggtitle("without R0 predictor") + 
  labs(x = "vaccination coverage", y = "final size") + 
  theme_bw()

# let's look at how well our model predicts the true observed distribution
bind_rows(fit_obs_noR0) %>%
  filter(vax_cov %in% c(0.3, 0.5), alpha != 0) %>%
  melt(c("vax_cov", "alpha")) %>%
  filter(variable != "fit") %>%
  mutate(quantile = ifelse(variable == "lwr", (1-alpha)/2, 1-(1-alpha)/2), 
         scenario_id = ifelse(vax_cov == 0.3, "S1", "S2")) %>%
  ggplot(aes(x = value, y = quantile)) + 
  geom_line(aes(color = "predicted")) + 
  geom_line(data = t_small$true_sims %>% filter(scenario_id %in% c("S1", "S2")) %>%
              reframe(quantile = quantiles, 
                     value = quantile(true_final_size, quantiles), .by = c("scenario_id")), 
             aes(color = "true")) +
  facet_wrap(vars(scenario_id)) + 
  labs(x = "final size", title = "final size distribution across locations, R0 not included in the model") + 
  scale_color_manual(values = c("black", "red")) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

bind_rows(fit_obs_noR0) %>%
  filter(vax_cov %in% c(0.3, 0.5), alpha != 0) %>%
  left_join(t_small$true_sims %>% filter(scenario_id %in% c("S1", "S2")), relationship = "many-to-many") %>%
  mutate(cov = ifelse(true_final_size <= upr & true_final_size >= lwr, 1, 0)) %>%
  summarize(cov = sum(cov)/n(), .by = c("scenario_id", "alpha")) %>%
  ggplot(aes(x = alpha, y = cov)) + 
  geom_abline(linetype = "dashed") + 
  geom_line() + 
  ggtitle("Model to prediction observations without R0 vs. actual observations") +
  facet_wrap(vars(scenario_id)) +
  theme_bw()

### FIT OBSERVATIONS INSTEAD OF ERRORS (WITH R0 COVARIATE) ---------------------
ggplot(data = obs_df_small, aes(x = vax_cov, y = location_R0)) + 
  geom_point()

ggplot(data = obs_df_small, aes(x = vax_cov, y = true_final_size, color = location_R0)) + 
  geom_point() +
  geom_line(data = t_small$true_sims, aes(group = location_id), alpha = 0.4) +
  scale_color_viridis_c() + 
  theme_bw()

bx_wR0 = boxcox(true_final_size ~ vax_cov + location_R0, 
                lambda = seq(-3, 3, 1/10), data = obs_df_small, plotit = TRUE)
lambda_wR0 = with(bx_wR0, x[which.max(y)])
lm_obs_wR0 = lm(boxcox_transform(true_final_size, lambda_wR0) ~ vax_cov + location_R0, 
                data = obs_df_small) 

ggplot(data = obs_df_small, aes(x = vax_cov, y = boxcox_transform(true_final_size, lambda_wR0), color = location_R0)) + 
  geom_point() +
  geom_line(data = t_small$true_sims, aes(group = location_id), alpha = 0.4) +
  scale_color_viridis_c() + 
  theme_bw()

# plot(lm_obs_wR0)

# get prediction intervals
pred_mat_wR0 = expand.grid(vax_cov = new_vax_cov, 
                       location_id = unique(obs_df_small$location_id)) %>%
  left_join(unique(obs_df_small %>% dplyr::select(location_id, location_R0)))

fit_obs_wR0 <- vector("list", length(alphas))
for(i in 1:length(alphas)){
  # if(alphas[i] == 0.5){browser()}
  fit_obs_wR0[[i]] <- predict(lm_obs_wR0, newdata = pred_mat_wR0, 
                       level = alphas[i], interval = "prediction") %>%
    as.data.frame() %>% 
    bind_cols(pred_mat_wR0) %>%
    mutate(
      alpha = alphas[i],
      fit = ifelse(is.na(inverse_boxcox(fit, lambda_wR0)), 0, inverse_boxcox(fit, lambda_wR0)),
      lwr = ifelse(is.na(inverse_boxcox(lwr, lambda_wR0)), 0, inverse_boxcox(lwr, lambda_wR0)), 
      upr = ifelse(is.na(inverse_boxcox(upr, lambda_wR0)), 0, inverse_boxcox(upr, lambda_wR0)))
}
fit_obs_wR0_long = bind_rows(fit_obs_wR0) %>% 
  reshape2::melt(c("vax_cov", "location_id", "location_R0", "fit", "alpha")) %>%
  mutate(quantile = ifelse(variable == "lwr", (1-alpha)/2, 1-(1-alpha)/2), 
         scenario_id = ifelse(vax_cov == 0.3, "S1", ifelse(vax_cov == 0.5, "S2", "E")))

# plot relationship vs. true relationship for each location
bind_rows(fit_obs_wR0) %>% filter(alpha == 0.95) %>% 
  mutate(lwr = ifelse(is.na(lwr), 0, lwr)) %>%
  ggplot(aes(x = vax_cov)) + 
  geom_line(aes(y = fit), color = "darkgray") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(data = t_small$true_sims, aes(y = true_final_size, color = location_R0)) +
  facet_wrap(vars(paste0("R0: ", round(location_R0, 2), ", location ", location_id))) + 
  labs(x = "vaccination coverage", y = "final size") + 
  scale_color_viridis_c() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = "none")

# so we have an estimate of the observation for each location, along the 
# entire scenario axis

## test how well these predictions capture the true observations (across locations)
# generate a predicted distribution of observations across locations by drawing
# an equal number of samples from each predicted distribution
fit_obs_wR0_samp = fit_obs_wR0_long %>%
  dplyr::select(-variable, -alpha, - fit) %>%
  unique() %>%
  reframe(obs = get_samps(quantile, value, 1e4),
          draw_id = 1:1e4, .by = c("vax_cov", "location_id")
  ) %>%
  reframe(quantile = quantiles, 
          value = quantile(obs, quantiles), .by = c("vax_cov"))

fit_obs_wR0_samp %>% filter(quantile %in% c(0.05, 0.5, 0.95)) %>% 
  mutate(quantile = paste0("Q", quantile*100)) %>%
  dcast(vax_cov ~ quantile) %>%
  ggplot(aes(x = vax_cov)) + 
  geom_line(aes(y = Q50)) +
  geom_ribbon(aes(ymin = Q5, ymax = Q95), alpha = 0.2) +
  geom_point(data = t_small$true_sims %>% filter(scenario_id == "T"), aes(y = true_final_size)) +
  ggtitle("with R0 predictor") +
  labs(x = "vaccination coverage", y = "final size") + 
  theme_bw() 

cov_obs_wR0 = bind_rows(fit_obs_wR0_samp) %>%
  filter(vax_cov %in% c(0.3, 0.5), quantile != 0.5) %>%
  mutate(alpha = round(ifelse(quantile < 0.5, 1-2*quantile, 1-2*(1-quantile)), 4), 
         range = ifelse(quantile < 0.5, "lwr", "upr"))  %>%
  dcast(vax_cov + alpha ~ range, value.var = "value") %>%
  left_join(t_small$true_sims %>% filter(scenario_id %in% c("S1", "S2")), relationship = "many-to-many") %>%
  mutate(cov = ifelse(true_final_size <= upr & true_final_size >= lwr, 1, 0)) %>%
  summarize(cov = sum(cov)/n(), .by = c("scenario_id", "alpha"))

bind_rows(cov_obs_noR0 %>% mutate(model = "without R0 predictor"), 
          cov_obs_wR0 %>% mutate(model = "with R0 predictor")) %>%
  ggplot(aes(x = alpha, y = cov)) + 
  geom_abline(linetype = "dashed") + 
  geom_line(aes(color = model)) + 
  ggtitle("Coverage of estimated observations and actual observations") + 
  facet_wrap(vars(scenario_id), ncol = 1) + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid = element_blank())


# but this is the distribution of observations across all locations, in the case
# of the model that includes R0, we can also assess how well the model captured
# the observation at each location
fit_obs_wR0_long %>%
  filter(scenario_id == "S1") %>%
  ggplot(aes(x = value, y = quantile)) + 
  geom_line() + 
  geom_vline(data = t_small$true_sims %>% filter(scenario_id == "S1"), 
             aes(xintercept = true_final_size), linetype = "dashed") +
  facet_wrap(vars(paste0("R0: ", round(location_R0, 2), ", location ", location_id)), scales = "free") + 
  labs(x = "final size", title = "Model predicted final size vs. actual final size, scenario 1") + 
  theme_bw() + 
  theme(legend.position = "none", 
  panel.grid = element_blank())

fit_obs_wR0_long %>%
  filter(scenario_id == "S2") %>%
  ggplot(aes(x = value, y = quantile)) + 
  geom_line() + 
  geom_vline(data = t_small$true_sims %>% filter(scenario_id == "S2"), 
             aes(xintercept = true_final_size), linetype = "dashed") +
  facet_wrap(vars(paste0("R0: ", round(location_R0, 2), ", location ", location_id)), scales = "free") + 
  labs(x = "final size", title = "Model predicted final size vs. actual final size, scenario 2") + 
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid = element_blank())

fit_obs_wR0_long %>%
  left_join( t_small$true_sims) %>%
  filter(scenario_id %in% c("S1", "S2")) %>%
  ggplot(aes(x = location_R0, y = fit - true_final_size)) + 
  geom_point() + 
  facet_wrap(vars(scenario_id)) + 
  theme_bw()

#### CALCULATE MODEL ERROR USING ESTIMATED OBSERVATIONS ------------------------
# model error = projection - estimated error distribution

## first estimated observations without covariates in the model
fit_obs_noR0_samp = fit_obs_noR0_long %>%
  filter(paste0(quantile, variable) != "0.5lwr") %>% # remove duplicates
  reframe(est_obs_samp = get_samps(quantile, value, 1e4),
          draw_id = 1:1e4, .by = c("vax_cov")
  )

# join samples with model projections for each scenario and calculate error
fit_obs_noR0_errors = fit_obs_noR0_samp %>% 
  filter(vax_cov %in% c(0.3, 0.5)) %>%
  mutate(scenario_id = ifelse(vax_cov == 0.3, "S1", ifelse(vax_cov == 0.5, "S2", NA))) %>%
  left_join(t_small$model_sims %>% filter(scenario_id %in% c("S1", "S2")), 
            relationship = "many-to-many") %>%
  mutate(est_error = final_size - est_obs_samp)

# summarize into a distribution of errors across locations (1 per scenario and model)
fit_obs_noR0_errorints = fit_obs_noR0_errors %>% 
  reframe(quantile = quantiles, 
          est_error = quantile(est_error, quantiles), .by = c("model_id", "scenario_id"))

## next estimated observations with covariates
fit_obs_wR0_samp = fit_obs_wR0_long %>%
  filter(paste0(quantile, variable) != "0.5lwr") %>% # remove duplicates
  reframe(est_obs_samp = get_samps(quantile, value, 1e4),
          draw_id = 1:1e4, .by = c("vax_cov", "location_id")
  )

# join samples with model projections for each scenario and calculate error
fit_obs_wR0_errors = fit_obs_wR0_samp %>% 
  filter(vax_cov %in% c(0.3, 0.5)) %>%
  mutate(scenario_id = ifelse(vax_cov == 0.3, "S1", ifelse(vax_cov == 0.5, "S2", NA))) %>%
  left_join(t_small$model_sims %>% filter(scenario_id %in% c("S1", "S2")), 
            relationship = "many-to-many") %>%
  mutate(est_error = final_size - est_obs_samp)

# summarize into a distribution of errors across locations (1 per scenario and model)
fit_obs_wR0_errorints = fit_obs_wR0_errors %>% 
  reframe(quantile = quantiles, 
          est_error = quantile(est_error, quantiles), .by = c("model_id", "scenario_id"))

# plot some outcomes
all_preds = fit_obs_noR0_errorints %>%
  mutate(method = "fit observations without R0") %>%
  bind_rows(
    fit_obs_wR0_errorints %>%
      mutate(method = "fit observations with R0")
  ) %>%
  bind_rows(
    pred_intervals_shrink_small %>%
      filter(vax_cov %in% c(0.3, 0.5)) %>% 
      melt(c("model_id", "vax_cov"), variable.name = "quantile", value.name = "est_error") %>%
      mutate(scenario_id = ifelse(vax_cov == 0.3, "S1", "S2"), 
             quantile = as.double(substr(quantile, 2, nchar(as.character(quantile))))/100,
             method = "estimate error") 
  ) %>%
  bind_rows(
    t_small$errors %>%
      reframe(quantile = quantiles, 
              est_error = quantile(error, quantiles), .by = c("model_id", "scenario_id", "vax_cov")) %>%
      mutate(method = "true error")
  )

ggplot(data = all_preds %>% filter(scenario_id %in% c("S1")), 
       aes(x = est_error, y = quantile, color = method)) + 
  geom_path() + 
  facet_wrap(vars(model_id), scales = "free") + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(3, "Dark2"), "black")) + 
  theme_bw() + 
  theme(legend.position = "bottom")

