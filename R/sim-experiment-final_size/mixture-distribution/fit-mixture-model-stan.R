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

summarize_predints <- function(fit, new_vax_cov){
  y_exclude = filter_posterior(as.data.frame(fit), "y_exclude\\[") %>%
    rename(exclude = value) %>%
    dplyr::select(-variable)
  y_nominal = filter_posterior(as.data.frame(fit), "y_new_nominal\\[") %>%
    left_join(y_exclude) %>%
    left_join(data.frame(index_id = 1:length(new_vax_cov), 
                         vax_cov = new_vax_cov)) %>%
    filter(exclude == 0) %>%
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
new_vax_cov = data.frame(vax_cov = seq(0.3, 0.5, length.out = 100))

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

write_rds(fit_errors, "R/sim-experiment-final_size/mixture-distribution/fit_shrinkage.rda")

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
  quant_fits = pred_intervals %>%
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
          
