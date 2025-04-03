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
library(bayesplot)

source("./R/final-size-functions.R")

scenario_labs = c("low vax scenario", "high vax scenario")
names(scenario_labs) = c("S1", "S2")

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

error_df_sub = error_df %>% 
  filter(model_id == "M10") %>%
  mutate(error_boxcox = ifelse(error < 0, boxcox_transform(-error, lambda), boxcox_transform(error, lambda)))

# do box-cox ahead of time (for now)
# bx = boxcox(error ~ vax_cov, data = error_df_sub, plotit = FALSE)
# lambda = with(bx, x[which.max(y)])
# 
error_df_sub$error_boxcox = ifelse(error_df_sub < 0, boxcox_transform(-error_df_sub$error, lambda), boxcox_transform(error_df_sub$error, lambda))
# 
# ggplot(data = error_df_sub, aes(x = vax_cov, y = error_boxcox)) + 
#   geom_point()
# 
# # try with lm
# fit_lm <- lm(error_boxcox ~ vax_cov, error_df_sub)

#### FIT WITH STAN -------------------------------------------------------------
mixture_model <- stan_model("R/sim-experiment-final_size/mixture-distribution/mixture_model.stan")

fit_errors <- sampling(
  mixture_model,
  data = list(
    N = nrow(error_df_sub),
    x = error_df_sub$vax_cov,
    y = error_df_sub$error, 
    N_new = length(new_vax_cov), 
    x_new = new_vax_cov
  ),
  seed = 7, 
  iter = 5000,
  chain = 4, 
  cores = 4
)
beepr::beep()

# some summary/diagnostics
posterior <- as.array(fit_errors)
posterior_df <- as.data.frame(fit_errors)
np <- nuts_params(fit_errors)

mcmc_pairs(fit_errors, pars = c("alpha_neg", "beta_neg","lambda_neg", "sigma_neg","lambda_pos", "alpha_pos", "beta_pos", "sigma_pos", "p"), np = np)

pars <- bind_cols(extract(fit_errors, c("alpha_pos", "beta_pos", "lambda_pos", "alpha_neg", "beta_neg", "lambda_neg", "p"))) %>% 
  mutate(draw_id = seq_len(n()))

#### GET PREDICTION INTERVALS --------------------------------------------------
y_nominal = filter_posterior(posterior_df, "y_new_nominal\\[") %>%
  left_join(data.frame(index_id = 1:length(new_vax_cov), 
                       vax_cov = new_vax_cov)) %>%
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

# nominal scale
ggplot(data = y_nominal, aes(x = vax_cov)) + 
  geom_point(data = error_df_sub, aes(x = vax_cov, y = error)) + 
  geom_ribbon(aes(ymin = Q5, ymax = Q95), alpha = 0.2) + 
  geom_ribbon(aes(ymin = Q25, ymax = Q75), alpha = 0.2) + 
  geom_line(aes(y = Q50), size = 1)

# check coverage
cov = calculate_coverage(
  # do some reshaping to match expected format for quant_fits
  quant_fits = y_nominal %>%
    mutate(model_id = "M10") %>%
    melt(c("model_id", "vax_cov")) %>% 
    filter(variable != "mean") %>%
    mutate(quantile = as.integer(gsub("Q", "", variable))/100) %>%
    dplyr::select(-variable), 
  error_df = t_large$errors,  # %>% mutate(error = abs(error))
  vax_cov_S1 = 0.3, vax_cov_S2 = 0.5, 
  summarize_by = "model"
)

ggplot(data = cov, aes(x = alpha, y = cov, group = model_id)) + 
  geom_line(aes(color = model_id)) + 
  geom_abline(size = 1) + 
  facet_wrap(vars(scenario_id), labeller = labeller(scenario_id = scenario_labs), ncol = 1) + 
  labs(x = "expected coverage", y = "actual coverage") + 
  theme_bw() + 
  theme(legend.position = "none")

