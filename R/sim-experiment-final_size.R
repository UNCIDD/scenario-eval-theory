### simulation experiment to test scenario evaluation theory
### assume projections of final outbreak size across varying vaccination rates

#### SETUP ---------------------------------------------------------------------
library(finalsize)
library(dplyr)
library(reshape2)
library(ggplot2)
library(quantreg)
library(gamlss)
library(bestNormalize)
source("./R/final-size-functions.R")

#### GENERATE SIMULATIONS TO TEST ----------------------------------------------
set.seed(10)
seed_id = sample(1:1E6, reps)
n_models = 2

# quantiles of interest for summary error distribution
quantiles = c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)

# new vaccination coverage
new_vax_cov = data.frame(vax_cov = seq(0.3, 0.5, length.out = 100))

# run simulation - using 500 locations for now
t_large <- full_sim(n_locations = 500, n_models = n_models,
                    vax_cov_S1 = 0.3, vax_cov_S2 = 0.5, 
                    R0_lwr = 2, R0_upr = 3.25,
                    seed = seed_id, fit_outcomes = FALSE)

# note: scenario_id = T returns only observed errors (not true errors to test against)
error_df = t_large$errors %>% filter(scenario_id == "T") 

#### PLOT ERRORS FOR EACH MODEL ------------------------------------------------
ggplot(data = error_df, aes(x = vax_cov, y = error)) + 
  geom_point() + 
  facet_wrap(vars(model_id), nrow = 1, scales = "free") + 
  theme_bw()

#### TEST NORMALITY AND TRY SOME TRANSFORMATIONS -------------------------------
# choose model 2 to test (negative values and skewed)
dat_filt = error_df %>% filter(model_id == "M2")

# check qq plot
qqnorm(dat_filt$error)
qqline(dat_filt$error)

# try some transforms# testing some transforms
# yeo-johnson seems to work with neg values
yj_errors = yeojohnson(dat_filt$error)
ordnorm_errors = orderNorm(dat_filt$error)
# from https://www.listendata.com/2015/09/regression-transform-negative-values.html
# note: cube root and log only allows you to trasnform abs error (cannot recover signed error)
cuberoot_errors = sign(dat_filt$error)*abs(dat_filt$error)^(1/3)
log_errors = sign(dat_filt$error)*log(abs(dat_filt$error))
adjlog_errors = log(1 + dat_filt$error - min(dat_filt$error))
sinh_errors = asinh(dat_filt$error)

par(mfrow = c(2,3))
MASS::truehist(yj_errors$x.t, main = "yeo-johnson")
MASS::truehist(yj_errors$x.t, main = "ordered norm")
MASS::truehist(adjlog_errors, main = "adjusted log")
MASS::truehist(cuberoot_errors, main = "sign(x) abs(x)^(1/3)")
MASS::truehist(log_errors, main = "sign(x) log(abs(x))")
MASS::truehist(sinh_errors, main = "asinh")

dat_filt = dat_filt %>% 
  mutate(abs_log_trans = log(abs(error)))

dat_filt %>% 
  select(model_id, vax_cov, error, abs_log_trans) %>% 
  melt(c("model_id", "vax_cov")) %>% 
  ggplot(aes(x = vax_cov, y = value)) + 
  
  
#### TRY SOME FITTING FOR ONE MODEL --------------------------------------------
# start by using GAM with log transform on abs error
# note, we can only fit/predict absolute error
gam_fit <- mgcv::gam(abs_log_trans ~ s(vax_cov), data = dat_filt)
gam_pred <- predict(gam_fit, newdata = new_vax_cov, se.fit = TRUE)

gam_PIs <- get_gam_PIs(mod = gam_fit, xvals = new_vax_cov, invfn = exp)

gam_PIs %>% 
  filter(round(quantile,3) %in% c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
  mutate(quantile = paste0("Q", quantile*1000)) %>% 
  reshape2::dcast(vax_cov ~ quantile) %>% 
  ggplot(aes(x = vax_cov)) + 
  geom_ribbon(aes(ymin = Q25, ymax = Q975), fill = "blue", alpha = 0.2) +
  geom_ribbon(aes(ymin = Q250, ymax = Q750), fill = "blue", alpha = 0.2) +
  geom_line(data = data.frame(vax_cov = new_vax_cov, pred = exp(gam_pred$fit)), 
            aes(y = pred), color = "blue", size = 1) + 
  geom_point(data = dat_filt, aes(y = abs(error))) +
  labs(x = "vaccination coverage", y = "absolute error") + 
  theme_bw()


#### REPEAT FOR BOTH MODELS ----------------------------------------------------
all_mod_fits <- vector("list", n_models)
for(i in 1:n_models){
  dat_filt = error_df %>% 
    filter(model_id == paste0("M", i)) %>% 
    mutate(abs_log_trans = log(abs(error)))
  gam_fit <- mgcv::gam(abs_log_trans ~ s(vax_cov), data = dat_filt)
  gam_PIs <- get_gam_PIs(mod = gam_fit, xvals = new_vax_cov, invfn = exp)
  all_mod_fits[[i]] <- gam_PIs %>% 
    mutate(quantile = paste0("Q", quantile*1000)) %>% 
    reshape2::dcast(vax_cov ~ quantile) %>% 
    cbind(mean = exp(predict(gam_fit, newdata = new_vax_cov))) # using exp transform here
}

# plot predictions
bind_rows(all_mod_fits, .id = "model_id") %>% 
  mutate(model_id = paste0("M", model_id)) %>%
  ggplot(aes(x = vax_cov)) + 
  geom_ribbon(aes(ymin = Q25, ymax = Q975), fill = "blue", alpha = 0.2) +
  geom_ribbon(aes(ymin = Q250, ymax = Q750), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = mean), color = "blue", size = 1) + 
  geom_point(data = error_df, aes(y = abs(error))) +
  facet_wrap(vars(model_id), scales = "free") + 
  labs(x = "vaccination coverage", y = "absolute error") + 
  theme_bw()
ggsave("sim-exp_prediction-intervals.pdf", width = 10, height = 5)

# find coverage and plot
fits_cov = calculate_coverage(
  # do some reshaping to match expected format for quant_fits
  quant_fits = bind_rows(all_mod_fits, .id = "model_id") %>% 
    mutate(model_id = paste0("M", model_id)) %>%
    melt(c("model_id", "vax_cov")) %>% 
    filter(variable != "mean") %>%
    mutate(quantile = as.integer(gsub("Q", "", variable))/1000) %>%
    select(-variable), 
  error_df = t_large$errors %>% mutate(error = abs(error)), 
  vax_cov_S1 = 0.3, vax_cov_S2 = 0.5, 
  summarize_by = "model"
  )

scenario_labs = c("low vax scenario", "high vax scenario")
names(scenario_labs) = c("S1", "S2")

ggplot(data = fits_cov, aes(x = alpha, y = cov, group = model_id)) + 
  geom_line(color = "blue", size = 1, alpha = 0.5) + 
  geom_abline() + 
  facet_wrap(vars(scenario_id), labeller = labeller(scenario_id = scenario_labs)) + 
  labs(x = "expected coverage", y = "actual coverage") + 
  theme_bw()
ggsave("sim-exp_coverage.pdf", width = 10, height = 5)

