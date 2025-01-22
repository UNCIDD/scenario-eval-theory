library(finalsize)
library(dplyr)
library(mgcv)
library(ggplot2)
library(gratia)

n_locations = 25
n_models = 6

set.seed(100)
vax_cov_S1 = 0.3
vax_cov_S2 = 0.5
vax_cov_T = runif(n_locations, vax_cov_S1, vax_cov_S2)

R0_T = runif(n_locations, 2, 3.5)
model_bias_R0 = c(rnorm(n_models, 0, 0.05), 0) # 0 for true model
names(model_bias_R0) = c(paste0("M", 1:n_models), "T")

sims <- expand.grid(model_id = c(paste0("M", 1:n_models), "T"),
                    location_id = 1:n_locations, 
                    scenario_id = c("S1", "T", "S2"))
sims$vax_cov = with(sims, ifelse(scenario_id == "T", vax_cov_T[location_id], 
                                 ifelse(scenario_id == "S1", vax_cov_S1, vax_cov_S2)))

sims_full_relationship <- expand.grid(model_id =  c(paste0("M", 1:n_models), "T"), 
                                      location_id = 1:n_locations, 
                                      vax_cov = seq(vax_cov_S1, vax_cov_S2, length.out = 20)) %>%
  mutate(scenario_id = "E")

sims <- bind_rows(sims, sims_full_relationship)

sims$R0 = R0_T[sims$location_id] + model_bias_R0[sims$model_id]

sims$final_size = NA

# final size variables
susc_immunised <- cbind(1,0)
colnames(susc_immunised) <- c("novax", "vax")

for(i in 1:nrow(sims)){
  p_susc_immunised <- cbind(
    susceptible = 1 - sims$vax_cov[i],
    immunised = sims$vax_cov[i]
  )
  fs <- final_size(
    r0 = sims$R0[i],
    contact_matrix = matrix(1.0)/1E3,
    demography_vector = 1E3,
    susceptibility = susc_immunised,
    p_susceptibility = p_susc_immunised
  )
  sims$final_size[i] <- fs[1,4]
}

ggplot(data = sims, aes(x = vax_cov, y = final_size, color = as.factor(R0))) +
  geom_line() +
  theme(legend.position = "none")

model_sims = sims %>% 
  filter(model_id != "T")
true_sims = sims %>%
  filter(model_id == "T") %>%
  rename(true_final_size = final_size) %>%
  select(location_id, scenario_id, vax_cov, true_final_size)

ggplot(data = sims %>% filter(scenario_id != "E"), 
       aes(x = vax_cov, y = final_size, color = model_id)) + 
  geom_point() + 
  geom_line(alpha = 0.5) + 
  facet_wrap(vars(location_id)) + 
  scale_color_manual(values = c(rep("black", n_models), "red")) + 
  theme_bw()

errors = model_sims %>% 
  left_join(true_sims) %>%
  mutate(error = final_size - true_final_size, 
         relerror = (final_size - true_final_size)/true_final_size)

new_vax_cov = data.frame(vax_cov = seq(vax_cov_S1, vax_cov_S2, length.out = 100))
lin_fits <- list()
gam_fits <- list()
alpha_vals <- c(seq(0.1, 0.9, 0.1), 0.95, 0.99)
quantiles = c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)
for(i in 1:n_models){
  dat_filt = errors %>%
    filter(model_id == paste0("M", i))
  lin_fits[[i]] = lm(error ~ vax_cov, data = dat_filt)
  alpha_vals = c(seq(0.1, 0.9, 0.1), 0.95, 0.99)
  preds_tmp <- vector("list", length(alpha_vals))
  for(j in 1:length(alpha_vals)){
    preds_tmp[[j]] = predict(lin_fits[[i]], newdata = new_vax_cov,
                             interval = "prediction", level = alpha_vals[[j]])
    preds_tmp[[j]] = data.frame(x = new_vax_cov$vax_cov, 
                                y = preds_tmp[[j]][,1], 
                                alpha = alpha_vals[j],
                                lwr = preds_tmp[[j]][,2],
                                upr = preds_tmp[[j]][,3])
  }
  lin_fits[[i]] = bind_rows(preds_tmp)
  # gam 
  gf = gam(error ~ s(vax_cov), data = dat_filt)
  ps <- posterior_samples(gf, n = 10000, data = new_vax_cov, seed = 24,
                          unconditional = TRUE) |>
    left_join(new_vax_cov %>% mutate(.row = row_number()), by = join_by(.row == .row))
  gam_fits[[i]] <- ps %>%
    reframe(quantile = quantiles, 
            value = quantile(.response, quantiles), .by = "vax_cov")
}



ggplot() + 
  geom_hline(yintercept = 0) + 
  geom_line(data = errors, 
            aes(x = vax_cov, y = error, group = location_id), color = "gray", alpha = 0.6) + 
  geom_line(data = errors %>% filter(scenario_id != "T") %>%
              group_by(model_id, vax_cov) %>% summarize(mean_error = mean(error)), 
            aes(x = vax_cov, y = mean_error), color = "black", size = 1.5)  +
  geom_point(data = errors %>% filter(scenario_id != "E"), 
             aes(x = vax_cov, y = error, color = scenario_id, shape = scenario_id)) + 
  geom_ribbon(data = bind_rows(lin_fits, .id = "model_id") %>%
                mutate(model_id = paste0("M", model_id)) %>%
                filter(alpha == 0.95),
              aes(x = x, ymin = lwr, ymax = upr), fill = "blue", alpha = 0.2) +
  geom_line(data = bind_rows(lin_fits, .id = "model_id") %>%
              mutate(model_id = paste0("M", model_id)) %>%
              filter(alpha == 0.5),
            aes(x = x, y = y), color = "blue", linewidth = 1, linetype = "longdash") +
  facet_wrap(vars(model_id)) + 
  labs(x = "vaccination coverage", y = "error") +
  scale_color_manual(values = c("darkgray", "darkgray", "black")) + 
  scale_shape_manual(values = c(8, 8, 19)) +
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

ggsave("linear-fits-example.pdf", width = 6, height = 4)

ggplot() + 
  geom_hline(yintercept = 0) + 
  geom_line(data = errors, 
            aes(x = vax_cov, y = error, group = location_id), color = "gray", alpha = 0.6) + 
  geom_point(data = errors %>% filter(scenario_id != "E"), 
             aes(x = vax_cov, y = error, color = scenario_id, shape = scenario_id)) + 
  geom_line(data = errors %>% filter(scenario_id != "T") %>%
              group_by(model_id, vax_cov) %>% summarize(mean_error = mean(error)), 
            aes(x = vax_cov, y = mean_error), color = "black", size = 1.5)  +
  geom_ribbon(data = bind_rows(gam_fits, .id = "model_id") %>%
                mutate(model_id = paste0("M", model_id)) %>% 
                filter(quantile %in% c(0.025, 0.975)) %>%
                mutate(quantile = paste0("Q", quantile*1000)) %>%
                reshape2::dcast(model_id + vax_cov ~ quantile), 
              aes(x = vax_cov, ymin = Q25, ymax = Q975), fill = "red", alpha = 0.2) + 
  geom_line(data = bind_rows(gam_fits, .id = "model_id") %>%
              mutate(model_id = paste0("M", model_id)) %>% 
              filter(round(quantile,3) == round(0.5,3)), 
            aes(x = vax_cov, y = value), color = "red", linewidth = 1, linetype = "dashed") + 
  facet_wrap(vars(model_id)) + 
  labs(x = "vaccination coverage", y = "relative error") +
  scale_color_manual(values = c("darkgray", "darkgray", "black")) + 
  scale_shape_manual(values = c(8, 8, 19)) +
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

ggsave("gam-fits-example.pdf", width = 6, height = 4)

# check mean true errors vs. model estimates
errors %>%  
  filter(scenario_id %in% c("S1", "S2")) %>%
  summarize(mean_error = mean(error), .by = c("model_id")) %>%
  left_join(bind_rows(lin_fits, .id = "model_id") %>%
              mutate(model_id = paste0("M", model_id)) %>%
              filter(alpha == 0.5, x %in% c(vax_cov_S1, vax_cov_S2)) %>%
              mutate(scenario_id = ifelse(x == vax_cov_S1, "S1", "S2")) %>% 
              select(model_id, scenario_id, y) %>% 
              rename(pred_error_lin = y)) %>%
  left_join(bind_rows(gam_fits, .id = "model_id") %>%
              mutate(model_id = paste0("M", model_id)) %>%
              filter(quantile == 0.5, vax_cov %in% c(vax_cov_S1, vax_cov_S2)) %>%
              mutate(scenario_id = ifelse(vax_cov == vax_cov_S1, "S1", "S2")) %>%
              select(model_id, scenario_id, value) %>%
              rename(pred_error_gam = value)
  ) %>%
  reshape2::melt(c("model_id", "scenario_id", "mean_error")) %>%
  ggplot(aes(x = value, y = mean_error, color = variable)) + 
  geom_abline() + 
  geom_point(aes(shape = model_id)) +
  geom_line() + 
  facet_wrap(vars(scenario_id)) + 
  labs(x = "mean true error across locations", y = "predicted error") +
  scale_color_manual(values = c("blue", "red")) + 
  theme_bw()




cov <- bind_rows(lin_fits, .id = "model_id") %>%
  mutate(model_id = paste0("M", model_id)) %>%
  filter(x %in% c(vax_cov_S1,vax_cov_S2)) %>%
  mutate(scenario_id = ifelse(x == vax_cov_S1, "S1", "S2")) %>%
  select(model_id, scenario_id, alpha, lwr, upr) %>%
  left_join(
    errors %>%
      filter(scenario_id %in% c("S1", "S2")) %>%
      mutate(obs = relerror) %>%
      select(model_id, location_id, scenario_id, obs),
    relationship = "many-to-many"
  ) %>%
  mutate(cov = ifelse(obs <= upr & obs >= lwr, 1, 0))

cov_gam <- bind_rows(gam_fits, .id = "model_id") %>%
  mutate(model_id = paste0("M", model_id)) %>%
  filter(vax_cov %in% c(vax_cov_S1,vax_cov_S2)) %>%
  mutate(scenario_id = ifelse(vax_cov == vax_cov_S1, "S1", "S2")) %>%
  mutate(alpha = round(ifelse(quantile < 0.5, 1-quantile, quantile),3), 
         bound = ifelse(quantile < 0.5, "lwr", "upr")) %>%
  reshape2::dcast(model_id + vax_cov + scenario_id + alpha ~ bound, value.var = "value") %>%
  filter(alpha != 0.5) %>%
  select(model_id, scenario_id, alpha, lwr, upr) %>%
  left_join(
    errors %>%
      filter(scenario_id %in% c("S1", "S2")) %>%
      mutate(obs = relerror) %>%
      select(model_id, location_id, scenario_id, obs),
    relationship = "many-to-many"
  ) %>%
  mutate(cov = ifelse(obs <= upr & obs >= lwr, 1, 0))


cov %>%
  mutate(fitting_model = "linear") %>%
  bind_rows(cov_gam %>% mutate(fitting_model = "gam")) %>%
  summarize(cov = sum(cov)/n(), .by = c("fitting_model", "scenario_id", "alpha")) %>%
  ggplot(aes(x = alpha, y = cov, color = fitting_model)) + 
  geom_line(aes(group = interaction(scenario_id, fitting_model))) +
  geom_abline() + 
  facet_wrap(vars(scenario_id)) + 
  labs(x = "expected coverage", y = "actual coverage") + 
  scale_color_manual(values = c("red", "blue")) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid = element_blank())

ggsave("coverage.pdf", width = 6, height = 4)



# 
# 
# library(mgcv)
# set.seed(1)
# dat <- gamSim(eg = 1, n = 400, dist = "normal", scale = 2)
# fit_gam <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = dat)
# 
# beta_sim <- coef(fit_gam)
# V <- vcov(fit_gam)
# num_beta_vecs <- 10000
# Cv <- chol(V)
# 
# set.seed(1)
# nus <- rnorm(num_beta_vecs * length(beta_sim))
# beta_sims <- beta_sim + t(Cv) %*% matrix(nus, nrow = length(beta_sim), ncol = num_beta_vecs)
# 
# n_obs <- 100
# sim_idx <- sample.int(nrow(dat), size = n_obs, replace = TRUE)
# sim_dat <- dat[sim_idx, c("x0", "x1", "x2", "x3")]
# 
# covar_sim <- predict(fit_gam, newdata = sim_dat, type = "lpmatrix")
# linpred_sim <- covar_sim %*% beta_sims
# 
# invlink <- function(x) x
# exp_val_sim <- invlink(linpred_sim)
# 
# y_sim <- matrix(rnorm(n = prod(dim(exp_val_sim)), 
#                       mean = exp_val_sim, 
#                       sd = sqrt(summary(fit_gam)$scale)), 
#                 nrow = nrow(exp_val_sim), 
#                 ncol = ncol(exp_val_sim))