library(dplyr)
library(ggplot2)
library(reshape2)

# let's try to construct some examples where error functions are linear (in 
# order to avoid fitting complications)

n_mods <- 5
n_locs <- 25

set.seed(100)
# draw model-specific bias for slope and intercept terms
slope_bias <- c(rnorm(n_mods, 0, 0), 0) # 0 bias for truth
names(slope_bias) <- c(paste0("M", 1:n_mods), "T")
intercept_bias <- c(rnorm(n_mods, 0, 15), 0)
names(intercept_bias) <- c(paste0("M", 1:n_mods), "T")

# setup scenario values
scenario1 <- 10
scenario2 <- 20

# draw a true scenario value, slope and intercept for each location
true_scenario_vals <- runif(n_locs, scenario1, scenario2) # one true value per location
true_slope <- runif(n_locs, -5, 0)
true_intercept <- runif(n_locs, 60, 100)

# calculate model values for each
all_preds <- expand.grid(model_id = c(paste0("M", 1:n_mods), "T"), 
                           location_id = 1:n_locs, 
                           scenario_id = c("S1", "S2", "T"))

extra_scenarios <- expand.grid(model_id = c(paste0("M", 1:n_mods), "T"), 
                               location_id = 1:n_locs, 
                               scenario_value = seq(scenario1+1, scenario2-1, length.out = 20), 
                               scenario_id = "E")

all_preds <- all_preds %>%
  mutate(scenario_value = ifelse(scenario_id == "T", true_scenario_vals[location_id], 
                                 ifelse(scenario_id == "S1", scenario1, scenario2))) %>%
  bind_rows(extra_scenarios) %>%
  mutate(model_slope = true_slope[location_id] + slope_bias[model_id], 
         model_intercept = true_intercept[location_id] + intercept_bias[model_id], 
         value = model_intercept + model_slope*scenario_value)

# plot projections
ggplot(data = all_preds %>% filter(scenario_id != "E"), aes(x = scenario_value, y = value, group = model_id, color = model_id)) + 
  geom_point() + 
  geom_line(alpha = 0.5) + 
  facet_wrap(vars(location_id)) + 
  scale_color_manual(values = c(rep("black", n_mods), "red")) + 
  theme_bw()


# calculate errors
true_preds <- all_preds %>% filter(model_id == "T") %>%
  rename(true_value = value) %>%
  dplyr::select(location_id, scenario_id, scenario_value, true_value)
model_preds <-  all_preds %>% filter(model_id != "T")

errors <- left_join(model_preds, true_preds) %>%
  mutate(error = value - true_value, 
         relerror = (value - true_value)/true_value)

# plot errors
# note, the absolute errors are linear but the relative errors are not. 
# let's see if this causes problems. 
ggplot(data = errors %>% filter(scenario_id != "E"), 
       aes(x = scenario_value, y = error, color = scenario_id, shape = scenario_id)) + 
  geom_line(data = errors, aes(group = location_id), alpha = 0.4, color = "darkgray") + 
  geom_point() +  
  facet_wrap(vars(model_id)) + 
  scale_color_manual(values = c("darkgray", "darkgray", "black", "black")) + 
  scale_shape_manual(values = c(8, 8, 19, 19)) +
  theme_bw()

ggplot(data = errors %>% filter(scenario_id != "E"), 
       aes(x = scenario_value, y = relerror, color = scenario_id, shape = scenario_id)) + 
  geom_line(data = errors, aes(group = location_id), alpha = 0.4, color = "darkgray") + 
  geom_point() +  
  facet_wrap(vars(model_id)) + 
  scale_color_manual(values = c("darkgray", "darkgray", "black", "black")) + 
  scale_shape_manual(values = c(8, 8, 19, 19)) +
  theme_bw()

ggplot(data = errors %>% mutate(intercept_bias = intercept_bias[model_id]), 
       aes(x = intercept_bias, y = error)) + 
  geom_point()

# okay so this works as expected

# now let's add differences in scale across locations
set.seed(200)
location_scale = runif(n_locs, 0.5, 1.5)
all_preds_scaledlocation <- all_preds %>%
  mutate(value = value * location_scale[location_id])

ggplot(data = all_preds_scaledlocation %>% filter(scenario_id != "E"), 
       aes(x = scenario_value, y = value, group = model_id, color = model_id)) + 
  geom_point(aes(shape = model_id)) + 
  geom_line(alpha = 0.5) + 
  facet_wrap(vars(location_id)) + 
  scale_color_manual(values = c(rep("black", n_mods), "red")) + 
  theme_bw()

# calculate errors
true_preds_scaledlocation <- all_preds_scaledlocation %>% filter(model_id == "T") %>%
  rename(true_value = value) %>%
  dplyr::select(location_id, scenario_id, scenario_value, true_value)
model_preds_scaledlocation <-  all_preds_scaledlocation %>% filter(model_id != "T")

errors_scaledlocation <- left_join(model_preds_scaledlocation, true_preds_scaledlocation) %>%
  mutate(error = value - true_value, 
         relerror = (value - true_value)/true_value)

# estimate error distribution
preds = vector("list", length(unique(errors_scaledlocation$model_id)))
newdata = data.frame(scenario_value = seq(scenario1, scenario2, 0.1))
pred_level = c(seq(0.2, 0.9, 0.1), 0.95, 0.99)
for(i in 1:length(unique(errors_scaledlocation$model_id))){
  errors_scaledlocation_sub = errors_scaledlocation %>%
    filter(scenario_id == "T", model_id == paste0("M",i))
  fit = lm(error ~ scenario_value, data = errors_scaledlocation_sub)
  preds_tmp = list()
  for(j in 1:length(pred_level)){
    preds_tmp[[j]] =  as.data.frame(predict(fit, newdata, interval = "predict", level = pred_level[j])) %>%
      mutate(scenario_value = newdata$scenario_value, model_id = paste0("M",i)) %>%
      melt(c("scenario_value", "model_id", "fit")) %>%
      mutate(quantile = ifelse(variable == "lwr", paste0("Q", (0.5-pred_level[j]/2)*1000), paste0("Q", (0.5+pred_level[j]/2)*1000)))
  }
  preds[[i]] = bind_rows(preds_tmp)
}

preds = bind_rows(preds) %>% 
  dcast(scenario_value + model_id + fit ~ quantile, value.var = "value")

# so if we don't correct for variation in the scale across locations, we 
# will (probably) get the mean right, but we will overestimate variance
ggplot(data = errors_scaledlocation %>% filter(scenario_id != "E", model_id !="T"), 
       aes(x = scenario_value)) + 
  geom_line(data = errors_scaledlocation, aes(y = error, color = scenario_id, group = location_id), alpha = 0.4, color = "darkgray") + 
  geom_point(aes(y = error, color = scenario_id, shape = scenario_id)) +  
  geom_ribbon(data = preds, 
              aes(ymin = Q25, ymax = Q975), alpha = 0.2, fill = "blue") +
  geom_ribbon(data = preds, 
              aes(ymin = Q250, ymax = Q750), alpha = 0.3, fill = "blue") +
  geom_line(data = preds, aes(y = fit), color = 'blue', linewidth = 1) + 
  geom_hline(data = data.frame(int = intercept_bias[-6], 
                               model_id = names(intercept_bias)[-6]), 
             aes(yintercept = int), color = "red", linewidth = 1, linetype = 'dashed') + 
  facet_wrap(vars(model_id)) + 
  scale_color_manual(values = c("darkgray", "darkgray", "black", "black")) + 
  scale_shape_manual(values = c(8, 8, 19, 19)) +
  theme_bw()

# estimated distribution of errors
preds %>% 
  filter(scenario_value %in% c(scenario1, scenario2)) %>%
  melt(c("scenario_value", "model_id", "fit")) %>%
  mutate(quantile = as.integer(substr(variable, 2, nchar(as.character(variable))))/1000) %>%
  ggplot(aes(x = value, y = quantile, color = model_id)) + 
  geom_point() + 
  geom_line() + 
  geom_vline(data = data.frame(int = intercept_bias[-6], 
                               model_id = names(intercept_bias)[-6]), 
             aes(xintercept = int, color = model_id), linetype = 'dashed') +
  facet_wrap(vars(scenario_value), ncol = 1, labeller = label_both) + 
  theme_bw()
  
  
# get coverage
cov = preds %>%
  filter(scenario_value %in% c(scenario1, scenario2)) %>%
  dplyr::select(-fit) %>%
  melt(c("scenario_value", "model_id")) %>%
  mutate(quantile = as.integer(substr(variable, 2, nchar(as.character(variable))))/1000) %>%
  mutate(alpha = round(ifelse(quantile < 0.5, 1-2*quantile, 1-2*(1-quantile)),3), 
         bound = ifelse(quantile < 0.5, "lwr", "upr")) %>%
  reshape2::dcast(model_id + scenario_value + alpha ~ bound, value.var = "value") %>%
  left_join(
    errors_scaledlocation %>%
      filter(scenario_value %in% c(scenario1, scenario2)) %>%
      mutate(obs = error) %>%
      dplyr::select(model_id, location_id, scenario_value, obs),
    relationship = "many-to-many", by = join_by(model_id, scenario_value)
  ) %>%
  mutate(cov = ifelse(obs <= upr & obs >= lwr, 1, 0)) %>%
  summarize(cov = sum(cov)/n(), .by = c("scenario_value", "model_id", "alpha"))

ggplot(data = cov, aes(x = alpha, y = cov, color = model_id)) + 
  geom_abline() + 
  geom_point() + 
  geom_line() + 
  facet_wrap(vars(scenario_value), labeller = label_both, ncol = 1) +
  theme_bw()

# get pdf of estimated error distribution
get_samps = function(value, quantile, n_samps = 1e5, 
                     quantiles = c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)){
  s = approx(quantile, value, runif(n_samps))[[2]]
  return(data.frame(sample = s, 
                    sample_id = 1:length(s)))
  # return(data.frame(quantile = quantiles, 
                    # value = quantile(s, quantiles, na.rm = TRUE)))
}

calibration_error = preds %>%
  filter(scenario_value %in% c(scenario1, scenario2)) %>%
  mutate(scenario_id = ifelse(scenario_value == scenario1, "S1", "S2")) %>%
  dplyr::select(-fit) %>%
  dplyr::select(-scenario_value) %>%
  melt(c("scenario_id", "model_id")) %>%
  mutate(quantile = as.integer(substr(variable, 2, nchar(as.character(variable))))/1000) %>%
  reframe(get_samps(value, quantile), .by = c("scenario_id", "model_id")) %>%
  rename(calibration_error = sample)

scenario_error = calibration_error %>%
  left_join(model_preds_scaledlocation %>% 
              filter(scenario_id %in% c("S1", "S2")) %>%
              dplyr::select(model_id, location_id, scenario_id, value) %>%
              left_join(true_preds_scaledlocation %>% filter(scenario_id == c("T")) %>% dplyr::select(-scenario_value, -scenario_id)) %>%
              mutate(observed_error = value - true_value) %>% 
              dplyr::select(model_id, location_id, scenario_id, observed_error), 
            relationship = "many-to-many", by = join_by(model_id, scenario_id)) %>%
  mutate(scenario_error = observed_error - calibration_error)

ggplot(data = scenario_error) + 
  geom_density(aes(x = scenario_error, color = as.factor(model_id), linetype = as.factor(scenario_id))) + 
  facet_wrap(vars(location_id)) + 
  theme_bw()
  
scenario_error %>% 
  filter(!is.na(scenario_error)) %>%
  summarize(mean = mean(scenario_error), 
            lwr = quantile(scenario_error, 0.05), 
            upr = quantile(scenario_error, 0.95), 
            .by = c("scenario_id", "location_id", "model_id")) %>%
  mutate(true_slope = true_slope[location_id]) %>%
  ggplot(aes(x = true_slope, color = as.factor(location_id))) + 
  geom_hline(aes(yintercept = 0)) + 
  geom_point(aes(y = mean), size = 3) + 
  geom_segment(aes(xend = true_slope, y = lwr, yend = upr), size = 1) + 
  geom_text(aes(y = mean, label = location_id), size = 2.5, color = "white") + 
  facet_grid(cols = vars(model_id), rows = vars(scenario_id)) + 
  labs(x = "true slope of relationship in location", y = "scenario error") +
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid.minor = element_blank())

ggsave("scenario_error_estimate.pdf", width = 10, height = 6)

# now percent scenario error
scenario_error %>% 
  filter(!is.na(scenario_error)) %>%
  mutate(total_error = abs(calibration_error) + abs(scenario_error), 
         scenario_error_pct = abs(scenario_error)/total_error) %>%
  summarize(mean = mean(scenario_error_pct), 
            lwr = quantile(scenario_error_pct, 0.05), 
            upr = quantile(scenario_error_pct, 0.95), 
            .by = c("scenario_id", "location_id", "model_id")) %>%
  mutate(true_slope = true_slope[location_id]) %>%
  left_join(true_preds_scaledlocation %>% filter(scenario_id == "T") %>% select(-true_value, -scenario_id)) %>%
  ggplot(aes(x = scenario_value, color = true_slope)) + #, color = as.factor(location_id)
  geom_hline(aes(yintercept = 0)) + 
  geom_point(aes(y = mean), size = 3) + 
  geom_segment(aes(xend = scenario_value, y = lwr, yend = upr), size = 1) + 
  geom_text(aes(y = mean, label = location_id), size = 2.5, color = "white") + 
  facet_grid(cols = vars(model_id), rows = vars(scenario_id)) + 
  labs(x = "realized scenario value", y = "absolute scenario error percent") + 
  scale_color_viridis_c() +
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid.minor = element_blank())
ggsave("percent_scenario_error_estimate.pdf", width = 10, height = 6)

