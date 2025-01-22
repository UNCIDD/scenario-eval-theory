library(dplyr)
library(ggplot2)

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
  select(location_id, scenario_id, scenario_value, true_value)
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
  geom_point() + 
  geom_line(alpha = 0.5) + 
  facet_wrap(vars(location_id)) + 
  scale_color_manual(values = c(rep("black", n_mods), "red")) + 
  theme_bw()

# calculate errors
true_preds_scaledlocation <- all_preds_scaledlocation %>% filter(model_id == "T") %>%
  rename(true_value = value) %>%
  select(location_id, scenario_id, scenario_value, true_value)
model_preds_scaledlocation <-  all_preds_scaledlocation %>% filter(model_id != "T")

errors_scaledlocation <- left_join(model_preds_scaledlocation, true_preds_scaledlocation) %>%
  mutate(error = value - true_value, 
         relerror = (value - true_value)/true_value)

# so if we don't correct for variation in the scale across locations, we 
# will (probably) get the mean right, but we will overestimate variance
ggplot(data = errors_scaledlocation %>% filter(scenario_id != "E"), 
       aes(x = scenario_value, y = error, color = scenario_id, shape = scenario_id)) + 
  geom_line(data = errors_scaledlocation, aes(group = location_id), alpha = 0.4, color = "darkgray") + 
  geom_point() +  
  geom_hline(data = data.frame(int = intercept_bias, 
                               model_id = names(intercept_bias)), 
             aes(yintercept = int), color = "blue") + 
  facet_wrap(vars(model_id)) + 
  scale_color_manual(values = c("darkgray", "darkgray", "black", "black")) + 
  scale_shape_manual(values = c(8, 8, 19, 19)) +
  theme_bw()

# so let's try to correct for variation in location scale

