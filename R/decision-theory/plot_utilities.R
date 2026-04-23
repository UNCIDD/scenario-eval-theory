library(dplyr)
library(reshape2)
library(ggplot2)
library(MASS) 

sim_out = readRDS("output/simulation/sim_out.rds")

source("R/simulation/0-parameters.R")
source("R/simulation/0-helper-functions.R")

linear_utility = function(y, k, c){
  return(y*k + c) 
}

k_val = 10
c_val = 5

samp_size = 250

true_utilities = sim_out$true_sims %>%
  filter(scenario_id %in% c("S1", "S2")) %>%
  mutate(utility = linear_utility(true_final_size, k_val, c_val)) %>%
  dcast(location_id ~ scenario_id, value.var = "utility")

true_final_sizes = sim_out$true_sims %>%
  filter(scenario_id %in% c("S1", "S2")) %>%
  dcast(location_id ~ scenario_id, value.var = "true_final_size")
 
# add some artificial uncertainty (for now) to the projections
set.seed(1928)
model_scenario_sd = expand.grid(
  model_id = paste0("M", 1:10) 
  # scenario_id = c("S1", "S2")
) %>%
  mutate(projection_sd = runif(n(), 0, 0.0005))

set.seed(1928)
projected_utilities = sim_out$model_sims %>%
  filter(scenario_id %in% c("S1", "S2")) %>%
  left_join(model_scenario_sd) %>%
  reframe(final_size_samp = mvrnorm(samp_size, final_size, projection_sd), 
          draw_id = 1:samp_size, .by = c("scenario_id", "location_id", "model_id")) %>%
  mutate(utility = linear_utility(final_size_samp, k_val, c_val)) %>%
  dcast(location_id + model_id + draw_id ~ scenario_id, value.var = "utility")

set.seed(1928)
projected_final_sizes = sim_out$model_sims %>%
  filter(scenario_id %in% c("S1", "S2")) %>%
  left_join(model_scenario_sd) %>%
  reframe(final_size_samp = mvrnorm(samp_size, final_size, projection_sd), 
          draw_id = 1:samp_size, .by = c("scenario_id", "location_id", "model_id")) %>%
  dcast(location_id + model_id + draw_id ~ scenario_id, value.var = "final_size_samp")

ggplot(data = projected_final_sizes %>% filter(location_id == 1, model_id %in% c("M5", "M10")), # for now
       aes(x = S1, y = S2, color = as.factor(model_id))) +
  geom_abline() + 
  geom_point(size = 2, alpha = 0.2) + 
  geom_point(data = true_final_sizes %>% filter(location_id == 1), color = "black", alpha = 0.4, size = 8, stroke = NA) +
  geom_point(data = true_final_sizes %>% filter(location_id == 1), color = "black") +
  stat_ellipse() +
  theme_bw()

prior_point = c(0.482, 0.134)
prior_r = 0.035
prior_col = "lightblue3"

model_point = c(0.505, 0.144)
model_r = 0.02
model_col = "blueviolet"

text_size = 2.5

ggplot(data = true_final_sizes %>% filter(location_id == 1), aes(x = S1, y = S2)) + 
  geom_point(aes(x = prior_point[1], y = prior_point[2]), color = prior_col, size = 2) +
  geom_ellipse(aes(x0 = prior_point[1], y0 = prior_point[2], a = prior_r, b = 0.04, angle = 12), fill = prior_col, alpha = 0.2, color = NA) +
  # geom_circle(aes(x0 = prior_point[1], y0 = prior_point[2], r = prior_r), fill = prior_col, alpha = 0.2, color = NA) +
  geom_point(aes(x = model_point[1], y = model_point[2]), color = model_col, size = 2) +
  geom_ellipse(aes(x0 = model_point[1], y0 = model_point[2], a = model_r, b = model_r*0.85, angle = -45), fill = model_col, alpha = 0.2, color = NA) + 
  geom_point(data = true_final_sizes %>% filter(location_id == 1), color = "black", alpha = 0.3, size = 13, stroke = NA) +
  geom_point(color = "black", size = 2) +
  geom_text(data = true_final_sizes %>% filter(location_id == 1), aes(label = "truth\n"), vjust = 0.2, size = text_size) +
  geom_text(aes(x = model_point[1], y = model_point[2], label = "model\n"), color = model_col, vjust = 0.2, size = text_size) +
  geom_text(aes(x = prior_point[1], y = prior_point[2], label = "prior\n"), color = prior_col, vjust = 0.2, size = text_size) +
  # now add decision lines
  geom_abline(intercept = -0.3, linetype = "dashed") +
  geom_abline(intercept = -0.74, slope = 1.8, linetype = "dotted") +
  scale_x_continuous(name = "projected value (action 1)") + 
  scale_y_continuous(name = "projected value (action 2)") + 
  theme_bw(base_size = 7) + 
  theme(#axis.text = element_blank(), 
        #axis.ticks = element_blank(), 
        panel.grid = element_blank()) 
ggsave("R/decision-theory/figures/projection_space.pdf", width = 3, heigh = 3)  


