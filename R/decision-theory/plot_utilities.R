library(dplyr)
library(reshape2)
library(ggplot2)
library(MASS) 
library(ggforce)
library(tidyr)
library(cowplot)

sim_out = readRDS("output/simulation/sim_out.rds")

source("R/simulation/0-parameters.R")
source("R/simulation/0-helper-functions.R")

linear_utility = function(y, j, c){
  return(y*j + c) 
}

quadratic_utility = function(y, k, c){
  return(k*y^2 + c)
}

j_val = c(5.5,6)*1e3 #(population in millions)
k_val = c(4.5,5)*1e3 #(population in millions)
c_val = c(2, 4)*1e3

prior_col = "lightblue3"
model_col = "blueviolet"

samp_size = 500

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
  mutate(projection_sd = runif(n(), 0, 0.0025))

chosen_location = 1

sim_out$model_sims %>%
  filter(scenario_id %in% c("S1", "S2"), location_id == chosen_location) %>%
  left_join(model_scenario_sd) %>%
  reframe(final_size_samp = mvrnorm(samp_size, final_size, projection_sd), 
          draw_id = 1:samp_size, .by = c("scenario_id", "location_id", "model_id")) %>%
  dcast(location_id + model_id + draw_id ~ scenario_id, value.var = "final_size_samp") %>% 
  ggplot(aes(x = S1, y = S2)) + 
  geom_point() + 
  stat_ellipse() + 
  facet_wrap(vars(model_id))

chosen_models = c("M5", "M10")

model_sims = sim_out$model_sims %>%
  filter(scenario_id %in% c("S1", "S2"), location_id == chosen_location, model_id %in% chosen_models) %>%
  left_join(model_scenario_sd) %>%
  reframe(final_size_samp = mvrnorm(samp_size, final_size, projection_sd), 
          draw_id = 1:samp_size, .by = c("scenario_id", "location_id", "model_id"))

projected_utilities = model_sims %>%
  mutate(utility = linear_utility(final_size_samp, 
                                  ifelse(scenario_id == "S1", j_val[1], j_val[2]), 
                                  ifelse(scenario_id == "S1", c_val[1], c_val[2]))) %>%
  dcast(location_id + model_id + draw_id ~ scenario_id, value.var = "utility")

projected_utilities_quadratic = model_sims %>%
  mutate(utility = quadratic_utility(final_size_samp, 
                                  ifelse(scenario_id == "S1", k_val[1], k_val[2]), 
                                  ifelse(scenario_id == "S1", c_val[1], c_val[2]))) %>%
  dcast(location_id + model_id + draw_id ~ scenario_id, value.var = "utility")

projected_final_sizes = model_sims %>%
  dcast(location_id + model_id + draw_id ~ scenario_id, value.var = "final_size_samp")

l_lin = range(projected_utilities %>% filter(location_id == 1, model_id %in% c("M5", "M10")) %>% dplyr::select(S1, S2))
l_quad = range(projected_utilities_quadratic %>% filter(location_id == 1, model_id %in% c("M5", "M10")) %>% dplyr::select(S1, S2))
l_util = range(c(l_lin, l_quad))
l_out = range(projected_final_sizes %>% filter(location_id == 1, model_id %in% c("M5", "M10")) %>% dplyr::select(S1, S2))
decision_lines = data.frame(y2 = seq(0, 0.25, length.out = 100))
decision_lines = decision_lines %>% 
  mutate(linear = 1/j_val[1] * (j_val[2] * y2 + c_val[2] - c_val[1]), 
         quadratic = sqrt(1/k_val[1] * (k_val[2] * y2 + c_val[2] - c_val[1]))) %>% 
  melt(c("y2"), variable.name = "utility_fn", value.name = "y1")

p1 = ggplot(data = projected_utilities, #
       aes(x = S1, y = S2, fill = as.factor(model_id))) +
  # geom_point(size = 2, alpha = 0.2) + 
  stat_ellipse(geom = "polygon", alpha = 0.3) +
  geom_abline() + 
  geom_point(data = projected_utilities %>% summarize(S1 = mean(S1), S2 = mean(S2), .by = c("model_id")), 
             aes(color = model_id), size = 0.8) + 
  labs(subtitle = "linear utility") +
  scale_color_manual(values = c(prior_col, model_col)) + 
  scale_fill_manual(values = c(prior_col, model_col)) + 
  scale_x_continuous(limits = l_util, name = "utility of action 1") +
  scale_y_continuous(limits = l_util, name = "utility of action 2") +
  theme_bw(base_size = 7) + 
  theme(legend.position = "none",
        panel.grid = element_blank())
p2 = ggplot(data = projected_utilities_quadratic, 
            aes(x = S1, y = S2, fill = as.factor(model_id))) +
  geom_abline(linetype = "dashed") + 
  geom_point(data = projected_utilities_quadratic %>% 
               summarize(S1 = mean(S1), S2 = mean(S2), .by = c("model_id")), 
             aes(color = model_id), size = 0.8) + 
  stat_ellipse(geom = "polygon", alpha = 0.3) +
  labs(subtitle = "quadratic utility") + 
  scale_color_manual(values = c(prior_col, model_col)) +
  scale_fill_manual(values = c(prior_col, model_col)) + 
  scale_x_continuous(expand = c(0,0), limits = l_util, name = "utility of action 1") +
  scale_y_continuous(expand = c(0,0), limits = l_util, name = "utility of action 2") +
  theme_bw(base_size = 7) + 
  theme(legend.position = "none",
        panel.grid = element_blank())
p3 = ggplot(data = projected_final_sizes,
       aes(x = S1, y = S2)) +
  # geom_point(size = 2, alpha = 0.2) + 
  stat_ellipse(aes(fill = as.factor(model_id)), geom = "polygon", alpha = 0.3) +
  geom_point(data = projected_final_sizes %>% 
               summarize(S1 = mean(S1), S2 = mean(S2), .by = c("model_id")), 
             aes(color = model_id), size = 0.8) + 
  geom_point(data = true_final_sizes %>% filter(location_id == 1), color = "black", alpha = 0.3, size = 8, stroke = NA) +
  geom_point(data = true_final_sizes %>% filter(location_id == 1), color = "black", size = 0.8) +
  geom_line(data = decision_lines, aes(x = y1, y = y2, linetype = utility_fn), color = "black") + 
  guides(color = FALSE) +
  scale_color_manual(values = c(prior_col, model_col), labels = c("prior", "model")) + 
  scale_fill_manual(values = c(prior_col, model_col), labels = c("prior", "model")) + 
  scale_linetype_discrete(labels = paste0(c("linear", "quadratic"), " utility")) + 
  scale_x_continuous(expand = c(0,0), name = "projected outcome of action 1") +
  scale_y_continuous(expand = c(0,0), name = "projected outcome of action 2") +
  theme_bw(base_size = 7) + 
  theme(legend.title = element_blank(), 
        legend.position = "bottom",
        panel.grid = element_blank())
l = cowplot::get_legend(p3)
plot_grid(
  plot_grid(plot_grid(p1, p2, ncol = 1), p3 + theme(legend.position = "noen"), 
            rel_widths = c(0.4, 0.6), align = "v", axis = "tb"), 
  l, rel_heights = c(0.9, 0.1), ncol = 1
)

ggsave("R/decision-theory/figures/projection_space.pdf", width = 6, heigh = 4)  

